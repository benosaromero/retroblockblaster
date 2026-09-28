def lookup_synthon_building_blocks(synthon_id, limit=50, reference_smiles=""):
    import json
    import arthor

    synthons_path = "/nfs/exj/Fe/scripts/common_db/export/Synthons.json"
    db_path = "/mnt/nfs/exj/arthor_db_backup/public_arthor/bb_instock.smi"

    with open(synthons_path, "r") as f:
        synthons = json.load(f)

    synthon = None

    if isinstance(synthons, dict):
        synthon = synthons.get(str(synthon_id))

        if synthon is None:
            synthon = synthons.get(int(synthon_id)) if str(synthon_id).isdigit() else None

    elif isinstance(synthons, list):
        for item in synthons:
            if isinstance(item, dict) and str(item.get("id")) == str(synthon_id):
                synthon = item
                break

    if synthon is None:
        return {
            "error": "synthon not found",
            "synthon_id": str(synthon_id),
            "results": [],
        }

    smarts_list = synthon.get("inSMARTS") or synthon.get("smarts") or synthon.get("SMARTS") or synthon.get("query")
    smarts = smarts_list[0] if isinstance(smarts_list, list) else smarts_list

    if not smarts:
        return {
            "error": "synthon has no SMARTS query",
            "synthon_id": str(synthon_id),
            "synthon_name": synthon.get("Name", ""),
            "in_smarts": smarts_list,
            "results": [],
        }

    db = arthor.SubDb("{}.atdb".format(db_path), smifile=db_path)
    qry = arthor.Query(smarts, "Sma:QI")

    ref_fp = None

    if reference_smiles:
        try:
            from rdkit import Chem, DataStructs
            from rdkit.Chem import AllChem

            ref_mol = Chem.MolFromSmiles(reference_smiles)

            if ref_mol is not None:
                ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=2048)
        except Exception:
            ref_fp = None

    results = []
    for hit in db.search(qry):
        line = hit.decode("utf-8") if hasattr(hit, "decode") else str(hit)
        parts = line.strip().split()

        smiles = parts[0] if len(parts) > 0 else ""
        bb_id = parts[1] if len(parts) > 1 else ""

        tc_score = ""

        if ref_fp is not None and smiles:
            try:
                hit_mol = Chem.MolFromSmiles(smiles)

                if hit_mol is not None:
                    hit_fp = AllChem.GetMorganFingerprintAsBitVect(hit_mol, 2, nBits=2048)
                    tc_score = round(DataStructs.TanimotoSimilarity(ref_fp, hit_fp), 3)
            except Exception:
                tc_score = ""

        if ref_fp is not None:
            if tc_score == "" or tc_score < 0.2:
                continue

        results.append({
            "smiles": smiles,
            "bb_id": bb_id,
            "synthon_id": str(synthon_id),
            "synthon_name": synthon.get("Name", ""),
            "db": "bb-stock",
            "tc_score": tc_score,
        })

        if len(results) >= limit:
            break

    if ref_fp is not None:
        results.sort(key=lambda item: item.get("tc_score", 0), reverse=True)

