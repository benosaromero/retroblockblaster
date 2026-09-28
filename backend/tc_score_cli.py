import json
import sys

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


def fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def main():
    payload = json.loads(sys.stdin.read())

    reference_smiles = payload.get("reference_smiles", "")
    candidates = payload.get("candidates", [])

    ref_fp = fingerprint(reference_smiles)

    results = []

    for item in candidates:
        smiles = item.get("smiles", "")
        score = ""

        if ref_fp is not None:
            fp = fingerprint(smiles)

            if fp is not None:
                score = round(DataStructs.TanimotoSimilarity(ref_fp, fp), 3)

        item["tc_score"] = score
        results.append(item)

    print(json.dumps(results))


if __name__ == "__main__":
    main()
