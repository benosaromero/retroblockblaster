import csv
import json
import os
import sys
import uuid

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


REACTIONS_JSON = os.environ.get(
    "REACTIONS_JSON",
    "/nfs/exj/Fe/scripts/common_db/export/Reactions.json"
)

SYNTHONS_JSON = os.environ.get(
    "SYNTHONS_JSON",
    "/nfs/exj/Fe/scripts/common_db/export/Synthons.json"
)


class BaseMolecule:
    def __init__(self, smiles, cid, mol):
        self.smi = smiles
        self.cid = cid
        self.mol = mol


class Reaction:
    def __init__(self, name, rxn_id, synthons, breakdown_smarts,
                 breakdown_rxn, build_smarts, build_rxn, num_reactants):
        self.name = name
        self.id = rxn_id
        self.synthons = synthons
        self.breakdown_smarts = breakdown_smarts
        self.breakdown_rxn = breakdown_rxn
        self.build_smarts = build_smarts
        self.build_rxn = build_rxn
        self.num_reactants = num_reactants


def extract_rxns(rxn_ids, rxn_map):
    rxns = []

    for rxn_id in rxn_ids:
        rxn_data = rxn_map["Reactions"][rxn_id]

        name = rxn_data["Name"]
        synthon_ids = rxn_data["SynthonIDs"].split(",")
        breakdown_smarts = rxn_data["Breakdown_SMARTS"]
        build_smarts = rxn_data["Build_SMARTS"]
        num_reactants = rxn_data["NumReactants"]

        try:
            breakdown_rxn = AllChem.ReactionFromSmarts(breakdown_smarts)
            build_rxn = AllChem.ReactionFromSmarts(build_smarts)
        except Exception:
            continue

        rxns.append(
            Reaction(
                name,
                rxn_id,
                synthon_ids,
                breakdown_smarts,
                breakdown_rxn,
                build_smarts,
                build_rxn,
                num_reactants,
            )
        )

    return rxns


def check_bbs_hac_min(bb_mols):
    return True


def run_analoging_assess(base_mol, rxns):
    analog_dict = {}
    utilized_rxns = []

    for rxn in rxns:
        num_reactants = int(rxn.num_reactants)
        try:
            products = rxn.breakdown_rxn.RunReactants((base_mol.mol,))
        except Exception:
            continue

        try:
            bbs = products[0]
        except IndexError:
            continue

        if len(bbs) != num_reactants:
            continue

        if not check_bbs_hac_min(bbs):
            continue

        synthon_ids = rxn.synthons

        if rxn.id == "CC-10.41":
            synthon_ids = list(reversed(synthon_ids))

        bb_mols = []
        for i in range(num_reactants):
            bb_mols.append([bbs[i], synthon_ids[i]])

        analog_dict[rxn.id] = bb_mols
        utilized_rxns.append(rxn)

    return analog_dict, utilized_rxns

def safe_filename(text):
    keep = []
    for char in text:
        if char.isalnum() or char in ("-", "_"):
            keep.append(char)
        else:
            keep.append("_")
    return "".join(keep)


def write_molecule_svg(smiles, image_dir, filename, width, height):
    if not smiles:
        return ""

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""

    if not os.path.exists(image_dir):
        os.makedirs(image_dir)

    safe_name = safe_filename(filename)
    path = os.path.join(image_dir, safe_name + ".svg")

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    svg = drawer.GetDrawingText()

    with open(path, "w") as f:
        f.write(svg)

    return path

def find_retro_reactions(smiles, cid, reaction_id=None):
    with open(REACTIONS_JSON) as json_file:
        rxn_map = json.load(json_file)

    with open(SYNTHONS_JSON) as json_file:
        synthon_map = json.load(json_file)

    if reaction_id:
        rxn_ids = [reaction_id]
    else:
        rxn_ids = rxn_map["Reaction_IDs"]

    rxns = extract_rxns(rxn_ids, rxn_map)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES: {}".format(smiles))

    base_mol = BaseMolecule(smiles, cid, mol)

    analog_dict, utilized_rxns = run_analoging_assess(base_mol, rxns)

    rows = []

    for rxn in utilized_rxns:
        bb_list = analog_dict[rxn.id]
        pair_id = str(uuid.uuid4())
        reagent_count = len(bb_list)

        bb1_smiles = ""
        bb1_synthon_id = ""
        bb1_synthon_name = ""

        bb2_smiles = ""
        bb2_synthon_id = ""
        bb2_synthon_name = ""

        bb3_smiles = ""
        bb3_synthon_id = ""
        bb3_synthon_name = ""

        bb4_smiles = ""
        bb4_synthon_id = ""
        bb4_synthon_name = ""

        if reagent_count >= 1:
            bb1_smiles = Chem.MolToSmiles(bb_list[0][0])
            bb1_synthon_id = bb_list[0][1]
            bb1_synthon_name = synthon_map.get(str(bb1_synthon_id), {}).get("Name", "")

        if reagent_count >= 2:
            bb2_smiles = Chem.MolToSmiles(bb_list[1][0])
            bb2_synthon_id = bb_list[1][1]
            bb2_synthon_name = synthon_map.get(str(bb2_synthon_id), {}).get("Name", "")

        if reagent_count >= 3:
            bb3_smiles = Chem.MolToSmiles(bb_list[2][0])
            bb3_synthon_id = bb_list[2][1]
            bb3_synthon_name = synthon_map.get(str(bb3_synthon_id), {}).get("Name", "")

        if reagent_count >= 4:
            bb4_smiles = Chem.MolToSmiles(bb_list[3][0])
            bb4_synthon_id = bb_list[3][1]
            bb4_synthon_name = synthon_map.get(str(bb4_synthon_id), {}).get("Name", "")

        rows.append([
            smiles,
            rxn.id,
            rxn.name,
            rxn.build_smarts,
            rxn.breakdown_smarts,
            pair_id,
            reagent_count,
            bb1_smiles,
            bb1_synthon_id,
            bb1_synthon_name,
            bb2_smiles,
            bb2_synthon_id,
            bb2_synthon_name,
            bb3_smiles,
            bb3_synthon_id,
            bb3_synthon_name,
            bb4_smiles,
            bb4_synthon_id,
            bb4_synthon_name,
        ])

    rows.sort(key=lambda row: row[4])

    return rows

def read_input(input_path):
    with open(input_path) as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split()
                smiles = parts[0]
                cid = parts[1] if len(parts) > 1 else "input_molecule"
                return smiles, cid

    raise ValueError("Input file is empty")


def main():
    input_path = sys.argv[1]
    output_path = sys.argv[2]

    job_dir = os.path.dirname(output_path)
    image_dir = os.path.join(job_dir, "images")

    smiles, cid = read_input(input_path)
    rows = find_retro_reactions(smiles, cid)

    with open(output_path, "w") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "original_smiles",
            "original_svg",
            "reaction_id",
            "reaction_name",
            "build_smarts",
            "breakdown_smarts",
            "pair_id",
            "reagent_count",
            "bb1_smiles",
            "bb1_svg",
            "bb1_synthon_id",
            "bb1_synthon_name",
            "bb2_smiles",
            "bb2_svg",
            "bb2_synthon_id",
            "bb2_synthon_name",
            "bb3_smiles",
            "bb3_svg",
            "bb3_synthon_id",
            "bb3_synthon_name",
            "bb4_smiles",
            "bb4_svg",
            "bb4_synthon_id",
            "bb4_synthon_name",
        ])

        for row in rows:
            original_smiles = row[0]
            reaction_id = row[1]
            reaction_name = row[2]
            build_smarts = row[3]
            breakdown_smarts = row[4]
            pair_id = row[5]
            reagent_count = row[6]

            bb1_smiles = row[7]
            bb1_synthon_id = row[8]
            bb1_synthon_name = row[9]

            bb2_smiles = row[10]
            bb2_synthon_id = row[11]
            bb2_synthon_name = row[12]

            bb3_smiles = row[13]
            bb3_synthon_id = row[14]
            bb3_synthon_name = row[15]

            bb4_smiles = row[16]
            bb4_synthon_id = row[17]
            bb4_synthon_name = row[18]

            original_svg = write_molecule_svg(
                original_smiles, image_dir, "original", 500, 350
            )
            bb1_svg = write_molecule_svg(
                bb1_smiles, image_dir, "bb1_" + pair_id, 240, 180
            )
            bb2_svg = write_molecule_svg(
                bb2_smiles, image_dir, "bb2_" + pair_id, 240, 180
            )
            bb3_svg = write_molecule_svg(
                bb3_smiles, image_dir, "bb3_" + pair_id, 240, 180
            )
            bb4_svg = write_molecule_svg(
                bb4_smiles, image_dir, "bb4_" + pair_id, 240, 180
            )

            writer.writerow([
                original_smiles,
                original_svg,
                reaction_id,
                reaction_name,
                build_smarts,
                breakdown_smarts,
                pair_id,
                reagent_count,
                bb1_smiles,
                bb1_svg,
                bb1_synthon_id,
                bb1_synthon_name,
                bb2_smiles,
                bb2_svg,
                bb2_synthon_id,
                bb2_synthon_name,
                bb3_smiles,
                bb3_svg,
                bb3_synthon_id,
                bb3_synthon_name,
                bb4_smiles,
                bb4_svg,
                bb4_synthon_id,
                bb4_synthon_name,
            ])

        if not rows:
            out.write("# No matching retro reactions found\n")


if __name__ == "__main__":
    main()
