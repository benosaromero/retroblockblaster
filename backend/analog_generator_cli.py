import json
import sys
from pathlib import Path

from bksltk.analogs import make_analogs_combinations, write_numbered_parent_png
from PIL import ImageFont


_original_truetype = ImageFont.truetype


def safe_truetype(font, size=10, index=0, encoding="", layout_engine=None):
    try:
        return _original_truetype(font, size=size, index=index, encoding=encoding, layout_engine=layout_engine)
    except Exception:
        return ImageFont.load_default()


ImageFont.truetype = safe_truetype

def read_json(path):
    with open(path, "r") as handle:
        return json.load(handle)


def read_csv_smiles(path):
    smiles = []

    if not Path(path).exists():
        return smiles

    with open(path, "r") as handle:
        header = handle.readline().strip().split(",")

        smiles_index = 0
        for possible in ["smiles", "SMILES", "analog_smiles", "canonical_smiles"]:
            if possible in header:
                smiles_index = header.index(possible)
                break

        for line in handle:
            parts = line.strip().split(",")

            if len(parts) > smiles_index and parts[smiles_index]:
                smiles.append(parts[smiles_index])

    return smiles


def main():
    parent_smiles = sys.argv[1]
    modifications_json = sys.argv[2]
    output_prefix = sys.argv[3]

    modifications_dict = read_json(modifications_json)
    modifications_dict = {int(key): value for key, value in modifications_dict.items()}

    numbered_png = "{}_parent_numbered.png".format(output_prefix)
    write_numbered_parent_png(parent_smiles, numbered_png)

    make_analogs_combinations(
        modifications_dict,
        parent_smiles,
        [1, 2],
        output_prefix,
    )

    analogs = read_csv_smiles("{}.csv".format(output_prefix))

    print(json.dumps({
        "parent_smiles": parent_smiles,
        "numbered_png": numbered_png,
        "csv": "{}.csv".format(output_prefix),
        "png": "{}.png".format(output_prefix),
        "analogs": analogs,
    }))


if __name__ == "__main__":
    main()
