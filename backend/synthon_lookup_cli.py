import json
import sys

from server_old_lookup import lookup_synthon_building_blocks


def main():
    synthon_id = sys.argv[1]
    reference_smiles = sys.argv[2] if len(sys.argv) > 2 else ""

    result = lookup_synthon_building_blocks(
        synthon_id,
        reference_smiles=reference_smiles,
    )

    print(json.dumps(result))


if __name__ == "__main__":
    main()
