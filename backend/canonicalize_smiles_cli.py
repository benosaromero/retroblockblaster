import sys
from rdkit import Chem


def main():
    smiles = sys.argv[1] if len(sys.argv) > 1 else ""
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print("")
        return

    print(Chem.MolToSmiles(mol, canonical=True))


if __name__ == "__main__":
    main()
