import sys

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


def main():
    smiles = sys.argv[1] if len(sys.argv) > 1 else ""
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print("")
        return

    drawer = rdMolDraw2D.MolDraw2DSVG(320, 220)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    print(drawer.GetDrawingText())


if __name__ == "__main__":
    main()
