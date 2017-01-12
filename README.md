# pdbrotate
Quick tool to rotate atoms in pdb files


## align.py
Aligns the molecule according to its principial axis.

    python3 align.py -f peptide.pdb -o out.pdb


## rotate.py
Rotates molecule around axis.

    python3 rotate.py -f peptide.pdb -o out.pdb --angle=90 --axis=0,0,1
