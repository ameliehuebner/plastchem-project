# rdkit tutorial steal
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from IPython.display import SVG

mfp_generator = AllChem.GetMorganGenerator(radius=3)

for i in range(1,20):
    m1 = Chem.MolFromSmiles(df_plast_chem["Properties"]["canonical_smiles"][i])
    fp1 = mfp_generator.GetSparseCountFingerprint(m1)

    ao = AllChem.AdditionalOutput()

    ao.CollectBitInfoMap()

    fp = mfp_generator.GetFingerprint(m1,additionalOutput=ao)
    bi = ao.GetBitInfoMap()
    # mfp2_svg = Draw.DrawMorganBit(m1, list(bi.keys())[1], bi, useSVG=True)
    # drawer = Draw.rdMolDraw2D.MolDraw2DSVG(450, 150)
    # #draw the molecule
    # drawer.DrawMolecule(m1)
    # drawer.FinishDrawing()
    # # get the SVG string
    # svg = drawer.GetDrawingText()
    # # fix the svg string and display it
    # display(SVG(svg.replace('svg:','')))

#DataStructs.DiceSimilarity(fp1,fp2)