import numpy as np
from rdkit import Chem 
from rdkit.Chem import rdMolAlign, AllChem
from tqdm import tqdm

def ConformersFromSmiles(smiles, nConfs=10, maxIters=500, verbose=False):
    """
    Generate conformers for a molecule specified by its SMILES string.

    Parameters:
    - smiles (str): SMILES representation of the molecule.
    - nConfs (int, optional): Number of conformers to generate. Default is 10.
    - maxIters (int, optional): Maximum number of force-field optimization iterations for each conformer. Default is 500.
    - verbose (bool, optional): If True, print information about conformer minimization convergence. Default is False.

    Returns:
    Chem.Mol: The molecule with generated conformers.
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = AllChem.AddHs(mol)

    confs = AllChem.EmbedMultipleConfs(mol, nConfs, randomSeed=42, 
                                        useRandomCoords=False,
                                        enforceChirality=True,
                                        useExpTorsionAnglePrefs=True, 
                                        useBasicKnowledge=True,
                                        useSmallRingTorsions=True, 
                                        useMacrocycleTorsions=True,  
                                        )
    converged = 0
    for conf in confs:
        opt = AllChem.MMFFOptimizeMolecule(mol, confId=conf, maxIters=maxIters)
        converged+=opt
        if (opt == -1) and verbose:
            print(f'Force field could not be setup for conformer {conf}!')
        else:
            converged += opt
    if verbose:
        print(f'{converged} conformer minimisations failed to converge')
    return mol


def ConfsToAlignedMolsList(multiConfMol):
    """
    Align conformers within a molecule and convert them to a list of aligned molecules.

    Parameters:
    - multiConfMol (Chem.Mol): A molecule with multiple conformers.

    Returns:
    List[Chem.Mol]: A list of aligned molecules derived from the conformers.
    """
    AllChem.AlignMolConformers(multiConfMol)
    mols = []
    cids = [x.GetId() for x in multiConfMol.GetConformers()]
    for cid in cids:
        mol = Chem.MolToMolBlock(multiConfMol, confId=cid)
        mol = Chem.MolFromMolBlock(mol, removeHs=False)
        mols.append(mol)
    return mols


def GetBestO3AScorePair(molList1, molList2, verbose=False):
    """
    Find the pair of molecules with the highest O3A score between two lists of molecules.

    Parameters:
    - molList1 (List[Chem.Mol]): List of molecules for comparison.
    - molList2 (List[Chem.Mol]): Another list of molecules for comparison.
    - verbose (bool, optional): If True, display progress using tqdm. Default is False.

    Returns:
    Tuple[Chem.Mol, Chem.Mol]: The pair of molecules with the highest O3A score.
    """
    scores = []
    coords = []
    for i in tqdm(range(0, len(molList1)), disable= not verbose):
        for j in range(i, len(molList2)):
            pyO3A = AllChem.GetO3A(molList1[i], molList2[j])
            scores.append(pyO3A.Score())
            coords.append((i, j))
    max_score = np.argmax(scores)
    idx1, idx2 = coords[max_score]
    
    return molList1[idx1], molList2[idx2]


def GetBestOverlapPairFromSmiles(smiles1, smiles2, nConfs=100, RemoveHs=True, verbose=False):
    """
    Generate conformers from SMILES strings, align and compare conformers, and return the best-matching pair.

    Parameters:
    - smiles1 (str): SMILES representation of the first molecule.
    - smiles2 (str): SMILES representation of the second molecule.
    - nConfs (int, optional): Number of conformers to generate. Default is 100.
    - RemoveHs (bool, optional): If True, remove hydrogens from the aligned molecules. Default is True.
    - verbose (bool, optional): If True, display progress and O3A score. Default is False.

    Returns:
    Tuple[Chem.Mol, Chem.Mol]: The best-matching pair of molecules after conformer generation, alignment, and comparison.
    """
    confs1 = ConformersFromSmiles(smiles1, nConfs)
    confs2 = ConformersFromSmiles(smiles2, nConfs)

    mols1 = ConfsToAlignedMolsList(confs1)
    mols2 = ConfsToAlignedMolsList(confs2)
    
    best1, best2 = GetBestO3AScorePair(mols1, mols2, verbose=verbose)

    pyO3A = rdMolAlign.GetO3A(best1, best2)
    rmsd, trans_matrix = pyO3A.Trans()
    AllChem.TransformMol(best1, trans_matrix)

    if verbose:
        print(f'O3A Score: {pyO3A.Score()}')

    if RemoveHs:
        best1 = AllChem.RemoveAllHs(best1)
        best2 = AllChem.RemoveAllHs(best2)

    return best1, best2


def GetBestOverlapPairFromMollist(mols1, mols2, RemoveHs=True, verbose=False):
    """
    Compare and align conformers from two lists of molecules and return the best-matching pair.

    Parameters:
    - mols1 (List[Chem.Mol]): List of molecules for comparison.
    - mols2 (List[Chem.Mol]): Another list of molecules for comparison.
    - RemoveHs (bool, optional): If True, remove hydrogens from the aligned molecules. Default is True.
    - verbose (bool, optional): If True, display progress and O3A score. Default is False.

    Returns:
    Tuple[Chem.Mol, Chem.Mol]: The best-matching pair of molecules after alignment and comparison.
    """
    best1, best2 = GetBestO3AScorePair(mols1, mols2, verbose=verbose)

    pyO3A = rdMolAlign.GetO3A(best1, best2)
    rmsd, trans_matrix = pyO3A.Trans()
    AllChem.TransformMol(best1, trans_matrix)

    if verbose:
        print(f'O3A Score: {pyO3A.Score()}')

    if RemoveHs:
        best1 = AllChem.RemoveAllHs(best1)
        best2 = AllChem.RemoveAllHs(best2)

    return best1, best2


def SaveBestOverlapAsPDB(smiles1, smiles2, nConfs=100, RemoveHs=True, name1='mol1', name2='mol2'):
    """
    Generate conformers, align and compare them, and save the best-matching pair as PDB files.

    Parameters:
    - smiles1 (str): SMILES representation of the first molecule.
    - smiles2 (str): SMILES representation of the second molecule.
    - nConfs (int, optional): Number of conformers to generate. Default is 100.
    - RemoveHs (bool, optional): If True, remove hydrogens from the aligned molecules. Default is True.
    - name1 (str, optional): Name for the first PDB file. Default is 'mol1'.
    - name2 (str, optional): Name for the second PDB file. Default is 'mol2'.

    Returns:
    None

    Saves:
    - PDB files for the best-matching pair with specified names.
    """
    confs1 = ConformersFromSmiles(smiles1, nConfs)
    confs2 = ConformersFromSmiles(smiles2, nConfs)

    mols1 = ConfsToAlignedMolsList(confs1)
    mols2 = ConfsToAlignedMolsList(confs2)

    best1, best2 = GetBestO3AScorePair(mols1, mols2)

    pyO3A = rdMolAlign.GetO3A(best1, best2)
    rmsd, trans_matrix = pyO3A.Trans()
    AllChem.TransformMol(best1, trans_matrix)

    print(f'O3A Score: {pyO3A.Score()}')

    if RemoveHs:
        best1 = AllChem.RemoveAllHs(best1)
        best2 = AllChem.RemoveAllHs(best2)
    
    Chem.MolToPDBFile(best1, f'{name1}.pdb')
    Chem.MolToPDBFile(best2, f'{name2}.pdb')