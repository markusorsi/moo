import numpy as np
from sklearn.decomposition import PCA
from rdkit import Chem 

xLOS_smarts = {
    'HYD' : ['[C]', '[a]', '[S]', '[F]', '[Cl]', '[Br]', '[I]'],
    'HBA' : ['[!$([#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]'],
    'HBD' : ['[!$([#6,H0,-,-2,-3])]']
}

xLOS_labels = xLOS_smarts.keys()
pca = PCA(n_components=3)

def get_propmat(mol):
        propmat = np.zeros((len(xLOS_labels), mol.GetNumAtoms()))
        for i, label in enumerate(xLOS_labels):
            category = xLOS_smarts[label]
            for smarts in category:
                pattern = Chem.MolFromSmarts(smarts)
                matched = False
                for match in mol.GetSubstructMatches(pattern, uniquify=True):
                    for j in match:
                        propmat[i, j] = 1
                    matched = True
                if matched: break
        
        return propmat
    

def get_crossed_propmat(seed, query):
    crossed_propmat = np.zeros([len(xLOS_labels), seed.GetNumAtoms(), query.GetNumAtoms()])
    seed_propmat = get_propmat(seed)
    query_propmat = get_propmat(query)
    for i in range(len(xLOS_labels)):
        crossed_propmat[i] = np.outer(seed_propmat[i], query_propmat[i])
    return crossed_propmat


def get_distmat(seed_coord, query_coord):
    distmat = np.zeros([len(seed_coord), len(query_coord)])
    for i, coords1 in enumerate(seed_coord):
        for j, coords2 in enumerate(query_coord):
            distmat[i, j] = np.linalg.norm(coords1-coords2)
    return distmat


def calc_xlos(seed, query):

    seed_conf = seed.GetConformer()
    query_conf = query.GetConformer()

    seed_coord = seed_conf.GetPositions()
    query_coord = query_conf.GetPositions()
        
    propmat = get_crossed_propmat(seed, query)
    distmat = get_distmat(seed_coord, query_coord)

    seed_propmat = get_propmat(seed)
    query_propmat = get_propmat(query)

    xlos_score = np.zeros(3)

    for i in range(len(xlos_score)):
        if np.count_nonzero(seed_propmat[i])+np.count_nonzero(query_propmat[i]) == 0:
            xlos_score[i] = 0
        else:
            xlos_score[i] = np.sum(np.exp(-np.square(distmat))*propmat[i])/(np.count_nonzero(seed_propmat[i])+np.count_nonzero(query_propmat[i]))

    return np.sum(xlos_score)