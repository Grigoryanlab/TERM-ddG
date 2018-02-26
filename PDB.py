# from prody import *
from General import *
found = carefulImport('prody')
if found:
    from prody import *
    confProDy(verbosity='critical')

a2aaa = {
'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 
'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 
'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 
'Y': 'TYR', 'V': 'VAL', 'X': 'XXX'
}

aaa2a = {aaa:a for a, aaa in a2aaa.iteritems()}
# add unnatural amino acids (the most common ones)
aaa2a['ASX'] = aaa2a['ASN']
aaa2a['CSO'] = aaa2a['CYS']
aaa2a['GLX'] = aaa2a['GLU'] # or GLN
aaa2a['HIP'] = aaa2a['HIS']
aaa2a['HSC'] = aaa2a['HIS']
aaa2a['HSD'] = aaa2a['HIS']
aaa2a['HSE'] = aaa2a['HIS']
aaa2a['HSP'] = aaa2a['HIS']
aaa2a['MSE'] = aaa2a['MET']
aaa2a['SEC'] = aaa2a['CYS']
aaa2a['SEP'] = aaa2a['SER']
aaa2a['TPO'] = aaa2a['THR']
aaa2a['PTR'] = aaa2a['TYR']
aaa2a['XLE'] = aaa2a['LEU'] # or ILE

natAA = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
         'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'MSE']

d3_l3 = {
'DAL':'ALA', 'DAR':'ARG', 'DAN':'ASN', 'DAS':'ASP', 'DCS':'CYS', 'DGN':'GLN', 'DGU':'GLU', 'DHI':'HIS', 
'DIL':'ILE', 'DLE':'LEU', 'DLY':'LYS', 'DME':'MET', 'DPH':'PHE', 'DPR':'PRO', 'DSE':'SER', 'DTH':'THR',
'DTR':'TRP', 'DTY':'TYR', 'DVA':'VAL'
} # three letter D-AA to the name of natural AA
d3_d4 = { d3: 'D'+l3 for d3, l3 in d3_l3.iteritems()} # 4-letter D-AA to the name of natural AA

aatypes = 'A C D E F G H I K L M N P Q R S T V W Y'
aatypes = aatypes.split()
aaindex = {aatypes[x] : x for x in range(20)}
aaindex['X'] = -1

aaFreq = [0.0795, 0.0133, 0.0589, 0.0684, 0.0410, 0.0692, 0.0234, 0.0582, 0.0584, 0.0955, 0.0217, 0.0433, 0.0451, 0.0382, 0.0520, 0.0602, 0.0542, 0.0699, 0.0139, 0.0357]

aaMW = [71, 103, 115, 129, 147, 57, 137, 113, 128, 113, 131, 114, 97, 128, 156, 87, 101, 99, 186, 163]


def s2t (res):
    assert len(res) == 1, 'Invalid input, should be one letter...'
    assert res in a2aaa, 'Invalid input, no such letter...'
    return a2aaa[res]


def t2s (res):
    assert len(res) == 3, 'Invalid input, should be three letter...'
    if not res in aaa2a: # this may not be an error, but an unrecognize unnatural AA
        return 'X'
    return aaa2a[res]


def getResid(res):
    '''<res> an Prody object'''
    return res.getChid() + str(res.getResnum()) + res.getIcode()

def ConRes(pdbf):
    '''return a list of all residues'''
    mol = parsePDB(pdbf, model=1)
    residues = []
    for res in mol.iterResidues():
        if res.getResname() in aaa2a:
            residues.append(res)
    return residues


def ConResDict(pdbf):
    '''return a dictionary of all residues'''
    mol = parsePDB(pdbf, model=1)
    residues = {}
    for res in mol.iterResidues():
        if res.getResname() in aaa2a:
            resid = getResid(res)
            residues[resid] = res
    return residues

    
def getResByInd(pdbf, cid, resn, mode = 1):
    '''mode 1 uses resnum, mode 2 uses iresnum (the i-th residue in the input structure); 
    if not find, will return None'''
    atoms = parsePDB(pdbf, model=1)
    res = None # if not find, will return None
    for ires in atoms.iterResidues():
        if mode == 1:
            if (str(ires.getResnum()) + ires.getIcode() == resn) and (ires.getChid() == cid):
                res = ires
                break
        if mode == 2:
            resn = int(resn)
            if (ires.getResindex() == resn-1) and (mode == 2):
                res = ires
                break
    return res


def printAllResid(pdbf):
    ''' print all residue names to a file '''
    atoms = parsePDB(pdbf, model=1)
    out = changeExt(pdbf, 'res')
    with open(out, 'w') as of:
        for ires in atoms.iterResidues():
            outstr = str(ires.getChid()) + ',' + str(ires.getResnum()) + ires.getIcode()
            of.write(outstr + '\n')


def findPositionInPDB(pdbf, resn, cid = None):
    '''return the index of given residue (iresnum) in the given pdb fragment (start from 1)
    <cid> If not provided, look for the first residue with the correct residue number.
    return -1 if not find
    '''
    atoms = parsePDB(pdbf, model=1)
    i = 0
    for ires in atoms.iterResidues():
        if str(ires.getResnum()) + ires.getIcode() == str(resn):
            if cid == None:
                return i+1
            elif ires.getChid() == cid:
                return i+1
        i += 1
    return -1 # if not find, return -1


def pdb2seq(pdbf, seqres = False):
    '''give a pdb file, return sequences as a dictionary of chain : sequence '''
    lines = open(pdbf).readlines()
    seqs = {}
    if seqres:
        for l in lines:
            if not re.match('SEQRES', l):
                continue
            arr = l.rstrip('\n').split()
            (cid, ress) = (arr[2], arr[4:])
            if not cid in seqs: # initialize a new chain
                seqs[cid] = ''
            for resname in ress:
                seqs[cid] += t2s(res)
    else:
        mol = parsePDB(pdbf, model=1)
        # protein = mol.select('protein').copy() # this is dangerous; some residue may not get selected; known by MASTER, but not known by Prody, could have index problem
        # but in this way, need to make sure the processed PDB is clean (water, etc.); processed by MASTER. 
        for res in mol.iterResidues():
            cid, resname = res.getChid(), res.getResname()
            if not cid in seqs:
                seqs[cid] = ''
            seqs[cid] += t2s(resname)
    return seqs
  

def renumber(pdbf, opdbf):
    p = parsePDB(pdbf).select('protein').copy()
    ind = 0
    lastCh = ''
    for r in p.iterResidues():
        Chid = r.getChid()
        if Chid != lastCh:
            ind = 0
            lastCh = Chid
        ind += 1
        r.setResnum(ind)
        r.setIcode('')
    writePDB(opdbf, p)