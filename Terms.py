# includes the functions required for TERM-associated analysis
# the format of residue is like 'A10', no comma

from General import *
from PDB import *


class Term:
    def __init__(self, parent = None, seed = None, contact = []):
        # it is ok to initialize empty term
        self.parent = parent # the original PDB file
        self.seed = seed
        self.contact = contact # a list if not none

        self.residues = []
        self.frag = None # here is the prody atom selection
        self.pdbf = None
        # create a name for associated PDB file, based on the name of seed and contact

    # the second method of initialization; from an existing PDB file
    # only PDB files with certain names can be read
    def readFromPDB(self, pdbf):
        pdbfname = removePath(pdbf)
        pdbinfo = pdbfname.split('_')[1:-1]
        assert len(pdbinfo) > 0, 'the input file cannot be used to create a TERM object, quit...'
        self.seed = pdbinfo[0]
        self.contact = pdbinfo[1:]
        self.pdbf = pdbf
        conres = ConRes(pdbf)
        self.residues = [r.getChid() + str(r.getResnum()) for r in conres]

    # only Term object can call this method
    def _adjacentWindow(self, atoms, cid, resnum, flank):
        assert isinstance(resnum, int)
        assert isinstance(flank, int)
        resList = []
        for ires in atoms.iterResidues():
            if (ires.getIcode() != '') or (ires.getChid() != cid):
                continue
            iresnum = ires.getResnum()
            if abs(resnum - iresnum)  <= flank:
                resList.append(str(iresnum))
        return resList

    def makeFragment(self, flank = 1, dry = False):
        assert self.seed != None, 'the seed of TERM has not been defined...'
        assert isinstance(self.contact, list), 'the contact list of TERM has not been initialized...'
        assert self.pdbf == None, 'the pdb file for this term already exists as ' + self.pdbf
        assert os.path.isfile(self.parent), 'the template pdb file does not exist...'
        self.pdbf = '_'.join([getBase(removePath(self.parent))] + [self.seed] + self.contact) + '.pdb'

        # name output pdf file systematically

        atoms = parsePDB(self.parent, model = 1)
        frag = False
        for r in [self.seed] + self.contact:
            resList = self._adjacentWindow(atoms, r[0], int(r[1:]), flank)
            if not resList:
                print("fragments can't be made for %s in %s because its residue list is empty" % (self.seed, self.pdbf))
                return -1
            self.residues.extend([r[0] + x for x in resList])
            selectstr = 'chain ' + r[0] + ' and resnum ' + ' '.join(resList)
            selection = atoms.select(selectstr)
            if not frag:
                frag = selection
            else:
                frag = frag | selection # important: using | instead of + will autosort
        frag = frag.copy()
        self.residues = list(set(self.residues))
        self.residues = sorted(self.residues, key = lambda x : (x[0], int(x[1:])))
        if dry:
            return self.residues 
        else:
            writePDB(self.pdbf, frag)
        return 1

    def getResidues(self):
        return self.residues

    def findResidue(self, cid, resnum): # return the index of a residue in the TERM
        if cid + str(resnum) in self.residues:
            return self.residues.index(cid + str(resnum))
        else:
            return -1

    def getSeed(self):
        if self.seed == None:
            print('the seed of this term is not declared...')
            return -1
        else:
            return self.seed

    def getContacts(self):
        if self.contact == None:
            print('the contact of this term is not declared...')
            return -1
        else:
            return self.contact

    def getSegLen(self):
        assert self.pdbf != None
        conres = ConRes(self.pdbf)
        sort_res = sorted(conres,  key = lambda x: (x.getChid(), x.getResnum()))
        segment = [1]
        for i in range(1, len(sort_res)):
            if (sort_res[i].getChid() == sort_res[i-1].getChid()) and (sort_res[i].getResnum() -1 == sort_res[i-1].getResnum()):
                segment[-1] += 1
            else:
                segment.append(1)
        return segment


def contactList(profile, resnum, cid = '', outFile = None, dmin = 0.01, dmax = 1.0001, monomer = True):
    ''' <profile> a .profile of contacts
    <cid> chain id, can be empty for null chain name
    <resnum> residue number
    <outFile> writing output to this file. If None, then no output file
    <dmin, dmax> cutoffs for the contacts to consider
    <monomer> if true, only consider contacts from the same chain
    return: a list which has the chain id and residue number of all the contacting residues
    '''
    assert isinstance(dmin, float)
    assert isinstance(dmax, float)
    assert os.path.isfile(profile), 'the profile doesnt exist!'

    seed = ','.join(map(str, [cid, resnum]))
    cons = [] # list of contacts
    conress = [] # list of amino acids for contacts
    if outFile != None:
        ofh = open(outFile, 'w')
    N = 1 # counter of contacts
    with open(profile) as mapf:
        for line in mapf:
            if not re.match('contact', line):
                continue
            line = line.strip()
            if line.find('\t' + seed+'\t') >= 0:
                larr = line.split()
                res1, res2 = larr[1], larr[2]
                ind = [res1, res2].index(seed) # ind is 0 or 1
                contact, cenres, conres = larr[2-ind], larr[4+ind], larr[5-ind]
                cond = float(larr[3])
                if (cond <= dmin) or (cond > dmax):
                    continue
                if monomer and (not contact.startswith(cid)):
                    continue
                cons.append(contact.replace(',', ''))
                conress.append(conres)
                if outFile != None:
                    items = [str(x).replace(',', '') for x in [N, seed, contact, cond, cenres, conres]]
                    writeItems(items, ofh)
                N += 1
    if outFile != None:
        ofh.close()
    return cons, conress


# def replaceBfactor (pdbf, outf, dataf, res_col = 2, score_col = -1):
#     '''<pdbf> pdb file
#     <outf> output file, with B-factor filled with other metrics
#     <file> a file with numerical data by residues
#     <res_col> the column with residue information
#     <score_col> the column with scores to fill in
#     <cid> if False, only need to match residue number
#     '''
#     with open(pdbf) as pf, open(dataf) as df, open(outf, 'w') as of:
#
#         dat = {}
#         for dl in df:
#             dlspl = dl.split()
#             resinfo = dlspl[res_col]
#             dat[resinfo] = dlspl[score_col].strip()
#             break
#
#         for pl in pf:
#             pl = pl.rstrip('\n')
#             if pl[0:4] != 'ATOM':
#                 continue
#             # find cid and resnum
#             plspl = pl.split()
#             cid, resnum = pl[21].strip(), pl[22:27].strip()
#
#             if (cid + resnum in dat) == False:
#                 bfact = ''.rjust(6)
#             else:
#                 bfact = ('%.3f' % float(dat[cidc + resnum])).rjust(6)
#
#             # b-factor is from 61-66
#             if len(pl) > 60:
#                 left = pl[0:60]
#                 right = pl[66:]
#                 of.write(left + bfact + right + '\n')
#             if len(pl) <= 60:
#                 left = pl.ljust(60)
#                 of.write(left + bfact + '\n')

def adhocCentralResidue(dirname, head, rlist):
    nhits = []
    for r in rlist:
        seqf = head+'_'+dirname+'_'+r+'.seq'
        nhit = 0
        with open(seqf) as sf:
            for sl in sf:
                nhit+=1
        nhits.append(nhit)
    maxi = nhits.index(max(nhits))
    return rlist[maxi]

## functions in TERMANAL
## to standardize input, require Prody residue object as input for residues

# def _neighborInFrag(res1, res2, head, path):
#     '''test if res1 is the neighbor of res2. Definition of neighbor is as in Structure paper,
#     the central residue in a fragment is the neighbor of other residue included in the fragment.
#     '''
#     res1id = getResid(res1).strip()
#     res2c, res2n = res2.getChid(), str(res2.getResnum()) + res2.getIcode()
#     # looking for the .pdb file centred at res1
#     pdbf = path + '/' + head + '_' + res1id + '.pdb'
#     r = None
#     r = getResByInd(pdbf, res2c, res2n)
#     if r == None:
#         return 0
#     else:
#         return 1

#
# def neighborList(pdbf, res, path):
#     '''return all neighbors for res in pdbf, include res itself.
#     Will call function(_neighborInFrag)
#     '''
#     allres = ConRes(pdbf)
#     head = getBase( removePath(pdbf) )
#     path = absPath(path)
#     nb = []
#     for r in allres:
#         res2id = getResid(r)
#         if _neighborInFrag(res, res2, head, path) == 1:
#             nb.append(res2)
#     return nb


# def AsNeighborList(allres, pdbf, path):
#     '''return a dictionary describing how many time a residue receives the score from other residue in a structure (act as neighbor)
#     This is useful in averaging the raw score to get the smoothed score
#     '''
#     nn = {}
#     head = removePath( getBase(pdbf) )
#     for r in allres:
#         resid = getResid(r)
#         nn[resid] = 0
#     path = absPath(path)
#     for r in allres:
#         rid = getResid(r).strip()
#         # make a pseudofragment (not actually a pdb file)
#         listf = path + '/' + head + '_' + rid + '.list'
#         # read seeds from the list file
#         cons = readConsFromList(listf)
#         # declare the term
#         term = Term(parent = pdbf, seed = r.getChid() + str(r.getResnum()), contact = cons)
#         pseudofrag = term.makeFragments(pdbf, seeds, flank = 2, dry = True)
#         for res in pseudofrag:
#             nn[res] += 1
#     return nn
#
# def readConsFromList(listf):


# some old functions to calcultate GDT (needs identical atom numbers)

## gdt
# structure manipulation is from prody
# algorithm is adapted from clusco, the result should be identical to that of clusco
# Atoms: atoms of the structure to rotate, is a numpy array
# rAtoms: atoms of the reference structure, is a numpy array
# dcut: the cutoff for selecting the largest set of default is 3.5
# wSize: window size for initial segment, default is 4
# gdtcut: 
# return: a gdt of these two sets of atoms
# def gdtTransformation(Atoms, rAtoms, dcut = 4, wSize = 3, gdtcut = [1.0, 2.0, 4.0, 8.0], iter = 5):
#     assert rAtoms.shape[1] == 3
#     assert Atoms.shape == rAtoms.shape
#     protlen = rAtoms.shape[0]
#     maxAlign = 0
#     bestTran = None
#     bestGDT = 0
#     bestgdt = []
#     for ii in range(protlen-wSize+1):
#         flag = [False for x in range(protlen)]
#         fraglen = wSize
#         for n in range(iter):
#             for j in range(ii, ii+wSize):
#                 flag[j] = True
#             # get the superposition between residues which are true
#             sub_Atoms = Atoms[ [x for x in range(len(flag)) if flag[x] == True] ]
#             sub_rAtoms = rAtoms[ [x for x in range(len(flag)) if flag[x] == True] ]
#             trans = calcTransformation(sub_Atoms, sub_rAtoms) #prody
#             t_Atoms = applyTransformation(trans, Atoms) #prody
#             # now calculated atom distance after superposition
#             for iii in range(protlen):
#                 dist = spatialDistance(t_Atoms[iii], rAtoms[iii])
#                 if dist < dcut:
#                     flag[iii] = True
#                 else:
#                     flag[iii] = False
#             newfraglen = len([x for x in range(len(flag)) if flag[x] == True])
#             if newfraglen == fraglen:
#                 break
#             fraglen = newfraglen
#         [tmpGDT, tmpgdt] = calculateGDT(t_Atoms, rAtoms, gdtcut)
#         if tmpGDT > bestGDT:
#             bestGDT = tmpGDT
#             bestgdt = tmpgdt
#             bestTran = trans
#
#     ft_Atoms = applyTransformation(bestTran, Atoms) #prody # final transformation
#     # the contribution of each atom to gdt
#     ctb = contributionGDT(ft_Atoms, rAtoms, gdtcut)
#     return [bestGDT, bestgdt, ctb]
#
#
# def calculateGDT(Atoms, rAtoms, gdtcut):
#     assert Atoms.shape == rAtoms.shape
#     protlen = Atoms.shape[0]
#     gdt_d = [0 for x in gdtcut]
#     gdt = []
#     for i in range(protlen):
#         for ii in range(len(gdtcut)):
#             if spatialDistance(Atoms[i], rAtoms[i]) < gdtcut[ii]:
#                 gdt_d[ii] += 1
#     for j in range(len(gdtcut)):
#         gdt.append( gdt_d[j]/float(protlen) )
#     GDT = mean(gdt)
#     return [GDT, gdt_d]
#
# def contributionGDT(atoms, ratoms, gdtcut):
#     assert atoms.shape == ratoms.shape
#     result = []
#     protlen = atoms.shape[0]
#     for i in range(protlen):
#         byres = []
#         for ii in range(len(gdtcut)):
#             if spatialDistance(atoms[i], ratoms[i]) < gdtcut[ii]:
#                 byres.append(1)
#             else:
#                 byres.append(0)
#         result.append(byres)
#     return result
#
# def spatialDistance(a, b):
#     assert len(a) == len(b)
#     d = 0
#     for i in range(len(a)):
#         d += (a[i]-b[i])**2
#     d = math.sqrt(d)
#     return d
