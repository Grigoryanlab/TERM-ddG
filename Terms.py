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
