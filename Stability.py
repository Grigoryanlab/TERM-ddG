# any functions specific for mustpress may be dumped here. 
# If a functon could be useful for multi-project, consider dump it into global function repository

from General import *
from PDB import *
from Analyze import *
from scipy.io import loadmat

class Mutation:
    def __init__(self, pid, cid, resn, wt, mt, expG):
        self.p = pid
        self.c = cid
        self.n = int(resn)
        self.w = wt
        self.m = mt
        self.e = expG
        self.dir = pid+'_'+cid+resn
    

def getMutation(line):
    order = ['pid', 'cid', 'wt', 'resn', 'mt', 'expG']
    arr = line.split()
    info = {}
    assert len(arr) == len(order), 'Line does not have the same length as order...'
    for i in range(len(arr)):
        info[order[i]] = arr[i]
    return Mutation(info['pid'], info['cid'], info['resn'], info['wt'], info['mt'], info['expG'])
            

# def getResidueSelfInfo(pid, cid, resn, loc, ext = '.clean.freedom'):
#     env = loc + pid.lower() + ext
#     with open(env) as fh:
#         for l in fh:
#             items = l.strip().split()
#             if items[0] != 'freedom':
#                 continue
#             if items[1] != cid + ',' + str(resn):
#                 continue
#             freedom, phi, psi = map(float, items[2:5])
#             return freedom, phi, psi
#     return 'NA', 'NA', 'NA'


def degFreedom(segLens, persLen):
    '''calculates the degrees of freedom for a fragment based on persistence length;
    <segLen>s is an array containing the length of each segment in an fragment - can be an int if there is only one segment;
    <persLen> is the persistence length;'''
    if not isinstance(segLens, list): # must be a one segment fragment
        segLens = [segLens]
    segLens = [float(segLen) for segLen in segLens]

    a = math.exp(-1.0/persLen) # persistence length between adjacent residues

    c = 0
    for segLen in segLens: # go through each segment
        d = a / (1.0 - a)
        c += d * (segLen - 1.0) - (d ** 2) * (1.0 - a ** (segLen - 1.0))
      
    fragLen = sum(segLens) # total number of residues in the fragment
    return fragLen * (1.0 - (2.0 / (fragLen * (fragLen - 1.0)))* c)

def rmsdEff(segLens, persLen=18, rmsdMax=1.0):
    '''calculates the effective rmsd for a fragment;
    <segLens> is an array containing the length of each segment in an fragment - can be an int if there is only one segment;
    <persLen> is the persistence length
    <rmsdMax> is the rmsd plateau at the limit (as the number of residues and/or segments increases)'''
    if not isinstance(segLens, list):
        segLens = [segLens]
    return round(rmsdMax / math.sqrt(sum(segLens) / degFreedom(segLens, persLen)), 3)


# create a .tab file for a user-defined range in a PDB file
# def scanPositions(pdbf, cid, start, end, outf):
#     mol = parsePDB(pdbf).select('chain ' + cid).copy()
#     pid = removePath(pdbf).split('.')[0]
#     ofh = open(outf, 'w')
#     for res in mol.iterResidues():
#         if (res.getResnum() >= start) and (res.getResnum() <= end):
#             outstr = '\t'.join([pid, cid, t2s(res.getResname()), str(res.getResnum()), 'A', '0']) + '\n'
#             ofh.write(outstr)


def removeRedundancy(matchf, head, crind, seqdb, nr=0.4, wd=15, usedRows = None):
    '''
    Return the indices in the original .seq files remained after redundancy removal
    :param matchf:
    :param head:
    :param crind:
    :param db:
    :param nr:
    :param wd:
    :param usedRows:
    :return: used_indices:
    '''
    # db = '/home/anthill/fzheng/home/searchDB/statistics/bc-30-sc-20141022.peprm2.db'
    # database = shelve.open(db)
    tempfile = changeExt(matchf, 'seqcontext.fasta')
    nr_tempfile =changeExt(matchf, 'nr.fasta')
    lines = open(matchf).readlines()
    if usedRows == None:
        usedRows = list(range(len(lines)))
    with open(tempfile, 'w') as tempfh:
        for r in usedRows:
            ml = lines[r]
            match_region_indices = index_from_match(ml)
            cindex = match_region_indices[crind]
            matchid = getBase(removePath(ml.split()[1]))
            fullseq = seqdb[matchid]
            if cindex - wd < 1:
                seqcontext = fullseq[0:(2 * wd + 1)]
            elif cindex + wd > len(fullseq):
                seqcontext = fullseq[-(2 * wd + 1):]
            else:
                seqcontext = fullseq[(cindex - wd - 1):(cindex + wd)]
            tempfh.write('>match:'+str(r)+'\n'+seqcontext+'\n')
    used_indices = []

    if len(usedRows) > 0:
        sub.call([PATH_usearch, '-cluster_fast', tempfile, '-id', str(nr), '-centroids', nr_tempfile, '-fulldp', '-query_cov', '1', '-target_cov', '1', '-maxgaps', '0'])
        with open(nr_tempfile) as nrf:
            for nrs in nrf:
                if nrs[0] == '>':
                    nrid = nrs.strip().split(':')[-1]
                    used_indices.append(int(nrid))
    return used_indices

def whetherIncludeHighOrder(pos, head, cons, minmatch=100):
    use = 1
    for i in range(len(cons)):
        cons_minus = cons[0:i] + cons[i+1:]
        minus = '_'.join(cons_minus)
        seq_minus = head + '_'+ pos + '_' + minus +'.seq'
        if not os.path.isfile(seq_minus):
            use = 0
        else:
            with open(seq_minus) as sm:
                nhit = 0
                for lm in sm:
                    nhit +=1
                if nhit < minmatch:
                    use = 0
    return use

# def parseStride(stridef, chain, resn, field, refarea = None):
#     '''Usage:
#     <stridef> the output of stride
#     <chain> chain id
#     <resn> residue number in the pdb file
#     <field> the field to parse. Must be one of the following: ss (2nd structure), phi, psi, area_abs, area_rel
#     <refarea> a dictionary of all the reference asa values
#     '''
#     assert field in ['ss', 'phi', 'psi', 'area_abs', 'area_rel']
#     for line in open(stridef):
#         if not line.startswith('ASG'):
#             continue
#         cid, num = line[9], line[11:15]
#         if (cid != chain) or (str(resn) != num.strip()):
#             continue
#         if field == 'ss':
#             return line[24]
#         if field == 'phi':
#             return line[42:49].strip()
#         if field == 'psi':
#             return line[52:59].strip()
#         if field == 'area_abs':
#             return line[61:69].strip()
#         if field == 'area_rel':
#             assert refarea != None
#             res = t2s(line[5:8])
#             return float(line[61:69].strip()) / refarea[res]
#     return -1


def loadParams(matf, var = 'optParams'):
    '''
    read the file generated by MATLAB as numpy file
    :param matf:
    :return:
    '''
    optparams = loadmat(matf)[var]
    return optparams


# def flatPositionScanning(file, col, resi, skiprow=1):
#     '''
#     processing the results of mutation scannning. In input file 20 mutations for one position are on the same row, and in the output file one mutation takes one row
#     :param file: input file
#     :param col: starting index for prediction values
#     :param resi: the column to find residue number
#     :param skiprow: skip the header in the input file
#     :return: a dictionary with format like {'1A': -0.03 ...}
#     '''
#     mutdict = {}
#     with open(file) as fh:
#         rowi = 0
#         for l in fh:
#             rowi += 1
#             if rowi <= skiprow:
#                 continue
#             info = l.strip().split()
#             assert len(info) >= col + 20, 'not enough elements in each row'
#             vals = info[col: col+20]
#             pos = info[resi]
#             for i in range(20):
#                 aa = aatypes[i]
#                 key = pos + aa
#                 mutdict[key] = vals[i]
#     return mutdict


# def parseMutScan(mutList, predict, startcol = 6):
#     '''
#
#     :param mutList:
#     :param predict:
#     :return:
#     '''
#     mlist =[]
#     out = []
#     with open(mutList) as mls:
#         for ml in mls:
#             mlist.append(ml.strip())
#     with open(predict) as pre:
#         for pl in pre:
#             info = pl.strip().split()
#             wt = info[2]
#             wtind = aaindex[wt]
#             pos = info[3]
#             scores = info[startcol:]
#             for i in range(20):
#                 if pos + aatypes[i] in mlist:
#                     out.append(float(scores[i]) - float(scores[wtind]))
#     return out


# def parseTermMaster(line):
#     line = line.strip().split()
#     termid = line[2]
#     segments = line[4:]
#     return termid, segments


# def TERMmatchSize(tid):
#     loc = '/home/ironfs/scratch/cmack2357/minCover/coverSL/20_1d0/terms/'
#     matchf = loc + tid[0:3] + '/' + tid + '.match'
#     n = 0
#     for l in open(matchf):
#         n += 1
#     return n
