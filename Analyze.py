from General import *
from PDB import *
from Master import *

# calculate the contact degree between two residues
# def getConD (cmap, res1, res2):
#     '''get contact degree from .cmap files
#     if nothing is found, return 0
#     the format of residues: A,1'''
#     mapf = open(cmap)
#     for l in mapf:
#         if not re.match('contact', l):
#             continue
#         l = l.rstrip('\n')
#         if ( l.find( res1+'\t' ) >= 0 ) and ( l.find( res2+'\t' ) >=0 ):
#             arr = l.split()
#             return float(arr[3])
#     mapf.close()
#     return 0

# from one line in the .match file, parse the PDB_chain id
def pdbFromMatch(line):
    name = removePath(line).split('.')[0]
    return name

# read a .homo file and return a dictionary
# def readHomoDict(homo):
#     hdict = {}
#     lines = open(homo).readlines()
#     for l in lines:
#         items = l.strip().split()
#     hdict[items[0]] = items[1:]
#     return hdict

def removeHomo(matchf, homo, usedRows = None, dry = True, inhead = None, outhead = None, log = False):
    '''<matchf> a .match file
    <homo> a list of pdbid to avoid
    <usedRows> the row number to consider in the input .match files. If None, use all rows
    <dry> if dry is true, do not output new file, only return updated usedRows
    <inhead> the header of input file
    <outhead> the header of output file
    '''
    assert getExt(matchf) == 'match', 'a .match file should be provided! Error in ' + matchf
    seqf = changeExt(matchf, 'seq')
    assert os.path.isfile(matchf) and os.path.isfile(seqf), '.seq or .match files do not exist! Error in ' + matchf

    lines = open(matchf).readlines()
    if usedRows == None:
        usedRows = list(range(len(lines)))

    updatedRows = []
    n_removed = 0
    for r in usedRows:
        line = lines[r]
        pdbid = pdbFromMatch(line)
        if pdbid in homo:
            n_removed += 1
        else:
            updatedRows.append(r)

    if not dry:
        assert (inhead != None) and (outhead != None)
        omatchf = matchf.replace(inhead, outhead)
        oseqf = changeExt(omatchf, 'seq')
        with open(omatchf, 'w') as omf, open(oseqf, 'w') as osf:
            slines = open(seqf).readlines()
            for r in updatedRows:
                omf.write(lines[r])
                osf.write(slines[r])
        if log:
            logf = changeExt(omatchf, 'log')
            with open(logf, 'w') as log:
                log.write('removed ' + str(n_removed) + ' sequences')
    return updatedRows


# def trimByRMSD (inp, coln, rcut, usedRows = None, dry = True, inhead = None, outhead = None):
#     '''<inp> the input file, can be either .match file or .seq file
#     <coln> the column in the input file to look for rmsd values, from zero
#     <rcut> the rmsd cutoff
#     <usedRows> the row number to consider in the input .match files. If None, use all rows
#     <inhead> the header of the input file name, which is unique to procedure.
#         for example, the header for 'uniq_bbrmsd1_rmsd0.8_1A2P_AA32_1.match' is 'uniq_bbrmsd1_rmsd0.8'
#     <outhead> the header of the output file to replace the input header
#     '''
#     assert getExt(inp) != ('match' or 'seq'), 'input not correctly provided! Error in ' + inp
#     if getExt(inp) == 'match':
#         match = inp
#         seq = changeExt(inp,'seq')
#     else:
#         seq = inp
#         match = changeExt(inp, 'match')
#     # need both .match and .seq files
#     assert os.path.isfile(match) and os.path.isfile(seq), '.seq or .match files do not exist! Error in ' + inp
#
#     rmsds = readColumn(seq, coln)
#     for x in range(len(rmsds)):
#         if float(rmsds) > rcut:
#             cutrow = x
#             break
#     updatedRows = [i for i in usedRows if i < cutrow]
#
#     if not dry:
#         assert (inhead != None) and (outhead != None)
#         omatchf = match.replace(inhead, outhead)
#         oseqf = changeExt(omatchf, 'seq')
#         with open(omatchf, 'w') as omf, open(oseqf, 'w') as osf, open(match) as mf, open(seq) as sf:
#             mlines = mf.readlines()
#             slines = sf.readlines()
#             for r in updatedRows:
#                 omf.write(mlines[r])
#                 osf.write(slines[r])
#     return updatedRows


# def underRMSD(col, cut, sorted = True):
#     '''return the number of sequence under a certain rmsd in a .seq file
#     '''
#     i = -1
#     if sorted == False:
#         col.sort(key = float)
#     for i in range(len(col)):
#         if float(col[i]) > cut:
#             return i
#     return i+1 # if all sequence in the file are under the cutoff


# def rmsdOfnseq(col, n, sorted = True):
#     '''return the maximum rmsd in the top n sequences
#     also returns n. if n is larger than the length of sequence file, modify the value of n
#     '''
#     if sorted == False:
#         col.sort(key = float)
#     n = min(n, len(col))
#     return [col[n-1], n]


def readColumn (file, coln, top = -1, skiprow = 0):
    '''<file> a file to read
    <coln> a column number, start with 0
    <top> top number of lines to read 
    return: a list of column items
    '''
    assert os.path.isfile(file)
    assert  skiprow >=0
    lines = file2array(file)[skiprow:]
    col = []
    count = 0
    for line in lines:
        if len(line.split()) <= coln:
            col.append('NA')
        else:
            col.append(line.split()[coln])
        count += 1
        if top > 0:
            if count == top:
                break
    return col


def readMultiColumn (file, coln, top = -1, skiprow = 0):
    '''<file> a file to read
    <coln> a list of column numbers to read, start with 0
    <top> top number of lines to read
    return: a list of column items
    '''
    cols = []
    for n in coln:
        cols.append(readColumn(file, n, top, skiprow))
    return cols


# def readAADistribution(col):
#     counts = [0 for i in range(20)]
#     for r in col:
#         aaind = aaindex[t2s(r)]
#         if aaind > 0:
#             counts[aaind] += 1
#     return counts


# def informationContent(col, lowcount = True):
#     lencol = len(col)
#     H = {}
#     for item in col:
#         item = t2s(item)
#         if not item in H:
#             H[item] = 1
#         else:
#             H[item] += 1
#     I = math.log(20)/math.log(2)
#     for k in H.keys():
#         p = float(H[k])/lencol
#         I += p*math.log(p)/math.log(2)
#     if lowcount:
#         I -= 19.0/(2*lencol*math.log(2))
#         if I < 0:
#             I = 0
#     return I


# def informationContentQuick(col, norm = False):
#     I = math.log(20)/math.log(2)
#     col = map(float, col)
#     for p in col:
#         if not norm:
#             p /= sum(col)
#         I += p * math.log(p)/math.log(2)
#     return I


def index_from_match(line):
    '''return a list of residue indices of match regions in a given line from match files'''
    ilist = []
    idxrange = re.search('(\[\(.+\)\])', line).group(0)
    numbers = re.findall('\d+', idxrange)
    for i in range(0, len(numbers), 2):
        start, end = int(numbers[i]), int(numbers[i+1])
        for x in range(start, end+1):
            ilist.append(x)
    return ilist


# def determineBin(bins, val):
#     '''<bins> a list of numbers for bin borders, including the boundaries
#     <val> value to assign a bin index for'''
#     assert isinstance(val, float)
#     if (val < bins[0]) or (val > bins[-1]):  # for values outside the boundary, return 'NA'
#         return 'NA'
#     for i in range(1, len(bins)):
#         if val < bins[i]:
#             return i - 1
#     return len(bins) - 2


# parses a fasta file
# returns list of (identifier, seq) tuples
# def parseFasta(fastaFile):
#   lines = file2array(fastaFile)
#   chid = lines[0][1:]
#   seq = lines[1]
#   out = []
#   for line in lines[2:]:
#     if line[0] == '>': # name
#       out.append((chid, seq))
#       chid = line[1:]
#       seq = ''
#     else:
#       seq += line
#
#   out.append((chid, seq))
#   return out


# read a fasta file into a dictionary of pid -> sequence
# def fasttaDict(fastFile):
#   dic = {}
#   for seqname, seq in parseFasta(fastFile):
#       dic[seqname] = [aaindex[a] for a in seq]
#   return dic


# convert a sequence file with letters to a file in which sequences are represented as intergers
# def seq321(seqf, outf, start = 1):
#     lines = file2array(seqf)
#     with open(outf, 'w') as out:
#         for l in lines:
#             seq = ''.join([t2s(a) for a in l.split()[start:]])
#             out.write(seq + '\n')


def createHomoProfile(fastaf, outf):
    blastout = changeExt(fastaf, 'blastout')
    cmd = PATH_blast + '/blastpgp -d /home/grigoryanlab/home/jack/from-fan/Data/blast/db/pdbaa/pdbaa -i ' + fastaf + ' -b 0 -e 1.0 -v 10000 -o ' + blastout
    os.system(cmd)
    with open(blastout) as bo, open(outf, 'w') as hf:
        homologs =[]
        for bl in bo:
            if bl.startswith('Query'):
                if len(homologs) > 0:
                    hf.write(' '.join(homologs) + '\n')
                homologs = []
                homologs.append(bl.strip().split()[-1])
            elif bl.startswith('pdb'):
                pid, chain = bl.split()[0].split('|')[1:]
                pid = pid.lower()
                homologs.append(pid + '_' + chain)
            hf.write(' '.join(homologs) + '\n')


def findHomo(homof, type = 1):
    with open(homof) as hf:
        if type == 1:
            Homo = {}
            for hl in hf:
                items = hl.strip().split()
                Homo[items[0]] = items[1:]
        else:
            Homo = []
            for hl in hf:
                items = hl.strip().split()
                Homo.extend(items)
    return Homo
