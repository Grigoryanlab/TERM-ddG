from General import *
from PDB import *
from Master import *

# from one line in the .match file, parse the PDB_chain id
def pdbFromMatch(line):
    name = removePath(line).split('.')[0]
    return name

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
