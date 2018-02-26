# these are the basic supporting standard packages. 
# other files should start with "from General import *"
import os, sys, re, glob, random, time, shutil, argparse, math
import subprocess as sub
import numpy as np

# some path manipulation


def getBase(filename):
    part = filename.rpartition('.')
    if part[0] == '':
        return part[-1]
    return part[0]


def getExt(filename):
    part = filename.rpartition('.')
    return part[-1]


def removePath(filepath, keep=1):
    '''<keep> specifies the number of subdirectories to keep on the right'''
    if filepath[-1] == '/':
        filepath = filepath.rstrip('/')
    part = filepath.rsplit('/', keep)
    if len(part) <= keep:
        return filepath
    else:
        return '/'.join(part[1:])


def getPath(filepath):
    if filepath[-1] == '/':
        filepath = filepath.rstrip('/')
    part = filepath.rpartition('/')
    if part[0] == '':
        return '.'
    return part[0]


def absPath(path):
    '''full path name'''
    return os.path.realpath(path)


def selfbin(path):
    '''the directory of the current file'''
    return os.path.dirname(os.path.realpath(path))


def changeExt(filename, ext):
    return getBase(filename) + '.' + ext


def listDir(path):
    return sorted([d for d in os.listdir(path) if os.path.isdir(path + '/' + d)])


def file2array(file, strip=True):
    f = open(file)
    arr = f.readlines()
    if strip is True:
        arr = [l.rstrip('\n') for l in arr]
    return arr


def skipTo(string, fh, begins=False):
    '''<str> a string
    <fh> a opened file handle
    <begins> only find matches in the beginning of line
    returns the hit line and the filehandle'''
    for line in fh:
        if begins and re.match(string, line):
            return line, fh
        elif re.search(string, line):
            return line, fh
    return -1


def cleanFiles(string):
    '''use it when rm doesn't work (sometimes too many arguments)
    <string> used as in glob'''
    files = glob.glob(string)
    for f in files:
        if os.path.isfile(f) is True:
            os.remove(f)
        elif os.path.isdir(f) is True:
            shutil.rmtree(f)


def writeItems(alist, ofh, sep='\t'):
    string = sep.join(alist) + '\n'
    ofh.write(string)


def carefulImport(mod):
    import imp
    found = False
    try:
        imp.find_module(mod)
        found = True
    except ImportError:
        print('Warning:' + mod + ' is not an existing module')
    return found

### Interface ###
# provide interface to other programs at steady paths, such as MASTER, MD, confind, Rosetta, etc.

PATH_blast = '/home/grigoryanlab/library/blast/blast-2.2.26/bin'
PATH_confind = '/home/anthill/fzheng/home/confind/confind'
PATH_master = '/home/grigoryanlab/library/MASTER/bin'
PATH_rotLib = '/home/grigoryanlab/library/confind/rotlibs/DB-2010'
PATH_termaster = ''
PATH_usearch = '/home/anthill/fzheng/home/software/usearch8.0' # usearch is a software to cluster sequence based on sequence identity and remove redundancy


PATH_namd = '/home/grigoryanlab/library/bin/namd2'
PATH_charmm = '/home/grigoryanlab/library/bin/charmrun'
PATH_vmd = '/home/grigoryanlab/local/bin/vmd'


PATH_rosetta = '/home/anthill/fzheng/home/rosetta2014'
PATH_pyrosetta = '/home/anthill/fzheng/home/software/PyRosetta.ScientificLinux-r56316.64Bit'
PATH_fpddemo = '/home/anthill/fzheng/home/software/FlexPepDock_AbInitio'


# thesis directory
PATH_thesisData = '/home/anthill/fzheng/home/Thesis/Data'

