from General import *

def confind(mbin = PATH_confind, **kwargs): # only works in python3.3
    cmd  = [mbin]
    cmd.append('--pp')

    for key, value in kwargs.iteritems():
        cmd.extend(['--'+key, str(value)])
    sub.call(cmd)


def createPDS(mbin = PATH_master, dry = False, **kwargs):
    silent = False
    cmd = [mbin + '/createPDS']

    for key, value in kwargs.iteritems():
        cmd.extend(['--'+key, str(value)])

    if dry == True:
        return ' '.join(cmd)
    else:
        if (silent):
            devnull = open(os.devnull, 'w')
            sub.call(cmd, stdout=devnull)
            devnull.close()
        else:
            sub.call(cmd)
        

def masterSearch (mbin = PATH_master, dry = False, rmsdcut = 2.0, bbrmsd = True, **kwargs):
    '''most arguments is the same with master program; if dry is True, have a dry run, only return the command'''
        
    silent = False
    cmd = [mbin + '/master']
        
    cmd.extend(['--rmsdCut', str(rmsdcut)])
    if bbrmsd:
        cmd.append('--bbRMSD')

    for key, value in kwargs.iteritems():
        cmd.extend(['--'+key, str(value)])

    if dry == True:
        return ' '.join(cmd)
    else:
        if (silent):
            devnull = open(os.devnull, 'w')
            sub.call(cmd, stdout=devnull)
            devnull.close()
        else:
            sub.call(cmd)


def matchInFile(mbin = PATH_master, dry = False, otype = 'match', **kwargs):
    cmd = [mbin + '/master']
    if not 'outType' in kwargs.keys():
        cmd.extend(['--outType', otype])
    for key, value in kwargs.iteritems():
        cmd.extend(['--'+key, str(value)])
    if dry == True:
        return ' '.join(cmd)
    else:
        sub.call(cmd)

    
def extractPDB(mbin = PATH_master, dry = False, **kwargs):
    cmd = [mbin + '/extractPDB']
    for key, value in kwargs.iteritems():
        cmd.extend(['--'+key, str(value)])
    if dry == True:
        return ' '.join(cmd)
    else:
        sub.call(cmd)
