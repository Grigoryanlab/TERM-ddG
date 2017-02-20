__author__ = 'fanzheng'

import sys
sys.path.insert(1, '/home/ifsdata/grigoryanlab/library/FanPythonMods/') # needed to find General if FanPythonMods is not already in the user's path
from General import *
import PDB

par = argparse.ArgumentParser()
par.add_argument('--p', nargs= '+', help = 'input structures')
par.add_argument('--o', help = 'a file of mutation as output')
args = par.parse_args()

outf = open(args.o, 'w')

for p in args.p:
    resi = PDB.ConRes(p)
    renamep = getPath(p) + '/' + removePath(p).replace('_', '-')
    if renamep != p:
        os.system('cp ' + p  + ' ' + renamep)
    for r in resi:
        chain, pos, wt = r.getChid(), r.getResnum(), r.getResname()
        outstr = [removePath(getBase(renamep)), chain, PDB.t2s(wt), str(pos), 'A', '0']
        outstr = '\t'.join(outstr)
        outf.write(outstr + '\n')

outf.close()
