__author__ = 'fanzheng'

import sys
sys.path.insert(1, '/home/grigoryanlab/library/FanPythonMods/') # needed to find General if FanPythonMods is not already in the user's path
from General import *
import PDB
import scipy.io as sio

par  =argparse.ArgumentParser()
par.add_argument('--p', required=True, help='input structure')
par.add_argument('--m', required = True, help = 'the output of matlab')
par.add_argument('--i', required = True, help = 'a residue file')
par.add_argument('--ii', required = True, help = 'a contact file')
par.add_argument('--o', required=True, help= 'output file')
args = par.parse_args()

vals = sio.loadmat(args.m)['currentParamsValues']

residueNames = {}
residues = PDB.ConRes(args.p)
for r in residues:
    rid = r.getChid() + str(r.getResnum())
    rname = r.getResname()
    residueNames[rid] = rname

tscore = 0
with open(args.o, 'w') as of, open(args.i) as rf, open(args.ii) as cf:
    rfl = rf.readlines()
    cfl = cf.readlines()
    for i in range(len(rfl)):
        rid = rfl[i].strip()
        rname = residueNames[rid]
        rind = PDB.aaindex[PDB.t2s(rname)]
        try:
            score = vals[i][0].flatten()[rind]
        except:
            pass
        tscore += score
        of.write('\t'.join([format(score, '.3f'), rid, rname, 'NA', 'NA']) + '\n')
    for j in range(len(rfl), len(rfl) + len(cfl)):
        cid1, cid2 = cfl[j-len(rfl)].strip().split()
        cname1, cname2 = residueNames[cid1], residueNames[cid2]
        cind1, cind2 = PDB.aaindex[PDB.t2s(cname1)], PDB.aaindex[PDB.t2s(cname2)]
        score = vals[j][0].reshape(20,20)[cind2, cind1] # note the order, the "main" residue should be put behind
        of.write('\t'.join([format(score, '.3f'), cid1, cname1, cid2, cname2]) + '\n')
        tscore += score
print(tscore)
