import os, sys, shelve
import General, PDB

if len(sys.argv) - 1 != 2:
        print '<usage> [path to MASTER database list] [output FASTA file]'
        exit(0)

lst, fasta_db = sys.argv[1:]

out = open(fasta_db, 'w')
tmppdbf = '/tmp/tmp.%d.pdb' % os.getpid()
for l in open(lst):
        pdsf = l.strip()
        name = General.getBase(General.removePath(pdsf))
        os.system(General.PATH_master + "/parsePDS --pds " + pdsf + " --pdb " + tmppdbf)
        seqs = PDB.pdb2seq(tmppdbf)
        for c in seqs:
                out.write('>' + name+'_' + ('' if (len(seqs) == 1) else c) +'\n')
                out.write(seqs[c]+'\n')
out.close()
os.remove(tmppdbf)
