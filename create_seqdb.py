import os, sys, shelve
import General, PDB

if len(sys.argv) - 1 != 2:
        print '<usage> [a list of pdb files] [output shelve .db file]'
        exit(0)

lst, shelve_db = sys.argv[1:]

# out = open(fasta, 'w')
db = shelve.open(shelve_db)
for l in open(lst):
        pdbf = l.strip()
        name = General.getBase(General.removePath(pdbf))
        seqs = PDB.pdb2seq(pdbf)
        for c in seqs: # because only single chain
                db[name] = seqs[c]
db.close()