import os, sys, shelve
import General, PDB

if len(sys.argv) - 1 != 3:
        print '<usage> [a list of pdb files] [output .fasta file] [output shelve .db file]'
        exit(0)

lst, fasta, shelve_db = sys.argv[1:]

out = open(fasta, 'w')
for l in open(lst):
        pdbf = l.strip()
        name = General.removePath(pdbf)
        seqs = PDB.pdb2seq(pdbf)
        out.write('>'+pdbf+'\n')
        for c in seqs: # because only single chain
                out.write(seqs[c]+'\n')

db = shelve.open(shelve_db)

for line in open(fasta).readlines():
        if line.startswith('>'):
                key = line[1:7]
        else:
                db[key] = line.strip()
db.close()