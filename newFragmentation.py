__author__ = 'fanzheng'

from General import *
import PDB, Master, conGraph, Terms, Stability, Analyze, Cluster
import pickle

SB = selfbin(sys.argv[0])
par = argparse.ArgumentParser()
par.add_argument('--i', required=True, help='input PDB file')
par.add_argument('--homof', help='the file with homologous information')
par.add_argument('--db', default = '/home/anthill/fzheng/ironfs/searchDB/bc-30-sc-20141022-newpds', help = 'the directory of the searching database')
par.add_argument('--dbl', help = 'list of the searching database')
par.add_argument('--he', help = 'the header of the output files')
par.add_argument('--c1', default = 20000, type = int, help = 'cutoff for top N matches')
par.add_argument('--c2', default = 100, type = int, help = 'the criteria of increasing the order of sub-TERMs')
par.add_argument('--c3', default = 2, type = int, help = 'the maximum order of sub-TERMs to consider')
par.add_argument('--nr', default=0.4, type=float, help='the level of redundancy removal')
args = par.parse_args()

# save input argument
argdict = vars(args)
os.system('rm -f *.pkl')
pklname = str(os.getpid()) + '.shared.pkl'
pklpath = os.path.realpath(pklname)
pkl = open(pklname, 'w')
pickle.dump(argdict, pkl)
pkl.close()

# get a contact graph
confind_out = changeExt(args.i, 'conf')
if not os.path.isfile(confind_out):
    Master.confind(p=args.i, o=confind_out, rLib=PATH_rotLib)

G = conGraph.conGraph(confind_out)
edges = G.edges()
edges = sorted(edges, key = lambda x:(x[0][0], int(x[0][1:]), x[1][0], int(x[1][1:])))

residues = PDB.ConRes(args.i)
selflist = 'residue.txt'
pairlist = 'contact.txt'
with open(selflist, 'w') as sf:
    for r in residues:
    	residue_id = r.getChid() + str(r.getResnum())
       	sf.write(residue_id + '\n')
with open(pairlist, 'w') as pf:
    for e in edges:
    	residue1, residue2 = e[0], e[1]
       	pf.write(residue1 + ' ' + residue2 + '\n')

homos = []
if args.homof != None:
    homos =  Analyze.findHomo(args.homof, type =2)

Jobs = []
dirname = removePath(getBase(args.i))

selfts = []
# search all local fragments
for n in G.nodes():
	if n in target_residues:
	    t = Terms.Term(parent=args.i, seed=n)
	    t.makeFragment(flank=2)
	    selfts.append(t)
	    crind = t.findResidue(n[0], n[1:])
	    rmsd_eff = Stability.rmsdEff(t.getSegLen())
	    cmd = ['python', SB + '/EnsemblePreparation.py',
	           '--pkl', pklpath,
	           '--p', t.pdbf,
			   '--homo', ' '.join(homos),
			   '--rmsd', str(rmsd_eff),
			   '--crind', str(crind),
	           '--ncon', '0']
	    cmd = ' '.join(cmd)
	    jobid = n
	    job = Cluster.jobOnCluster([cmd], jobid, args.he + '_' + changeExt(t.pdbf, 'seq'))
	    Jobs.append(job)
	    job.submit(3)
	    time.sleep(0.5)

Cluster.waitJobs(Jobs, giveup_time=1, subdir=False)

# search all pair fragments
pairts = []
for e in G.edges():
    r1, r2 = e[0], e[1]
    t = Terms.Term(parent=args.i, seed=r1, contact=[r2])
    t.makeFragment(flank=1)
    pairts.append(t)
    cenr = Terms.adhocCentralResidue(dirname, args.he, [r1, r2])
    crind =  t.findResidue(cenr[0], cenr[1:]) # should consider two positions in redundancy removal; but worry about that later
    rmsd_eff = Stability.rmsdEff(t.getSegLen())
    cmd = ['python', SB + '/EnsemblePreparation.py',
           '--pkl', pklpath,
           '--p', t.pdbf,
		   '--homo', ' '.join(homos),
		   '--rmsd', str(rmsd_eff),
		   '--crind', str(crind),
           '--ncon', '1']
    cmd = ' '.join(cmd)
    jobid = r1 + '_' + r2
    job = Cluster.jobOnCluster([cmd], jobid, args.he + '_' + changeExt(t.pdbf, 'seq'))
    Jobs.append(job)
    job.submit(3)
    time.sleep(0.5)

Cluster.waitJobs(Jobs, giveup_time=1, subdir=False)

exit()
# for all possible triple fragments, the criteria of search is at least two sub-pair fragments have more than 100 hits
tripts = []
for n in G.nodes():
    contactsofn = G.neighbors(n)
    contactsofn = sorted(contactsofn, key = lambda x: (x[0], int(x[1:])))
    for j in range(len(contactsofn)-1):
        for k in range(j+1, len(contactsofn)):
            nn1 = contactsofn[j]
            nn2 = contactsofn[k]
            seqf1 = '_'.join([args.he, dirname, n, nn1]) + '.seq'
            seqf2 = '_'.join([args.he, dirname, n, nn2]) + '.seq'
        if (not os.path.isfile(seqf1)) or (not os.path.isfile(seqf2)):
            continue
        else:
            nhit1, nhit12 = 0, 0
            with open(seqf1) as sm1:
                nhit1 = 0
                for lm in sm1:
                    nhit1 +=1
            with open(seqf2) as sm2:
                nhit2 = 0
                for lm in sm2:
                    nhit2 +=1
            if min(nhit1, nhit2) >= args.c2:
                t = Terms.Term(parent=args.i, seed=n, contact=[nn1, nn2])
                t.makeFragment(flank=1)
                tripts.append(t)
                cenr = Terms.adhocCentralResidue(dirname, args.he, [n,nn1,nn2])
                crind =  t.findResidue(cenr[0], cenr[1:]) # should consider two positions in redundancy removal; but worry about that later
                rmsd_eff = Stability.rmsdEff(t.getSegLen())
                cmd = ['python', SB + '/EnsemblePreparation.py',
                       '--pkl', pklpath,
                       '--p', t.pdbf,
                       '--homo', ' '.join(homos),
                       '--rmsd', str(rmsd_eff),
                       '--crind', str(crind),
                       '--ncon', '2']
                cmd = ' '.join(cmd)
                jobid = ' '.join([n, nn1, nn2])
                job = Cluster.jobOnCluster([cmd], jobid, args.he + '_' + changeExt(t.pdbf, 'seq'))
                Jobs.append(job)
                job.submit(3)
                time.sleep(0.5)

Cluster.waitJobs(Jobs, giveup_time=1, subdir=False)

# in this specific code does not consider higher order terms
