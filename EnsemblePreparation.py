import sys
from General import *
import Master, Stability, Cluster, Analyze
import pickle, shelve #  for sharing params

from socket import gethostname
print('Node: %s' % gethostname())

SB = selfbin(sys.argv[0])
par = argparse.ArgumentParser()
# inherit the parameter space from mutationListIteration.py
par.add_argument('--pkl', required = True, help = 'a shared pickle file, to pass arguments from parent script')
# make sure the argument do not clash between the two name space
par.add_argument('--p', nargs = '+', required = True, help = 'term structure')
par.add_argument('--homo', nargs = '*', help = 'a list of PDB-id to avoid as homologs')
par.add_argument('--rmsd', nargs = '+', required = True, help = 'rmsd cutoff')
par.add_argument('--crind', nargs = '+', required=True, help='the index of the seed residue in the ensemble')
par.add_argument('--ncon', required=True, help='the maximum number of contacts in this job, only used to create the finish label')
args = par.parse_args()

odir = os.getcwd()
ldir = Cluster.createLocalSpace()

os.chdir(ldir)
os.system('cp ' + args.pkl + ' ./')
pklname = removePath(args.pkl)
pkl = open(pklname)
shared = pickle.load(pkl)
pkl.close()
argsdict = vars(args)
argsdict.update(shared)

# make sure the database has been copied, as usual
if args.dbl == None:
	dblist = args.db + '/list'
else:
	dblist = args.dbl
#if sub.call(['perl', '-w', SB + '/copyDBLocally.pl', '-n', removePath(args.db), '-l', args.db+'/list']):
#	print('Error copying the master DB locally - aborting!')
#	exit(1)

# open the sequence database
PATH_seqdb = '/home/grigoryanlab/home/jack/from-fan/Data/searchDB/statistics/bc-30-sc-20141022.peprm2.db'
seqdb = shelve.open(PATH_seqdb, 'r')


# run MASTER search
def runMaster(pdb, rmsdcut):
	Master.createPDS(type='query',
					 pdb=pdb,
				 	)
	Master.masterSearch(query=changeExt(pdb, 'pds'),
						targetList=dblist,
						seqOut=changeExt(pdb, 'seq'),
						matchOut=changeExt(pdb, 'match'),
						rmsdCut=rmsdcut,
						topN=args.c1
						)

if not isinstance(args.p, list):
	args.p = [args.p]

for i in range(len(args.p)):
	frag = args.p[i]
	os.system('cp ' + odir+'/'+frag + ' ./')
	runMaster(frag, float(args.rmsd[i]))
	seqf = changeExt(frag, 'seq')
	matchf = changeExt(frag, 'match')
	if len(args.homo) > 0:
		rowsremain_nh = Analyze.removeHomo(matchf, args.homo) # dry run
		rowsremain_nr = Stability.removeRedundancy(matchf, args.he, int(args.crind[i]), seqdb, usedRows=rowsremain_nh)
	else:
		rowsremain_nr = Stability.removeRedundancy(matchf, args.he, int(args.crind[i]), seqdb)

	# write result
	nseqf, nmatchf = args.he + '_' + seqf, args.he + '_' + matchf
	with open(nseqf, 'w') as nr_sfh, open(nmatchf, 'w') as nr_mfh, open(seqf) as sfh, open(matchf) as mfh:
		n = 0
		for ml in mfh:
			if n in rowsremain_nr:
				nr_mfh.write(ml)
			n+=1
		n = 0
		for sl in sfh:
			if n in rowsremain_nr:
				nr_sfh.write(sl)
			n+=1
	# move files
	shutil.move(nseqf, odir+'/'+nseqf)
	shutil.move(nmatchf, odir+'/'+nmatchf)

Cluster.destroyLocalSpace(ldir)
os.chdir(odir)
os.system('touch .finished.'+args.ncon)
