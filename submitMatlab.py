__author__ = 'fanzheng'

from General import *
import PDB, Cluster, Stability
from scipy.io import loadmat

SB = selfbin(sys.argv[0])
par = argparse.ArgumentParser()
par.add_argument('--l', required = True, help = 'a list of mutations')
par.add_argument('--path', default= SB +'/MATLAB', help = 'the path for MATLAB scripts')
par.add_argument('--head', default = 'nr', help = 'header of the .seq files')
par.add_argument('--ext', required = True, help = 'the extension of the output .mat file')
par.add_argument('--lamb', default = '[]', help = 'control the stength of constraint')
par.add_argument('--noprior', action = 'store_true', help = 'whether use prior')
par.add_argument('--puse', action = 'store_true', help = 'if true, only load the data and output parameter usage')
par.add_argument('--conpot', default = '', help = 'file for contact potential')
par.add_argument('--o', required = True, help = 'output file, for the result of prediction')
par.add_argument('--oo', action = 'store_true', help = 'if true, the matlab files already exists and only output the results')
par.add_argument('--scan', action = 'store_true', help = 'if true, will output all amno acids' )
par.add_argument('--bg', action='store_true', help = 'if true, will modify EPs by background parameters')
args = par.parse_args()

def outputscore(mutlist, output, scan =0, bg =0):
    if bg:
        aafreq= np.array(PDB.aaFreq)
    with open(mutlist) as mutf, open(output, 'w') as ofh:
        for l in mutf:
            l = l.strip()
            mut = Stability.getMutation(l)
            consf = mut.dir + '/' + mut.dir + '.conlist'
            cons =[]
            with open(consf) as cf:
                for cfl in cf:
                    cinfo = cfl.strip().split()
                    conname, conres, cond = cinfo[2], cinfo[-1], float(cinfo[3])
                    if not os.path.isfile(mut.dir+ '/'+mut.dir + '_' + conname.replace(',', '') + '.pdb'):
                        continue
                    if cond >= 0.02:
                        cons.append([conname, conres, cond])
            cons = sorted(cons, key = lambda x :x[0])

            matf = mut.dir + '/' + mut.dir + '.' + args.ext + '.mat'
            if not os.path.isfile(matf):
                ofh.write(l + '\tNA\n')
                continue
            optparams = loadmat(matf)['optParams']
            assert (optparams.shape[0] - 20) / 400 == len(cons), mut.dir
            # if (optparams.shape[0] - 20) / 400 != len(cons):
            #     ofh.write(l + '\tNA\n')
            #     continue
            selfdiff = optparams[PDB.aaindex[mut.m]] - optparams[PDB.aaindex[mut.w]]
            if bg:
                selfdiff += np.log(aafreq[PDB.aaindex[mut.m]]) - np.log(aafreq[PDB.aaindex[mut.w]])
            pairdiff = []
            if scan:
                scanscore = optparams[0:20].T[0]
                if bg:
                    scanscore += np.log(aafreq)
                scanscore -= scanscore[PDB.aaindex[mut.w]]
            for i in range(len(cons)):
                paramsi = optparams[20 + 400*i: 20+400*(i+1)].reshape(20,20)
                conind = PDB.aaindex[PDB.t2s(cons[i][1])]
                pairdiff_i = paramsi[conind, PDB.aaindex[mut.m]] - paramsi[conind, PDB.aaindex[mut.w]]
                pairdiff.append(pairdiff_i)
                # outstr = l + '\t' + cons[i][0] + '\t' + format(float(pairdiff_i), '.3f') # test
                # ofh.write(outstr + '\n') # test
                if scan:
                    scancon = np.array([paramsi[conind, x] for x in range(20)])
                    # use wt aa as the reference state
                    ref = paramsi[conind, PDB.aaindex[mut.w]]
                    scancon -= ref
                    scanscore += scancon
            if not scan:
                outstr = '\t'.join([l, format(float(selfdiff), '.3f'), format(float(sum(pairdiff)), '.3f')])
            else:
                outstr = '\t'.join([l] + [format(x, '.3f') for x in scanscore] )
            ofh.write(outstr + '\n')
    # exit the program
    exit(0)

scan, bg = 0, 0
if args.scan:
    scan = 1
if args.bg:
    bg = 1

if args.oo:
    outputscore(args.l, args.o, scan, bg) # this will ignore all the remaining lines

assert args.conpot != None
assert args.head != None

dirs = os.walk('.').next()[1]
dirs.sort()

dirs =[]
with open(args.l) as mutf:
    for l in mutf:
        mut = Stability.getMutation(l.strip())
        if not mut.dir in dirs:
            dirs.append(mut.dir)

noprior, puse = '[]', '[]'
if args.noprior:
    noprior = '1'
if args.puse:
    puse = '1'

jobs = []
for d in dirs:
    if os.path.isfile(d + '/' + d + '.' + args.ext + '.mat'):
        continue
    cmds = []
    cmds.append(' '.join(['matlab', '-nodisplay', '-nojvm', '-r', '"addpath(\'' + args.path + '\');main(\'' + d + '\',\'' + args.head + '\',\'' + d + '/' + d + '.' + args.ext + '\',' + args.lamb + ',' + noprior + ',' + puse + ',\'' + args.conpot + '\');\"']))

    if not args.puse:
        job = Cluster.jobOnCluster(cmds, d, d + '/' + d + '.' +args.ext+'.mat')
    else:
        job = Cluster.jobOnCluster(cmds, d, d +'/' +d +'.'+args.ext+'.puse')
    jobs.append(job)
    job.submit(3)

Cluster.waitJobs(jobs, giveup_time=0, sleep_time=5)

# score
outputscore(args.l, args.o, scan, bg)

