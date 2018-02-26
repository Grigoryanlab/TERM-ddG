import time
from getpass import getuser
# functions about using anthill cluster
from General import *

class jobOnCluster:
    def __init__(self, cmds, myid='', label=''):  # can skip myid and label if want to be simple
        self.cmds = cmds  # this is a list
        self.running = 0
        self.jobid = ''
        self.tried = 0
        self.myid = myid  # this is an easy identifier (usually a string) for user to determine where problems come from
        self.label = label  # this is a file indicate whether the job has finished

    def submit(self, time):
        self.jobid = qsub(self.cmds, hrs=int(time))
        self.running = 1
        self.tried += 1

    def checkjob(self):
        if checkJobRun(self.jobid) == 0:
            self.running = 0

    def checkfinish(self):  # label is a file when it exists it means the job has finished
        if os.path.isfile(self.label) == True:
            return 1
        else:
            return 0


def qsub(cmds, user='', fileName=None, mem=2, hrs=3, ironfs=False, opts=[], maxJobs=2500, avoidNode = []):
    assert mem > 0, 'memory must be greater than 0'
    assert isinstance(hrs, int) and hrs >= 1, 'hrs must be and integer greater than 0'

    if fileName == None:
        fileName = 'script-%.15f.sh' % time.time()

    with open(fileName, 'w') as fH:
        fH.write('#!bin/bash\n\n')
        fH.write('#$ -cwd\n')
        fH.write('#$ -j y\n')
        fH.write('#$ -V\n')
        fH.write('#$ -l vf=' + str(mem) + 'G\n')
        fH.write('#$ -l h_rt=' + str(hrs - 1) + ':59:59\n')

        if len(avoidNode) > 0:
            fH.write('#$ -l hostname=\'!' + '|'.join(avoidNode) + '\'\n')

        if ironfs == True:
            fH.write('#$ -l ironfs\n')

        for opt in opts:
            if opt != None:
                fH.write('\n' + opt + '\n')

        for cmd in cmds:
            fH.write('\n' + cmd)

    qsubCall = ['qsub', fileName]
    if isinstance(maxJobs, int):
        while numJobs(user) >= maxJobs:
            time.sleep(2)

    p = sub.Popen(qsubCall, stdout=sub.PIPE)
    std_out = p.communicate()[0]
    print(std_out)
    jobid = parseJobID(std_out)
    return jobid  # will return job id


def waitJobs(jobs, type = 'list', sleep_time=120, rerun_time=24, giveup_time=3, subdir = True):
    success = True
    odir = os.getcwd()
    if type == 'list':
        assert isinstance(jobs, list)
        while len(jobs)>0:
            for j in jobs:
                j.checkjob()
                if j.running == 0:
                    jobs.remove(j)
                    if j.checkfinish() == 0:
                        print('Job about ' + j.myid + ' may have died ...')
                        if j.tried > giveup_time:
                            print('have failed %d times, give up ...' % (giveup_time+1))
                            success = False
                            continue
                        print('resubmitting ... ')
                        if subdir:
                            os.chdir(j.myid)
                        j.submit(rerun_time)
                        jobs.append(j)
                        os.chdir(odir)
            print('Running, '+ str(len(jobs)) + ' jobs left ...')
            time.sleep(sleep_time)
    if type == 'dict':
        assert isinstance(jobs, dict)
        while len(jobs)>0:
            for k in jobs.keys():
                j = jobs[k]
                j.checkjob()
                if j.running == 0:
                    jobs.pop(k)
                    if j.checkfinish() == 0:
                        print('Job about ' + j.myid + ' may have died ...')
                        if j.tried > giveup_time:
                            print('have failed %d times, give up ...' % (giveup_time+1))
                            success = False
                            continue
                        print('resubmitting ... ')
                        if subdir:
                            os.chdir(j.myid)
                        j.submit(rerun_time)
                        jobs[k] = j
                        os.chdir(odir)
            print('Running, '+ str(len(jobs)) + ' jobs left ...')
            time.sleep(sleep_time)
    return success



def numJobs(user=''):
    if user == '':
        user = getuser()
    while True:
        out = sub.check_output('qstat -u ' + user + ' | grep "" | wc -l', stderr=sub.STDOUT, shell=True)
        try:
            njobs = int(out) - 2
            break
        except:
            time.sleep(5)
            continue
    return max(0, njobs)


def parseJobID(string):
    id = re.search('job\s(\d+)', string)
    if id != None:
        return id.group(1)


def checkJobRun(jid):
    while True:
        try:
            out = sub.check_output('qstat', stderr=sub.STDOUT, shell=True)
            break
        except:
            time.sleep(5)
            continue
    if out == None:
        return 0
    out = out.split('\n')
    for l in out:
        arr = l.split()
        if arr == []:
            continue
        if (arr[0] == jid) and (not arr[4] in ['dr']):
            return 1
    return 0


def createLocalSpace(user=''):
    if user == '':
        user = getuser()
    while True:
        rint = str(random.randint(0, 1000000)).rjust(6, '0')
        ldir = '/data/scratch/'+ user + '/' + rint
        if os.path.isdir(ldir):
            continue
        else:
            break
    os.makedirs(ldir)
    return ldir


def destroyLocalSpace(ldir):
    shutil.rmtree(ldir)
