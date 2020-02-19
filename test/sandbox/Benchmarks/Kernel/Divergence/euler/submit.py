#!/usr/bin/env python
# Generated with: /cluster/home/fabianw/bin/configJob.py --euler -N 1 -n 1 --fullnode 36 -J euler-div-kernel-512 -W 00:30
import argparse
import sys
import os
import json
import subprocess
import datetime

def parseArgs():
    parser = argparse.ArgumentParser()

    # strings args
    parser.add_argument('--cwd', type=str, default=os.getcwd(), help='Current working directory')

    return parser.parse_known_args()

if __name__ == "__main__":
    args, optional = parseArgs()
    now = datetime.datetime.now()

    os.environ['JC_JOBWD'] = args.cwd

    if os.path.isfile('./common.sh'):
        dump = 'python -c "import os, json;print(json.dumps(dict(os.environ)))"'
        pipe = subprocess.Popen(['bash', '-c', 'source ./common.sh && %s' % (dump)], stdout=subprocess.PIPE)
        env  = json.loads(pipe.stdout.read())
        os.environ = env

        # get job ID for JC submissions
        os.environ['JC_LOGHOME'] = "."
        countfile = "%s/.count" % (os.environ['JC_LOGHOME'])
        if os.path.isfile(countfile):
            with open(countfile, 'r') as cf:
                count = int(cf.readline())
            with open(countfile, 'w') as cf:
                cf.write('%d' % (count+1))
        else:
            with open(countfile, 'w') as cf:
                cf.write('2')
            count = 1
        os.environ['JC_JOBID']   = '%04d' % count
        os.environ['JC_JOBDATE'] = now.strftime('%Y-%m-%d')
        os.environ['JC_JOBTIME'] = now.strftime('%H:%M')

        tmplogfile = "%s/.<euler-div-kernel-512_%s_%s_%s>.tmp" % (
            os.environ['JC_LOGHOME'],
            os.environ['JC_JOBDATE'],
            os.environ['JC_JOBTIME'],
            os.environ['JC_JOBID']
            )
        with open(tmplogfile, 'w') as log:
            log.write('%s\n%s\n%s\n' % (
                os.environ['JC_JOBDATE'],
                os.environ['JC_JOBTIME'],
                os.environ['JC_JOBID']
                )
            )
        os.environ['JC_JOBTMP'] = tmplogfile
        submission = "bsub -cwd ${JC_JOBWD} -E \"${JC_PROLOG} '${JC_JOBTMP}'\" -Ep \"${JC_EPILOG} '${JC_JOBTMP}'\" -R fullnode < ${JC_JOB}"
        subprocess.Popen(submission, shell=True, env=dict(os.environ)).wait()
    else:
        submission = "bsub -cwd ${JC_JOBWD} -R fullnode < ./job.sh"
        subprocess.Popen(submission, shell=True, env=dict(os.environ)).wait()
