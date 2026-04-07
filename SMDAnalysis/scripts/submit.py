import subprocess
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument("-f", help="job submission script")
parser.add_argument("-N", help="number of job you want to submit")
parser.add_argument("-crj", help="Currently Running Job number")
parser.set_defaults(crj = 0)
args = parser.parse_args()

def get_jobnumber(lines):
    for line in lines:
        if line.startswith('Submitted'):
            slines = line.split()
            assert slines[0] == "Submitted", "Submitted"
            assert slines[1] == "batch", "batch"
            assert slines[2] == "job", "job"
            job = int(slines[3])
    return job


if int(args.crj) == 0:
    lines = subprocess.run(['sbatch', args.f], stdout=subprocess.PIPE).stdout.decode('utf-8').strip().split('\n')
else:
    lines = subprocess.run(['sbatch', '--dependency=afterany:{}'.format(args.crj), args.f], stdout=subprocess.PIPE).stdout.decode('utf-8').strip().split('\n')

job = get_jobnumber(lines)


for i in range(int(args.N) - 1):
    lines = subprocess.run(['sbatch', '--dependency=afterany:{}'.format(job), args.f], stdout=subprocess.PIPE).stdout.decode('utf-8').strip().split('\n')
    job = get_jobnumber(lines)

