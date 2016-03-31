import subprocess
import shlex
import sys
import os


def run_subproc(args, verbose=True):
    """ args should be the string that you would normally run from bash """
    print("Running (via python): {0}".format(args))
    sargs = shlex.split(args)
    #FNULL = open(os.devnull, 'w')
    p = subprocess.Popen(sargs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = []
    for nextline in iter(p.stdout.readline, ""):
        sys.stdout.write(nextline)
        output.append(nextline)
        sys.stdout.flush()
    poutput = p.stdout.read()
    perr = p.stderr.read()
    preturncode = p.wait()
    if(preturncode != 0):
        print("{0} exit status: {1}".format(args,preturncode))
        print("{0} failed: {1}".format(args,perr))
    return ''.join(output)


def main():
    import sys
    command = ' '.join(sys.argv[1:])
    run_subproc(command)


if __name__ == "__main__":
    main()
