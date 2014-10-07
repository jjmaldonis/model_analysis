def run_subproc(args):
    """ args should be the string that you would normally run from bash """
    import subprocess
    import shlex
    import sys
    print("Running (via python): {0}".format(args))
    sargs = shlex.split(args)
    p = subprocess.Popen(sargs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for nextline in iter(p.stdout.readline, ""):
        sys.stdout.write(nextline)
        sys.stdout.flush()
    poutput = p.stdout.read()
    perr = p.stderr.read()
    preturncode = p.wait()
    if(preturncode != 0):
        print("{0} exit status: {1}".format(args,preturncode))
        print("{0} failed: {1}".format(args,perr))

def main():
    import sys
    command = ' '.join(sys.argv[1:])
    run_subproc(command)

if __name__ == "__main__":
    main()
