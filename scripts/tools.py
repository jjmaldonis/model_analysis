#import cPickle

def drange(start, stop, step):
    r = start
    if(stop > start):
        while r < stop:
            yield r
            r += step
    else:
        while r > stop:
            yield r
            r += step

def run_subproc(args):
    """ args should be the string that you would normally run from bash """
    import subprocess
    import shlex
    import sys
    #print("Running (via python): {0}".format(args))
    sargs = shlex.split(args)
    p = subprocess.Popen(sargs)#, stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
    #output = []
    #for nextline in iter(p.stdout.readline, ""):
    #    sys.stdout.write(nextline)
    #    output.append(nextline)
    #    sys.stdout.flush()
    #poutput = p.stdout.read()
    #perr = p.stderr.read()
    #preturncode = p.wait()
    #if(preturncode != 0):
    #    print("{0} exit status: {1}".format(args,preturncode))
    #    print("{0} failed: {1}".format(args,perr))
    #return ''.join(output)

def save_obj(obj, filename ):
    with open(filename, 'wb') as f:
        cPickle.dump(obj, f, cPickle.HIGHEST_PROTOCOL)

def load_obj(filename):
    with open(filename, 'rb') as f:
        return cPickle.load(f)

def makewave():
    with open(sys.argv[1]) as f:
        content = [line.strip().split() for line in f]
    content = zip(*content)
    
    of = open(sys.argv[1][:-4]+'.new.txt','w')
    of.write('IGOR\n')
    for row in content:
        of.write('WAVES/N=({0})\t {1}\nBEGIN\n'.format(len(row),row[0]))
        for i,val in enumerate(row):
            if(i > 0):
                of.write("{0}\n".format(val))
        of.write('END\n')
    of.close()

def printTable(matrix,transpose=False):
    # The input should be a matrix with col,row. E.g. matrix[0] is the first column to print
    # The second parameter allows you to transpose your input matrix
    if(not transpose): # then by column
        matrix = [[x if type(x) == type('a') else str(x) for x in col] for col in matrix]
        rows = len(matrix[0])
        cols = len(matrix)
        colSpacing = [ 3+max([len(x) for x in matrix[i]]) for i in range(cols)]
        table = ''
        for r in range(rows):
            for c in range(cols):
                table = table + matrix[c][r] + ' '*(colSpacing[c]-len(matrix[c][r]))
            table = table + '\n'
        return table.strip()
    else:
        matrix = [[row[i] for row in matrix] for i in range(len(matrix[0]))]
        matrix = [[x if type(x) == type('a') else str(x) for x in col] for col in matrix]
        rows = len(matrix[0])
        cols = len(matrix)
        colSpacing = [ 3+max([len(x) for x in matrix[i]]) for i in range(cols)]
        table = ''
        for r in range(rows):
            for c in range(cols):
                table = table + matrix[c][r] + ' '*(colSpacing[c]-len(matrix[c][r]))
            table = table + '\n'
        #matrix = [[x if type(x) == type('a') else str(x) for x in row] for row in matrix]
        #cols = len(matrix[0])
        #rows = len(matrix)
        #rowSpacing = [ 3+max([len(x) for x in matrix[i]]) for i in range(rows)]
        #table = ''
        #for c in range(rows):
        #    for r in range(cols):
        #        table = table + matrix[c][r] + ' '*(rowSpacing[c]-len(matrix[c][r]))
        #    table = table + '\n'
        return table.strip()
