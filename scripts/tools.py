import cPickle

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

def printTable(matrix):
	# The input should be a matrix with col,row. E.g. matrix[0] is the first column to print
	matrix = [[x if type(x) == type('a') else str(x) for x in col] for col in matrix]
	rows = len(matrix[0])
	cols = len(matrix)
	colSpacing = [ 3+max([len(x) for x in matrix[i]]) for i in range(cols)]
	table = ''
	for r in range(rows):
		for c in range(cols):
			table = table + matrix[c][r] + ' '*(colSpacing[c]-len(matrix[c][r]))
		table = table + '\n'
	return table
