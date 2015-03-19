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


