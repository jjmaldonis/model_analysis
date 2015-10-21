from subproc import run_subproc

def normal_distribution(mu, sigma):
    return lambda x: 1.0/(sigma*sqrt(2*pi))*exp(-(x-mu)*(x-mu)/(2*sigma*sigma))

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

def save_obj(obj, filename ):
    try:
        cPickle
    except:
        import cPickle
    with open(filename, 'wb') as f:
        cPickle.dump(obj, f, cPickle.HIGHEST_PROTOCOL)

def load_obj(filename):
    try:
        cPickle
    except:
        import cPickle
    with open(filename, 'rb') as f:
        return cPickle.load(f)
