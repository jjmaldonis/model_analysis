import sys
import os
import subprocess
import shutil
import numpy as np
from collections import defaultdict
from model import Model

def main():
    models_prefix = sys.argv[1]
    print("Collecting paths and models...")
    end = models_prefix.rfind('/')+1
    if(end == 0):
        path = './'
    else:
        path = models_prefix[:models_prefix.rfind('/')+1]
    print(path)
    files = os.listdir(path)
    modelfiles = []
    for file in files:
        if models_prefix[models_prefix.rfind('/')+1:] in file:
            modelfiles.append(path+file)
    print("Sorting by filename...")
    modelfiles.sort()
    print('Last model to be counted is {0}'.format(modelfiles[-1]))

    os.chdir(path)
    i = 0
    for model in modelfiles:
        #if(i%100 == 0):
        shutil.copyfile(model,'/home/jjmaldonis/model.xyz')
        model = model[model.rfind('/')+1:]
        print(model)
        #try:
        #    subprocess.call(['femsim', '/home/jjmaldonis/model.xyz'])
        #    os.rename('vk.txt','vk_{0}txt'.format(model.strip().split()[-1][:-4]))
        #except:
        #    pass
        i += 1
    os.remove('/home/jjmaldonis/model.xyz')

if __name__ == "__main__":
    main()
