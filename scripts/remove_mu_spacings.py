import sys, os

def main():
    models_prefix = sys.argv[1]

    print("Collecting paths and models...")
    if( '/' in models_prefix):
        path = models_prefix[:models_prefix.rfind('/')+1]
    else:
        path = './'
    files = os.listdir(path)
    modelfiles = []
    for file in files:
        if models_prefix[models_prefix.rfind('/')+1:] in file:
            modelfiles.append(path+file)
    for file in modelfiles:
        #print(file,''.join(file.strip().split()))
        os.rename(file,''.join(file.strip().split()))


if __name__ == '__main__':
    main()
