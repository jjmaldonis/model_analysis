import os
import subprocess

def main():
    for name in os.listdir(os.getcwd()):
        if name[len(name)-4:] == '.xyz':
            #os.system('sed -i "2 c\\%s" param_file.in' % name)
            os.system('qsub submit.sh '+name)
            #subprocess.call(['sed', "'2 c\name' param_file.in"])
            #subprocess.call(['qsub', 'femsim '+name])
            #subprocess.call(['echo', 'femsim '+name])

    



if __name__ == "__main__":
    main()
