

from model import Model
import sys




def main():
    modelfile = sys.argv[1]
    m = Model(modelfile)
    extension = modelfile[-3:]

    modelfile = modelfile[modelfile.rfind('/')+1:]
    modelfile = modelfile[:-4]+extension
    #print(modelfile)
    m.write_dat('temp.dat')

    of = open("cgm_temp.in",'w')
    of.write("""units		metal
boundary	p p p
atom_style	atomic

read_data temp.dat

pair_style	eam/alloy
pair_coeff	* * /home/jjmaldonis/bin/ZrCuAl2011.eam.alloy Zr Cu Al

thermo_style	custom step etotal fmax fnorm
thermo		1

minimize	1.0e-6 0 1000 10000
write_restart	lmp_emin_exp.restart """)
    of.close()


if __name__== '__main__':
    main()
