import sys, os
from model import Model
from voronoi_3d import voronoi_3d,print_data
from categorize_vor import vor_stats, categorize_atoms
from vor import fortran_voronoi_3d

def highlight(wholemodelfile,submodelfile,outfile=None):
    """ Sub-model is returned """
    wholem = Model(wholemodelfile)
    subm = Model(submodelfile)
    for atom in subm.atoms:
        if(atom in wholem.atoms):
            wholem.atoms[wholem.atoms.index(atom)].z = 0
        else:
            print("Couldn't find atom {0} in large model!".format(atom))
    if(outfile == None):
        outfile = 'highlight.real.xyz'
    #wholem.write_real_xyz(outfile)
    return sumb

#def vp_analyze(wholemodelfile,sm):
def vp_analyze(wm,sm,vpparamfile,outfile=None):
    """ sm is modified and returned """
    #wm = fortran_voronoi_3d(wholemodelfile,3.5)
    for i,atom in enumerate(sm.atoms):
        if(atom in wm.atoms):
            sm.atoms[i].vp = wm.atoms[wm.atoms.index(atom)].vp
        else:
            print("Couldn't find atom {0} in large model!".format(atom))

    categorize_atoms(sm,vpparamfile)
    vs = vor_stats(sm) # Prints what you probably want

    vp_list = list(sm.vp_dict)
    vp_list.append("Undef")
    for i,atom in enumerate(sm.atoms):
        if(atom in wm.atoms):
            #print(atom.vp.type, wm.atoms[wm.atoms.index(atom)].vp.type)
            sm.atoms[i].z = vp_list.index(atom.vp.type) + 1
    #print(vp_list)
    print("Fraction of xtal/ico: {0}".format(float(vs["Crystal-like"]["Total"])/float(vs["Icosahedra-like"]["Total"] + vs["Full-icosahedra"]["Total"])))
    if(outfile == None):
        outfile = 'vp_colorized.real.xyz'
    #sm.write_real_xyz(outfile)
    #print_data(sm)


def main():
    wholemodelfile = sys.argv[1]
    submodelfile = sys.argv[2]
    vpparamfile = sys.argv[3]

    wm = Model(wholemodelfile)
    voronoi_3d(wm,3.5)
    try:
        sm = Model(submodelfile)
        #highlight(wholemodelfile,submodelfile)
        vp_analyze(wm, sm, vpparamfile)
    except:
        for file in sorted(os.listdir(submodelfile)):
            print("Processing file {0}...".format(file))
            #highlight(wholemodelfile,file)
            sm = Model(file)
            vp_analyze(wm, sm, vpparamfile)



        




if __name__ == '__main__':
    main()
