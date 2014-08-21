import sys, os
from basic_model import Model
from hutch import Hutch
from atom import Atom
from atom_graph import AtomGraph
from voronoi_3d import voronoi_3d
from categorize_vor import load_param_file, set_atom_vp_types, vor_stats

def nnnic(atom,model,cluster):
    """ number of nearest neighbors in cluster 
    @param atom is the atom around which we want to find its neighbors
    @param model is the model that the cluster belongs to
    @param cluster is the cluster of atoms, in list form """
    n = 0
    for atomi in atom.neighs:
        if atomi in cluster:
            n += 1
    return n
    

def read_intensity_file(intensityfile):
    #maxx = 0.0
    print("Reading intensity file {0}".format(intensityfile))
    with open(intensityfile) as f:
        for l,line in enumerate(f):
            line = line.strip().split()
            if(l == 0):
                npixx, npixy, npixz = tuple([int(x) for x in line])
                if(npixx == npixy == npixz):
                    npix = npixx
                else:
                    pass
                intensities = [[[0.0 for i in range(npix)] for i in range(npix)] for i in range(npix)]
            else:
                for i,x in enumerate(line):
                    x = float(x)
                    #if(x > maxx): maxx = x
                    #intensities[i%npix][int(i/npix)][l-1] = x
                    intensities[int(i/npix)][i%npix][l-1] = x  # Correct
                print("Read in line {0} of file {1} in python.".format(l-1,intensityfile))
    #print("Max intensity = {0}".format(maxx))
    return intensities

def main():
    """ Columns that I should have:
        % atoms that are xtal-like (i.e. what % of all the xtal-like atoms are in this sub-model?)
        % atoms that are ico-like
        % atoms that are mixed
        % atoms that are undefined
        Can you see planes?
        Is the spot visible in the sub-model's FT?
        What is the ratio of the max intensity of the spot in the original FT / that in the sub-model's FT? (Do this on a per-atom basis to account for the different number of atoms.)"""

    # paramfile, jobid, original ft, main ft direc, VP categories paramfile
    paramfile = sys.argv[1] # Paramfile
    with open(paramfile) as f:
        params = f.readlines()
    params = [line.strip() for line in params]
    modelfile = params[0]
    num_spots = int(params[1])
    jobid = sys.argv[2] # jobid
    orig_ft = sys.argv[3] # original ft
    main_direc = sys.argv[4] # spot ft direc
    spot_ft_direc = main_direc+'/spot_fts/'
    spot_fts = sorted(os.listdir(spot_ft_direc))
    spot_fts = [spot_ft_direc+line for line in spot_fts]
    vp_paramfile = sys.argv[5] # VP paramfile for categorizing

    spot_ints = []
    orig_intensities = read_intensity_file(orig_ft)
    for intensityfile in spot_fts:
        try:
            i = int(intensityfile[intensityfile.index('_spot')+5:intensityfile.index('_spot')+6])
        except:
            i = int(intensityfile[intensityfile.index('_spot')+5])
        print("Intensity file: {0} => i={1}".format(intensityfile,i))
        j = i*7 +2
        x0,y0,z0    = tuple([float(x) for x in params[j+1].split()])
        sx,sy,sz    = tuple([float(x) for x in params[j+2].split()])
        cxy,cxz,cyz = tuple([float(x) for x in params[j+3].split()])
        xi,xf,xc    = tuple([float(x)*2 for x in params[j+4].split()][:3])
        yi,yf,yc    = tuple([float(x)*2 for x in params[j+5].split()][:3])
        zi,zf,zc    = tuple([float(x)*2 for x in params[j+6].split()][:3])
        ft_intensities = read_intensity_file(intensityfile)
        maxoi = 0.0
        maxi = 0.0
        for x in range(xi-1,xf):
            for y in range(yi-1,yf):
                for z in range(zi-1,zf):
                    #print(x,y,z,orig_intensities[x][y][z])
                    if( ft_intensities[x][y][z] > maxi ):
                        maxi = ft_intensities[x][y][z]
                    if( orig_intensities[x][y][z] > maxoi ):
                        maxoi = orig_intensities[x][y][z]
        print("{0}/{1}={2} I still need to rescale by number of atoms".format(maxi,maxoi,maxi/maxoi))
        spot_ints.append((i,maxi,maxoi,maxi/maxoi))
        print(spot_ints)

    mm = Model(modelfile)
    voronoi_3d(mm,3.5)
    vp_dict = load_param_file(vp_paramfile)
    set_atom_vp_types(mm,vp_dict)
    vor_stats(mm)
    mtype_dict = {}
    for atom in mm.atoms:
        mtype_dict[atom.vp.type] = mtype_dict.get(atom.vp.type,0) + 1
    for i in range(0,num_spots):
        j = i*7 +2
        submodelfile = params[j] + jobid + '.xyz'
        try:
            sm = Model(submodelfile)
            # Assign each atom its correct type
            smtype_dict = {}
            for atom in sm.atoms:
                #print(mm.atoms[mm.atoms.index(atom)])
                atom.vp.type = mm.atoms[mm.atoms.index(atom)].vp.type
                smtype_dict[atom.vp.type] = smtype_dict.get(atom.vp.type,0) + 1
            print(submodelfile)
            #vor_stats(sm)
            for key in smtype_dict:
                print("  Percent of {0} atoms: {1}/{2}= {3}".format(key,smtype_dict[key],mtype_dict[key],round(100.0*smtype_dict[key]/mtype_dict[key])))
                print("  Ratio of xtal-like to ico-like: {0}".format( (100.0*smtype_dict['Crystal-like']/mtype_dict['Crystal-like']) / (100.0*smtype_dict['Icosahedra-like']/mtype_dict['Icosahedra-like'])))
        except:
            pass
    print(spot_ints)

    return
    # sys.argv[1] should be the modelfile in the xyz format
    # sys.argv[2] should be the cutoff desired
    modelfile = sys.argv[1]
    cutoff = float(sys.argv[2])
    ag = AtomGraph(modelfile,cutoff)
    model = Model(modelfile)
    model.generate_neighbors(cutoff)
    submodelfile = sys.argv[3]

    mixedmodel = Model('Mixed atoms',model.lx,model.ly,model.lz, [atom for atom in ag.model.atoms if atom.vp.type == 'Mixed'])
    icolikemodel = Model('Ico-like atoms',model.lx,model.ly,model.lz, [atom for atom in ag.model.atoms if(atom.vp.type == 'Icosahedra-like' or atom.vp.type == 'Full-icosahedra')])
    fullicomodel = Model('Full-icosahedra atoms',model.lx,model.ly,model.lz, [atom for atom in ag.model.atoms if atom.vp.type == 'Full-icosahedra'])
    xtalmodel = Model('Xtal-like atoms',model.lx,model.ly,model.lz, [atom for atom in ag.model.atoms if atom.vp.type == 'Crystal-like'])
    undefmodel = Model('Undef atoms',model.lx,model.ly,model.lz, [atom for atom in ag.model.atoms if atom.vp.type == 'Undef'])
    #mixedmodel.write_cif('mixed.cif')
    #mixedmodel.write_our_xyz('mixed.xyz')
    #icolikemodel.write_cif('icolike.cif')
    #icolikemodel.write_our_xyz('icolike.xyz')
    #fullicomodel.write_cif('fullico.cif')
    #fullicomodel.write_our_xyz('fullico.xyz')
    #xtalmodel.write_cif('xtal.cif')
    #xtalmodel.write_our_xyz('xtal.xyz')
    #undefmodel.write_cif('undef.cif')
    #undefmodel.write_our_xyz('undef.xyz')
    icomixedmodel = Model('ico+mix atoms',model.lx,model.ly,model.lz, mixedmodel.atoms + icolikemodel.atoms)
    #mixedmodel.write_cif('icomixed.cif')
    #mixedmodel.write_our_xyz('icomixed.xyz')
    vpcoloredmodel = Model('vp colored atoms',model.lx,model.ly,model.lz, ag.model.atoms)
    for atom in vpcoloredmodel.atoms:
        if(atom.vp.type == 'Full-icosahedra'):
            atom.z = 1
        elif(atom.vp.type == 'Icosahedra-like'):
            atom.z = 2
        elif(atom.vp.type == 'Mixed'):
            atom.z = 3
        elif(atom.vp.type == 'Crystal-like'):
            atom.z = 4
        elif(atom.vp.type == 'Undef'):
            atom.z = 5
    #vpcoloredmodel.write_cif('vpcolored.cif')
    #vpcoloredmodel.write_our_xyz('vpcolored.xyz')
    subvpcoloredmodel = Model(submodelfile)
    for atom in subvpcoloredmodel.atoms:
        atom.z = vpcoloredmodel.atoms[ag.model.atoms.index(atom)].z
    subvpcoloredmodel.write_cif('subvpcolored.cif')
    subvpcoloredmodel.write_our_xyz('subvpcolored.xyz')
    return

    golden = False

    #cluster_prefix = 'ico.t3.'
    #cluster_types = 'Icosahedra-like', 'Full-icosahedra' # Need this for final/further analysis

    #cluster_prefix = 'fi.t3.'
    #cluster_types = ['Full-icosahedra'] # Need this for final/further analysis

    cluster_prefix = 'xtal.t3.'
    cluster_types = 'Crystal-like' # Need this for final/further analysis

    #cluster_prefix = 'mix.t3'
    #cluster_types = ['Mixed'] # Need this for final/further analysis

    #cluster_prefix = 'undef.t3'
    #cluster_types = ['Undef'] # Need this for final/further analysis

    # Decide what time of clustering you want to do
    #clusters = ag.get_clusters_with_n_numneighs(cutoff,5,cluster_types) #Vertex
    #clusters = ag.get_vertex_sharing_clusters(cutoff,cluster_types) #Vertex
    #clusters = ag.get_edge_sharing_clusters(cutoff,cluster_types) #Edge
    #clusters = ag.get_face_sharing_clusters(cutoff,cluster_types) #Face
    clusters = ag.get_interpenetrating_atoms(cutoff,cluster_types) #Interpenetrating
    #clusters = ag.get_interpenetrating_clusters_with_neighs(cutoff,cluster_types) #Interpenetrating+neighs
    #clusters = ag.get_connected_clusters_with_neighs(cutoff, cluster_types) #Connected (vertex) + neighs
    #v,e,f,i = ag.vefi_sharing(cluster_types)
    #print("V: {0}  E: {1}  F: {2}  I: {3}".format(int(v),int(e),int(f),int(i)))

    orig_clusters = clusters[:]
    # Print orig clusters
    j = 0
    for i,cluster in enumerate(clusters):
        print("Orig cluster {0} contains {1} atoms.".format(i,len(cluster)))
        if(golden):
            for atom in cluster:
                if(atom.vp.type in cluster_types):
                    atom.z = 0
        # Save cluster files
        cluster_model = Model("Orig cluster {0} contains {1} atoms.".format(i,len(cluster)),model.lx, model.ly, model.lz, cluster)
        cluster_model.write_cif('{1}cluster{0}.cif'.format(i,cluster_prefix))
        cluster_model.write_our_xyz('{1}cluster{0}.xyz'.format(i,cluster_prefix))

    allclusters = []
    for cluster in clusters:
        for atom in cluster:
            if(atom not in allclusters):
                allclusters.append(atom)
                #if(atom.vp.type in cluster_types): print('  {0}\t{1}'.format(atom,atom.vp.type))
    allclusters = Model("All clusters.",model.lx, model.ly, model.lz, allclusters)
    allclusters.write_cif('{0}allclusters.cif'.format(cluster_prefix))
    allclusters.write_our_xyz('{0}allclusters.xyz'.format(cluster_prefix))
    print("{0}allclusters.cif and {0}allclusters.xyz contain {1} atoms.".format(cluster_prefix, allclusters.natoms))

    if(not golden):
        x_cluster = []
        for i,atom in enumerate(model.atoms):
            if atom not in allclusters.atoms:
                x_cluster.append(atom)
        x_cluster = Model("Opposite cluster of {0}".format(cluster_prefix),model.lx, model.ly, model.lz, x_cluster)
        x_cluster.write_cif('{0}opposite.cif'.format(cluster_prefix))
        x_cluster.write_our_xyz('{0}opposite.xyz'.format(cluster_prefix))
        print('{0}opposite.cif and {0}opposite.xyz contain {1} atoms.'.format(cluster_prefix, x_cluster.natoms))
    
    if(False): # Further analysis
        cn = 0.0
        for atom in model.atoms:
            cn += atom.cn
        cn /= model.natoms

        vpcn = 0.0
        count = 0
        for atom in ag.model.atoms:
            if( atom.vp.type in cluster_types ):
                vpcn += atom.cn
                count += 1
        vpcn /= count

        natomsinVPclusters = allclusters.natoms # Number of atoms in VP clusters
        nVPatoms = count # Number of VP atoms
        numsepVPatoms = nVPatoms * vpcn # Number of atoms in VP clusters if all clusters were separated
        maxnumatoms = model.natoms # Max number of atoms in VP clusters if all clusters were separated but still within the model size

        print('Average CN is {0}'.format(cn))
        print('Average CN of VP atoms is {0}'.format(vpcn))
        print('# atoms in all clusters: {0}. # VP atoms * vpcn: {1}. # VP atoms: {2}'.format(natomsinVPclusters,numsepVPatoms,nVPatoms))
        print('~ Number of VP that can fit in the model: {0}'.format(maxnumatoms/vpcn))
        print('Ratio of: (# atoms involved in VP clusters)/(# atoms involved in VP clusters if all clusters were completely separated):                          {0}%  <--- Therefore {1}% sharing.'.format(round(float(natomsinVPclusters)/(numsepVPatoms)*100.0,3),100.0-round(float(natomsinVPclusters)/(numsepVPatoms)*100.0,3)))
        print('Ratio of: (# atoms involved in VP clusters)/(# atoms involved in VP clusters if all clusters were separated as much as possible within the model): {0}%  <--- Therefore {1}% sharing.'.format(round(float(natomsinVPclusters)/min(numsepVPatoms,maxnumatoms)*100.0,3),100.0-round(float(natomsinVPclusters)/min(numsepVPatoms,maxnumatoms)*100.0,3) if numsepVPatoms < maxnumatoms else round(float(natomsinVPclusters)/min(numsepVPatoms,maxnumatoms)*100.0,3)))

        vor_instance = Vor()
        vor_instance.runall(modelfile,cutoff)
        index = vor_instance.get_indexes()
        vor_instance.set_atom_vp_indexes(model)
        vp_dict = categorize_vor.load_param_file('/home/jjmaldonis/model_analysis/scripts/categorize_parameters_iso.txt')
        atom_dict = categorize_vor.generate_atom_dict(index,vp_dict)
        categorize_vor.set_atom_vp_types(model,vp_dict)
        # Count the number of common neighbors in each of the VP
        vp_atoms = []
        for atom in model.atoms:
            if(atom.vp.type in cluster_types):
                vp_atoms.append(atom)
        common_neighs = 0.0
        atom_pairs = []
        for atomi in vp_atoms:
            for atomj in vp_atoms:
                if(atomi != atomj):
                    if(atomi in atomj.neighs and [atomi,atomj] not in atom_pairs and [atomj,atomi] not in atom_pairs):
                        common_neighs += 1
                        atom_pairs.append([atomi,atomj])
                    #if(atomj in atomi.neighs): common_neighs += 0.5
                    for n in atomi.neighs:
                        if(n in atomj.neighs and [n,atomj] not in atom_pairs and [atomj,n] not in atom_pairs):
                            common_neighs += 1
                            atom_pairs.append([n,atomj])
                    #for n in atomj.neighs:
                    #    if(n in atomi.neighs): common_neighs += 0.5
        # Now common_neighs is the number of shared atoms
        #print(common_neighs)
        print('Percent shared based on common neighsbors: {0}'.format(100.0*common_neighs/natomsinVPclusters))
        
        # Debugging stuff, mainly for get_connected_clusters_with_neighs function
        #print(clusters[0][57], clusters[0][57].vesta())
        #for atom in ag.model.atoms:
        #cm = Model('temp',model.lx,model.ly,model.lz,clusters[0])
        #cm.generate_neighbors(3.5)
        #for i,atom in enumerate(cm.atoms):
        #    found = False
        #    for n in atom.neighs:
        #        if(n.vp.type in cluster_types):
        #            found = True
        #    if(not found and atom.vp.type not in cluster_types):
        #        #print('Found an atom that isnt connected to a VP type! {0} {1}'.format(allclusters.atoms.index(atom),atom.neighs))
        #        print('Found an atom that isnt connected to a VP type! {0} {1} {2} {3}'.format(i+1,atom,atom.neighs,atom.vp.type))
        #    #if(atom.neighs == []):
        #    #    for atom2 in ag.model.atoms:
        #    #        if(atom in atom2.neighs):
        #    #            print('Found an atom with empty neighbors as a neighbor of another! {0} {1} {2}'.format(atom,atom2,atom2.neighs))
        #    #if(clusters[0][57] in atom.neighs):
        #        #print(atom,atom.neighs)


if __name__ == "__main__":
    main()
