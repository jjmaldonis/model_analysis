There are two ways to run a StructOpt job.

1.  By enterning parameter values
2.  By uploading an input file

The steps below need to be followed in order to run a StructOpt job by entering the parameter values on the StructOpt web application.

1.  Click on the **Home** tab.
2.  Click on **Primary parameters** tab on the left pane of the web page. The primary parameters will be displayed on the right pane of the page.
    1.  The **Structure** parameter specifies the type of structure to be optimized. There are three options to choose for this parameter i.e. Cluster, Crystal and Defect. Select **Cluster** structure.
    2.  The **Optimizer type** parameter specifies the type of optimizer to be used for optimization. There are three options to choose for this parameter i.e. Random, Simulated annealing (SA), Basin hopping (BH) and Genetic algorithm (GA). Select **Genetic algorithm (GA)** optimizer type.
    3.  The **Atom list** parameter allows you to set the number and type of atoms in a simulation. The details include Type of atom, No. of atom, Mass of atom and Chemical potential of atom. One row corresponds to one atom type. In order to add a row, click on the **Add row** button. In order to delete a row, click on **Delete row** button. In the first row, enter **C** as the Type of atom, **6** as No. of atom, **12** as Mass of atom and **-4.000** as Chemical potential of atom. If you want to add a new row, click on the Add row button and enter the Type of atom, No. of atom, Mass of atom and Chemical potential of atom.
    4.  The **Number of atoms** parameter denotes the number of atoms to be used in initiating a simulation. The value of this parameter is auto populated and is equal to the sum of the values in the **No. of atom** text box of the Atom list parameter.
    5.  The **Number of individuals** parameter specifies the number of individuals in a population. Enter **20** as the value for this parameter.
3.  To enter the values for **Output parameters**, click on either the **Next** button on the bottom right of the center pane or on the Output parameters tab on the left pane of the web page. The output parameters will be displayed on the right pane of the page.
    1.  The **Email address** parameter specifies the email address to which the notifications and StructOpt job output results will be mailed. Enter your email address as the value for this parameter.
    2.  The **System name** parameter specifies the user-defined name for the StructOpt job run. Enter **Carbon** as the value for this parameter.
    3.  The **File name** parameter specifies the name of output folder and file for optimizer data. Enter **Carbon** as the value for this parameter.
    4.  The **Genealogy** parameter is a flag indicating if a genealogy of the structures in a population should be written. Select **True** as the value for this parameter.
    5.  The **Output format** parameter specifies the format for outputting fitness data to files. Enter **formationenergy** as the value for this parameter.
    6.  The **All energy file output** parameter is a flag to determine whether to output all energy file for simulation. Select **True** as the value for this parameter.
    7.  The **Best individuals output** parameter is a flag to determine whether to output and track best structures. Select **True** as the value for this parameter.
    8.  The **Number of bests** parameter specifies the number of individuals to keep track of in best list. Enter **40** as the value for this parameter.
    9.  The **Write individual defect** parameter is a flag indicating if atoms in region 1 and 2 for defect calculation should be written. Select **False** as the value for this parameter.
    10.  The **Output vacancies** parameter is a flag indicating if vacancy atoms should be written to the structure output. Select **False** as the value for this parameter.
4.  To enter the values for **General algorithm parameters**, click on either the **Next** button on the bottom right of the center pane or on the General algorithm parameters tab on the left pane of the web page. The general algorithm parameters will be displayed on the right pane of the page.
    1.  The **Genealogy tree** parameter is a flag that indicates if the optimizer should load structures from a previously started simulation. Select **False** as the value for this parameter.
    2.  The **Number of restart atoms** parameter specifies the number of atoms included in region 1 and 2 for defect simulation for loading into a restarted simulation. Enter **0** as the value for this parameter.
    3.  The **Seed** parameter specifies the random number seed generator. Enter **0** as the value for this parameter.
    4.  The **Forcing control** parameter specifies the method for controlling realistic configurations. There are three options to choose for this parameter i.e. Concentration, Chemical potential and Energy bias. Select **Concentration** forcing control.
    5.  The **Debug** parameter specifies the method for debugging a section. There are nine options to choose for this parameter i.e. Main optimizer algorithm, Selection schemes, Mutation schemes, Crossover schemes, Fitness schemes, Energy evaluation algorithm, Find defects function, Setup calculator function and None. Select **None** debug method.
    6.  The **Algorithm type** parameter specifies the type of algorithm. There are four options to choose for this parameter i.e. (λ+μ), Island_Method, (λ,μ) and Island_Method (λ,μ). Select **(λ+μ)** algorithm type.
    7.  The **Finger printing** parameter is a flag determining whether fingerprint function for structures should be calculated. Select **False** as the value for this parameter.
    8.  The **Fixed region** parameter is a flag determining if region 3 should be fixed and forces set to zero for a Defect structure optimization. Select **False** as the value for this parameter.
    9.  The **Constrain position** parameter is a flag to constrain the area in which a defect can be located. Select **False** as the value for this parameter.
    10.  The **Rattle atoms** parameter is a flag indicating if atoms should undergo additional shaking before local minimization. Useful for systems with strong potentials. Select **False** as the value for this parameter.
5.  To enter the values for **Population and individual generation parameters**, click on either the **Next** button on the bottom right of the center pane or on the Population and individual generation parameters tab on the left pane of the web page. The population and individual generation parameters will be displayed on the right pane of the page.
    1.  The **Initial structure size (Angstroms)** parameter specifies the size of volume for generating initial structure. It has a default value of **natoms**<sup>1/3</sup>×r<sub>ab</sub>, where **natoms** is the number of atoms and r<sub>ab</sub> is the **Average atom distance**. Enter **3.0** as the value for this parameter.
    2.  The **Average atom distance (Angstroms)** parameter specifies the average distance between atoms in structure. It has a default value of 2.5\. Enter **1.6** as the value for this parameter.
    3.  The **Custom generation types** parameter indicates the preferred method for generating initial random population. It has a default value of Box. Select **Box** as the value for this parameter.
    4.  The **Size of large box (Angstroms)** parameter specifies the size of box to embed non-periodic cluster structure in for energy evaluation. It has a default value of 500\. Enter **500** as the value for this parameter.
6.  To enter the values for **Post processing parameters**, click on either the **Next** button on the bottom right of the center pane or on the Post processing parameters tab on the left pane of the web page. The post processing parameters will be displayed on the right pane of the page.
    1.  The **Lattice concentration** parameter is a flag indicating if the concentration of lattice atoms versus defect atoms should be calculated. Select **True** as the value for this parameter.
    2.  The **Standard post-processing** parameter is a flag to determine whether to run basic post-processing analysis on optimization output. Select **False** as the value for this parameter.
    3.  The **Genealogy tree** parameter is a flag indicating whether or not to construct a full genealogy tree during post-processing. Select **False** as the value for this parameter.
7.  To enter the values for **Evaluation parameters**, click on either the **Next** button on the bottom right of the center pane or on the Evaluation parameters tab on the left pane of the web page. The evaluation parameters will be displayed on the right pane of the page.
    1.  The **Calculation method** parameter indicates the calculation method to be used for the StructOpt job. There are four options to choose for this parameter i.e. VASP, LAMMPS, MAST_VASP, MAST_LAMMPS. Select **MAST_VASP** as the value for this parameter.
    2.  On selecting either VASP or MAST_VASP as the value for the **Calculation method** parameter, the following parameters are displayed.
        1.  The **VASP Calculator** parameter specifies VASP function to be called. Enter **Vasp()** as the value for this parameter.
        2.  The **VASP parameters** parameter specifies VASP parameters to be used while running VASP. This parameter has a set of default VASP parameter values.Modify the parameter values if needed.
    3.  On selecting either MAST or MAST_LAMMPS as the value for the **Calculation method** parameter, the following parameters are displayed.
        1.  The **Pair style** parameter indicates the Pair style to be used in LAMMPS and MAST_LAMMPS. There are ten options to choose for this parameter i.e. None, BOP, buck, other, tersoff, eam, eam/fs, eam/cd, edip and lj/cut. Select **None** as the value for this parameter.  
            If you select **BOP** as the value for the **Pair style** parameter, the following parameter(s) are displayed.
            1.  The **BOP cutoff** parameter specifies the cutoff distance for BOP potential. It has a default value of 1.If you select **Buckingham** as the value for the **Pair style** parameter, the following parameter(s) are displayed.
            1.  The **Buckingham cutoff** parameter specifies the distance for Buckingham potential. It has a default value of 1.
            2.  The **Buckingham parameters** parameter specifies the cutoff for Buckingham potential. It has a default value of ['* * 100.00 1.5 200.0']If you select **lj/cut** as the value for the **Pair style** parameter, the following parameter(s) are displayed.
            1.  The **Pair coefficient** parameter specifies the pair coefficients for use in LAMMPS or MAST_LAMMPS. It has a default value of * * 1.0 1.0 2.5.If you select **Other** as the value for the **Pair style** parameter, the following parameter(s) are displayed.
            1.  The **Other pair style** parameter specifies the other parameters to be specified for LAMMPS or MAST_LAMMPS.
            2.  The **Pair style name** parameter specifies the pair style name for using in LAMMPS or MAST_LAMMPS 'other' style pair style.
            3.  The **Pair coefficient** parameter specifies the pair coefficients for use in LAMMPS or MAST_LAMMPS.
        2.  The **Potential files** parameter specifies the file for potential to be used in LAMMPS or MAST_LAMMPS.
        3.  The **Keep LAMMPS files** parameter is a flag indicating whether or not LAMMPS data, input, trajectory, and log files should be kept. The default value for this parameter is False.
        4.  The **LAMMPS minimization** parameter specifies the command to be used to locally minimize structure in LAMMPS or MAST_LAMMPS. It has a default value of 1e-8 1e-8 5000 10000.
        5.  The **LAMMPS minimizer style** parameter specifies the command for the style to be used to locally minimize structure in LAMMPS or MAST_LAMMPS. It has a default value of cg  
            nmin_modify line quadratic.
        6.  The **No. of LAMMPS thermo steps** parameter specifies the number of intervals between when states are written to LAMMPS trajectory output. Can be used to prevent large trajectory files. It has a default value of 1.
        7.  The **ASE Minimizer** parameter is a flag to determine whether to use built in ASE minimizers. The default value for this parameter is False.
        8.  The **Maximum solver force** parameter specifies the maximum value of force for ASE structure minimizer. It has a default value of 0.01.
        9.  The **Maximum solver steps** parameter specifies the maximum number of steps to complete for ASE structure minimizer. It has a default value of 2500.
        10.  The **Energy cutoff factor** parameter specifies the value for the energy per atom limit used to identify problems with energy calculations. It has a default value of 10.
8.  To enter the values for **Crossover parameters**, click on either the **Next** button on the bottom right of the center pane or on the Crossover parameters tab on the left pane of the web page. The crossover parameters will be displayed on the right pane of the page.
    1.  The **Crossover probability** parameter specifies the crossover probability. It has a default value of 0.8 for **Genetic algorithm(GA)** value for **Optimizer type** parameter and 0 otherwise. Enter **0.55** as the value for this parameter.
    2.  The **Crossover scheme** parameter indicates the scheme to be used for crossover. There are ten options to choose for this parameter i.e. clustbx, cxtp, cxtpa, cxtpc, newclus, randalloybox, rotct_rand_clus, rotct_rand and rotct. Select **cxtpa** as the value for this parameter.
9.  To enter the values for **Mutation parameters**, click on either the **Next** button on the bottom right of the center pane or on the Mutation parameters tab on the left pane of the web page. The mutation parameters will be displayed on the right pane of the page.
    1.  The **Mutation probability** parameter specifies the probability of a random move occurring. It has a default value of 0.2 for **Genetic algorithm(GA)** value for **Optimizer type** parameter and 1 otherwise. Enter **0.10** as the value for this parameter.
    2.  The **Mutant add** parameter is a flag to determine whether to use built in ASE minimizers. The default value for this parameter is False.
    3.  The **Mutation options** parameter indicates the types of mutations to be performed in optimization. There are thirty three options to choose for this parameter i.e. ase_minimization, atoms_add, atoms_remove, basin_hop_la, basin_hop_permute, basin_hop_ra_atoms, basin_hop_rattle, basin_hop_rotate, cell_relax_lammps, cell_shape, lattice_alteration_crystal, lattice_alteration_group, lattice_alteration_nn, lattice_alteration_rdrd, lattice_alteration_small, lattice_alteration, move_la, permutation_bulk, permutation_crystal, permutation_crystal_multi, permutation, quench, random_replacement, rattle, rotation_geo, rotation, scale_size, swap_int_local, swap_int, swap_vacancy, swap, zp_rotation_fixed, zp_rotation. It has a default value of rotation and lattice_alteration. Select **lattice_alteration, lattice_alteration_rdrd, rotation_geo, rotation and random_replacement'** as the value for this parameter.
10.  To enter the values for **Selection parameters**, click on either the **Next** button on the bottom right of the center pane or on the Selection parameters tab on the left pane of the web page. The selection parameters will be displayed on the right pane of the page.
    1.  The **Selection scheme** parameter indicates the scheme for selecting and pairing individuals for crossover operations. There are fifteen options to choose for this parameter i.e. Best, Cost, Fuss, Fuss1, Fuss1r, Fuss2, Fussf, Fussr, Metropolis, Multitournament, Random_pick, Rank, Tournament, Tournament1 and Tournament2\. It has a default value of tournament2\. Select **fuss** as the value for this parameter.
    2.  The **Natural selection scheme** parameter indicates the scheme for selecting which individuals will survive to the next generation. There are fifteen options to choose for this parameter i.e. Best, Cost, Fuss, Fuss1, Fuss1r, Fuss2, Fussf, Fussr, Metropolis, Multitournament, Random_pick, Rank, Tournament, Tournament1 and Tournament2\. It has a default value of Best. Select **Best** as the value for this parameter.
    3.  The **Fitness scheme** parameter indicates the method for determining the fitness of an individual. There are ten options to choose for this parameter i.e. Chemical potential swap energy, Energy per atom, Enthalpy, Exponential, Formation energy, Silicon bias, Stem cost, Stem cost rotation, Surface Energy and Total energy. It has a default value of Total energy. Select **Total energy** as the value for this parameter.
    4.  The **Convergence scheme** parameter indicates the scheme to determine when the algorithm has converged. There are four options to choose for this parameter i.e. max_gen, gen_rep_min, gen_rep_avg and std. It has a default value of max_gen. Select **gen_rep_min** as the value for this parameter.
    5.  The **Maxgen** parameter specifies the maximum number of generations to be run in an optimization. It has a default value of 5\. Enter **100** as the value for this parameter.
    6.  The **Reqrep** parameter specifies the required number of repetitive fitnesses in order for system to be considered converged. It has a default value of 0\. Enter **20** as the value for this parameter.
    7.  The **Tolerance** parameter specifies the fitness difference for considering a repetitive minimum fitness between generations. It has a default value of 0.001\. Enter **0.15** as the value for this parameter.
    8.  The **Convergence control scheme** parameter indicates the scheme to be used to eliminate duplicate structures and run selection method. There are eleven options to choose for this parameter i.e. adapting, energy_cluster, fingerprint_niche, fitpred_bests, fitpred_new, fitpred, mutation_dups, mutation_dups_energy, mutation_dups_quench, mutation_dups_zp. It has a default value of mutation_dups. Select **mutation_dups** as the value for this parameter.
    9.  The **Demin** parameter specifies the difference in fitness for considering a structure identical to another. It has a default value of 0.005\. Enter **0.1** as the value for this parameter.
11.  To enter the values for **Crystal specific parameters**, click on either the **Next** button on the bottom right of the center pane of the **Population and individual generation parameters** tab or on the Crystal specific parameters tab on the left pane of the web page. The crystal specific parameters will be displayed on the right pane of the page.
    1.  The **Cell shape options** parameter indicates the options for crystal structure optimizer to explore. There are six options to choose for this parameter i.e. Cubic, Hexagonal, Triclinic, Monoclinic, Orthorhombic and Tetragonal. It has a default value of Cubic.
12.  To enter the values for **Defect specific parameters**, click on either the **Next** button on the bottom right of the center pane of the **Population and individual generation parameters** tab or on the Defect specific parameters tab on the left pane of the web page. The defect specific parameters will be displayed on the right pane of the page.
    1.  The **Size factor** parameter specifies the parameter controlling how much of material surrounding region 1 to be included in region 2\. Given as a distance from the center of region 1\. It has a default value of 1.75.
    2.  The **Supercell** parameter specifies the size of supercell of bulk structure provided for defect calculation. It has a default value of (1, 1, 1).
    3.  The **Potential files** parameter specifies the bulk solid file for a Defect structure optimization.
    4.  The **Solid cell size** parameter specifies the cell size for the bulk solid file provided for a Defect structure optimization. It has a default value of (1, 1, 1).
    5.  The **Evaluation of solid** parameter is a flag for determining if a bulk solid should be evaluated for a given system. The default value for this parameter is False.
    6.  The **Find defect scheme** parameter is a flag for determining if the defects should be tracked through the evolution or assumed constant. The default value for this parameter is True.
    7.  The **Track vacancies** parameter is a flag indicating if the optimizer should keep track of vacancies on the lattice. The default value for this parameter is False.
    8.  The **Track swaps** parameter is a flag indicating if the optimizer should keep track of swaps between atom types on the lattice. The default value for this parameter is False.
    9.  The **Alloy** parameter is a flag indicating **********. The default value for this parameter is True.
    10.  The **Random vacancy location start** parameter is a flag indicating if a vacancy should be initialized at a random location or near the center of the bulk solid file. The default value for this parameter is False.
    11.  The **Random location start** parameter is a flag indicating f the defect should be initialized at a random location or at the center of the bulk solid file. The default value for this parameter is False.
13.  After entering all the parameter values, click on the **Run StructOpt** button at the bottom of the **Selection parameters** tab.
14.  If any validation errors are encountered on the entered parameter values, an error message is displayed to the user to modify the required parameter values.
15.  It there are no validation errors, the StructOpt job is successfully submitted and the user is redirected to a confirmation page.

StructOpt website allows user to track their following StructOpt jobs.

1.  Ongoing jobs
2.  Archived jobs
3.  Killed jobs

1.  Click on the **Job tracker** tab.
2.  To list the details of all the **Ongoing jobs**, click on the Ongoing jobs tab on the left pane of the web page. The ongoing jobs details table will be displayed on the right pane of the page.
    1.  The column **Job name** lists the recipe name or the job name of the ongoing StructOpt job.
    2.  The column **Start Time** lists the time at which of the ongoing StructOpt job was submitted.
    3.  The column **Input File** lists the link to download the input file corresponding to the ongoing StructOpt job.
    4.  The column **Kill Job** lists the button to kill/terminate the ongoing StructOpt job.

1.  Click on the **Job tracker** tab.
2.  To list the details of all the **Archived jobs**, click on the Archived jobs tab on the left pane of the web page. The archived jobs details table will be displayed on the right pane of the page.
    1.  The column **Job name** lists the recipe name or the job name of the completed StructOpt job.
    2.  The column **Start Time** lists the time at which of the StructOpt job was submitted.
    3.  The column **End Time** lists the time at which of the StructOpt job was completed.
    4.  The column **Input File** lists the link to download the input file corresponding to the completed StructOpt job.

1.  Click on the **Job tracker** tab.
2.  To list the details of all the **Killed jobs**, click on the Killed jobs tab on the left pane of the web page. The killed jobs details table will be displayed on the right pane of the page.
    1.  The column **Job name** lists the recipe name or the job name of the killed StructOpt job.
    2.  The column **Start Time** lists the time at which of the StructOpt job was submitted.
    3.  The column **End Time** lists the time at which of the StructOpt job was killed.
    4.  The column **Input File** lists the link to download the input file corresponding to the killed StructOpt job.

* * *

<small>File translated from T<sub><span class="small">E</span></sub>X by [T<sub><span class="small">T</span></sub>H](http://hutchinson.belmont.ma.us/tth/), version 4.08.  
On 24 Mar 2016, 23:52.</small>
