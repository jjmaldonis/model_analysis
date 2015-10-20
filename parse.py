import sys

def parse_t3_catvors(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x for x in content if "exit status" not in x]
    content = [x for x in content if "failed!" not in x]
    #content.pop(0)
    #content.pop(0)
    #content.pop(0)

    step = []
    icot = []
    icoz = []
    icoc = []
    icoa = []
    xtalt = []
    xtalz = []
    xtalc = []
    xtala = []
    mixt = []
    mixz = []
    mixc = []
    mixa = []
    undt = []
    undz = []
    undc = []
    unda = []
    fit = []
    fiz = []
    fic = []
    fia = []

    #type = 'None'
    stepi = 0
    i = 0
    for line in content:
        if "Step:" in line:
            #i = len(step) - 1
            #try:
            #    print("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18} {19} {20}".format(step[i], icot[i], icoz[i], icoc[i], icoa[i], xtalt[i], xtalz[i], xtalc[i], xtala[i], mixt[i], mixz[i], mixc[i], mixa[i], undt[i], undz[i], undc[i], unda[i], fit[i], fiz[i], fic[i], fia[i]))
            #except:
            #    pass
            stepi = int(line.strip().split()[-1][:-4])
            #Step: 1; Model: /home/jjmaldonis/aci/t3/model_update_9318_        1000.xyz
            #stepi = int(line.strip().split()[1][0:-1])
            #print(stepi)
            step.append(stepi)
        elif "Mixed" in line:
            type = "Mixed"
            mixt.append(line.strip().split()[1])
        elif "Undef" in line:
            type = "Undef"
            undt.append(line.strip().split()[1])
        elif "Icosahedra-like" in line:
            type = "Icosahedra-like"
            icot.append(line.strip().split()[1])
        elif "Crystal-like" in line:
            type = "Crystal-like"
            xtalt.append(line.strip().split()[1])
        elif "Full-icosahedra" in line:
            type = "Full-icosahedra"
            fit.append(line.strip().split()[1])
            while(len(fit) < len(icot)): fit.append(0)
        else:
            if type == "Mixed":
                if "Zr" in line:
                    mixz.append(line.strip().split()[1])
                elif "Cu" in line:
                    mixc.append(line.strip().split()[1])
                elif "Al" in line:
                    mixa.append(line.strip().split()[1])
            elif type == "Undef":
                if "Zr" in line:
                    undz.append(line.strip().split()[1])
                elif "Cu" in line:
                    undc.append(line.strip().split()[1])
                elif "Al" in line:
                    unda.append(line.strip().split()[1])
            elif type == "Icosahedra-like":
                if "Zr" in line:
                    icoz.append(line.strip().split()[1])
                elif "Cu" in line:
                    icoc.append(line.strip().split()[1])
                elif "Al" in line:
                    icoa.append(line.strip().split()[1])
            elif type == "Crystal-like":
                if "Zr" in line:
                    xtalz.append(line.strip().split()[1])
                elif "Cu" in line:
                    xtalc.append(line.strip().split()[1])
                elif "Al" in line:
                    xtala.append(line.strip().split()[1])
            elif type == "Full-icosahedra":
                if "Zr" in line:
                    fiz.append(line.strip().split()[1])
                    while(len(fiz) < len(icot)): fiz.append(0)
                elif "Cu" in line:
                    fic.append(line.strip().split()[1])
                    while(len(fic) < len(icot)): fic.append(0)
                elif "Al" in line:
                    fia.append(line.strip().split()[1])
                    while(len(fia) < len(icot)): fia.append(0)

    print("Step, icoT, icoZ, icoC, icoA, xtalT, xtalZ, xtalC, xtalA, mixT, mixZ, mixC, mixA, undT, undZ, undC, undA, fiT, fiZ, fiC, fiA")
    for i in range(0,len(step)):
        print("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18} {19} {20}".format(step[i], icot[i], icoz[i], icoc[i], icoa[i], xtalt[i], xtalz[i], xtalc[i], xtala[i], mixt[i], mixz[i], mixc[i], mixa[i], undt[i], undz[i], undc[i], unda[i], fit[i], fiz[i], fic[i], fia[i]))
                    

if __name__ == "__main__":
    parse_t3_catvors(sys.argv[1])
