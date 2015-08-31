#!/usr/bin/env python3

##### make 1. ".xyz" movie files with a given range and nlocks #####
#####      2. ".xyz" single snapshot at specified step         #####
#####      3. ".ltcp" single snapshot at specified step        #####
#####                                      By Sky on July 2015 #####

import sys
import os.path

if (len(sys.argv) is not 2) and (len(sys.argv) is not 4):
    print("Make .xyz or .ltcp files from history files")
    exit("Usage: %s <t0.lxyz/t0.ltcp> <history.def> <history.sol>" % sys.argv[0])

### Global variables ###
N_atoms, N_vac= (0, 0)
defect, solute= ([], [])
x, y, z= ([], [], [])
is1file, ismovie, isltcp= (False, False, False)
step_min, step_max, bk_steps, step_out= (0, 0, 0, 0)
### Global variables ###

### Functions ###
def writedata(timestep_, time_):
    global defect, solute, x, y, z
    print(N_vac, file=VAC);             print("time:", timestep_, time_, file=VAC);
    print(len(defect)-N_vac, file=ITL); print("time:", timestep_, time_, file=ITL);
    for tp, ltcp in defect:
        if tp == 0: print("0", x[ltcp], y[ltcp], z[ltcp], file=VAC)
        else:       print("2", x[ltcp], y[ltcp], z[ltcp], file=ITL)
    
    print(len(solute), file=A02);       print("time:", timestep_, time_, file=A02);
    for ltcp in solute:
        print("-1", x[ltcp], y[ltcp], z[ltcp], file=A02)

def writeltcp(timestep_, time_):
    global defect, solute, x, y, z
    states= [1] * N_atoms
    for tp, ltcp in defect:
        states[ltcp]= tp
    for ltcp in solute:
        states[ltcp]= -1

    print(N_atoms, file=LTC); print("T:", timestep_, time_, file=LTC);
    for s, a, b, c in zip(states, x, y, z):
        if (s == 1) or (s == -1): print(s, a, b, c, file=LTC)
        elif s == 0:              print(s, a, b, c, "0 0 0", file=LTC)
        else:                     print(s, a, b, c, "0 0 0 0 0", file=LTC)

count_step= -1
count_out= 0
isout= False
def check_exe(timestep_):
# 0: X output, 1: output, 2: exit
    global count_step, count_out, isout

    if ismovie:
        if timestep_ > step_max: 
            return 2 

        if timestep_ >= step_min:
            count_step += 1
            
            if count_step%bk_steps == 0:
                count_out += 1
                if count_out>1000: exit("no more than 1000 steps can be dumped to history") # check 
                sys.stdout.write("\r%d (timesteps)/%d (shots)" % (timestep_, count_out))
                isout= True
                return 1

    else:
        if timestep_ == step_out: 
            isout= True
            return 1
        if timestep_ >  step_out: 
            return 2

    return 0 # if doesn't match any conditions above
### Functions ###

if len(sys.argv) == 2: # Extract arguments 
    (dump, name_t0)= sys.argv
    is1file= True
    ismovie= False
    isltcp=  False
else:                   
    (dump, name_t0, name_def, name_sol)= sys.argv
    is1file= False
    if ".def" not in name_def or ".sol" not in name_sol: 
        exit("name %s or %s is not .def or .sol" % (name_def, name_sol))

if is1file: # Input information variables
    print("Input only one file.\nAre you sure to output only one snapshot? (yes/no)")
    yorn= str(input())
    if yorn == "yes":   pass
    elif yorn == "no":  exit("Good bye! (DO NOTHING)")
    else:               exit("ERROR: Please specify yes or no")
else:
    print("Output a movie or a snapshot? (1:movie)")
    if int(input()) == 1:
        print("<movie mode>(xyz)")
        ismovie= True; isltcp= False
        print("insert min step"); step_min= int(input());
        print("insert max step"); step_max= int(input());
        print("insert bk steps"); bk_steps= int(input());
    else:
        ismovie= False
        print("Output for calculations(.ltcp) or VMD(.xyz)? (1:ltcp)")
        if int(input()) == 1:   isltcp= True;  print("<sinlge shot>(ltcp)")
        else:                   isltcp= False; print("<sinlge shot>(xyz)")
        print("insert the output timestep"); step_out= int(input())

print("start calculating...");

if(isltcp): # Open output files
    if ".ltcp" not in name_t0: exit("name %s is not .ltcp" % name_t0)
    LTC= open(str(step_out)+".ltcp", "w")
else:
    if ".xyz"  not in name_t0: exit("name %s is not .xyz" % name_t0)
    VAC= open("vac.xyz", "w")
    ITL= open("itl.xyz", "w")
    A02= open("a02.xyz", "w")

if not os.path.isfile(name_t0): exit("%s doesn\'t exist!" % name_t0)
with open(name_t0, "r") as T0: # Reading t0 file
    line= T0.readline()
    if line.strip().isdigit(): N_atoms= int(line)
    next(T0) # skip line 2
    defect= []
    solute= []
    x, y, z= ([0]*N_atoms, [0]*N_atoms, [0]*N_atoms)
    for i in range(0, N_atoms):
#! some unknown problem
        tp, x[i], y[i], z[i], *others = T0.readline().split()
        if int(tp) == 1:    pass
        elif int(tp) == -1: solute.append(i)
        else:               defect.append([tp, i])

if is1file:
    writedata(0, 0)
    exit("Job Completed for %s (SINGLE SHOT)"  % name_t0)  # complete
else:
    if check_exe(0) == 1:   writedata(0, 0)
    else:                   pass

    if not os.path.isfile(name_def): exit("%s doesn\'t exist!" % name_def)
    if not os.path.isfile(name_sol): exit("%s doesn\'t exist!" % name_sol)
    with open(name_def, "r") as DEF, open(name_sol, "r") as SOL: # Reading history files
        while True:

            ########## read 1st line ##########
            dline= DEF.readline()
            sline= SOL.readline()
            if (not dline) or (not sline): break # check EOF
            if dline.strip().isdigit() and sline.strip().isdigit():
                dnl= int(dline) # number of lines of def
                snl= int(sline)
            else:
                exit("(read line1) supposed to be number of atoms in def && sol file:\n", dline, "\n", sline)

            ########## read 2nd line ##########
            (c_t1, timestep, time   )= DEF.readline().split(); timestep= int(timestep); time= float(time)
            (c_t2, ts_check, t_check)= SOL.readline().split(); ts_check= int(ts_check); t_check= float(t_check)
            if (c_t1 != "T:") or (c_t2 != "T:"): exit("(read line2) supposed to be T: at beginning: %s %s" % (c_t1, c_t2)) # check 
            if abs(timestep-ts_check)>1e-8 or abs(time-t_check)>1e-8: exit("(read line2) time diff:\ndef: %s %s\nsol: %s %s" % (timestep, time, ts_check, t_check)) # check

            ######### read coordinates ########
            defect= []
            solute= []
            N_vac= 0
            cal_status= check_exe(timestep)
            if cal_status == 0: # X output
                for i in range(0, dnl): next(DEF)
                for i in range(0, snl): next(SOL)
            elif cal_status == 1: # output
                for i in range(0, dnl):
                    (tp, ltcp, ix, iy, iz)= DEF.readline().split()
                    if int(tp) == 0: N_vac += 1
                    defect.append([int(tp), int(ltcp)])
                for i in range(0, snl):
                    ltcp= int(SOL.readline())
                    solute.append(ltcp)
                if isltcp:  writeltcp(timestep, time)
                else:       writedata(timestep, time)
            elif cal_status == 2: # exit
                break

    if not isout: exit("ERROR: DID NOT OUTPUT ANY DATA ")
    elif ismovie: exit("\nJob Completed for %d snapshots (MOVIE MODE)" % count_out) # complete
    else:         exit(  "Job Completed at timestep %d (SINGLE SHOT)"  % step_out)  # complete
