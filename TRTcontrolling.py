import os
import time
import math
import sys
from CoolProp.CoolProp import PropsSI

# Explanation:
# V_SCR           --> DAC1
# V_pump          --> DAC0 
# V_T_ref1        --> AD0
# V_T1            --> AD1
# V_T_ref2        -->
# V_T2            -->
# V_flow          -->
# V_deltaPressure -->
# V_Q_elec        -->


#########################    functions   #################################
def getTransitTime():
    global volumeFlow

    t_transit = total_length * math.pi * pow(diameter,2)/ (4 * volumeFlow)
    return t_transit

# write Logfile and DA values
def writeLogAndDAfile():
    global volumeFlow, V_pump, V_SCR, Q_tot_in, T_1, T_2
    
    file_obj = open("logfile.txt","a")
    file_obj.write("%s\t %f\t\t %f\t\t %f\t\t %f\n" % (time.strftime("%c"), volumeFlow, Q_tot_in, T_1, T_2))
    file_obj.close()

    DAfile_list = open("DA_vals.txt").readlines()     # necessary?? (for reading the reference voltage?)
    DAfile_obj = open("DA_vals.txt","w")
    DAfile_obj.write("%f\t %f\t %s\nDAC0\t\t DAC1\t\nV_SCR\t\t V_pump\t\t V_ref_set (5.0 or 3.3)" % (V_pump, V_SCR, DAfile_list[0].split()[2]))
    DAfile_obj.close()

    
def getParameters():
    global diameter, total_length, set_time

    file_list = open("Parameters.txt").readlines()
    diameter = float(file_list[1].split()[0])
    total_length = float(file_list[1].split()[1]) + float(file_list[1].split()[2])
    set_time = int(file_list[1].split()[3])

    
def readSetpoints(currentTs):
    global maxTimesteps, Q_set, volumeFlowrate_set

    file_list = open("Setpoints.txt").readlines() 
    maxTimesteps = int(file_list[len(file_list)-1].split()[0])
    
    for currentLine in range(1,len(file_list)):

        if currentTs < int(file_list[currentLine].split()[0]):

            Q_set = float(file_list[currentLine-1].split()[1])
            volumeFlowrate_set = float(file_list[currentLine-1].split()[2])
            break

def calculationOfDA(Q_set_,volumeFlowrate_set_):
    global V_SCR, V_pump, T_1, T_2, Q_tot_in, volumeFlow, counter
    counter =+ 1
    R_T_ref1 = 5000
    R_T_ref2 = 3000
    totalPressure = 101325
    a_stein = 2.1e-3
    b_stein = 0.7e-4
    c_stein = 6.5e-7
    d_stein = 4e-9
    a_flow = 2
    b_flow = 4.5
    
    AD_list = open("AD_vals.txt").readlines()
    V_T_ref1        = float(AD_list[0].split()[0])
    V_T1            = float(AD_list[0].split()[1])
    V_T_ref2        = float(AD_list[0].split()[2])
    V_T2            = float(AD_list[0].split()[3])
    V_flow          = float(AD_list[0].split()[4])
    V_deltaPressure = float(AD_list[0].split()[5])
    V_Q_elec        = float(AD_list[0].split()[6])

    # calculate resistance of thermistors
    R_T1 = V_T1/ (V_T_ref1/ R_T_ref1);
    R_T2 = V_T2/ (V_T_ref2/ R_T_ref2);    

    
    # calculate temperature with Steinhart-Hart equation
    T_1 = 1/(a_stein + b_stein * math.log(R_T1) + c_stein * pow(math.log(R_T1),2) + d_stein * pow(math.log(R_T1),3))
    T_2 = 1/(a_stein + b_stein * math.log(R_T2) + c_stein * pow(math.log(R_T2),2) + d_stein * pow(math.log(R_T2),3))


    volumeFlow = a_flow + b_flow * V_flow
#### add real equation
    delta_p = V_deltaPressure                


    T_mean = (T_1 + T_2)/2

    c_p = PropsSI('Cpmass','T', T_mean, 'P', totalPressure,'Water')
    rho = PropsSI('Dmass','T', T_mean, 'P', totalPressure,'Water')

    Q_tot_in = math.fabs(rho * volumeFlow * c_p * (T_2 - T_1)) + delta_p * volumeFlow

    error_Q = Q_tot_in - Q_set_
    error_Flow = volumeFlow - volumeFlowrate_set_
    
#### setzwerte aufsummieren und mittelwert bilden?
#### wie PID Regler verwirklichen

    
    
    #print("V_T_ref1: %f, V_T1: %f V_T_ref2: %f, V_T2: %f, V_flow: %f, V_deltaPressure: %f, V_Q_elec: %f, c_p: %f, T_mean: %f" % (V_T_ref1, V_T1, V_T_ref2, V_T2, V_flow, V_deltaPressure, V_Q_elec, c_p, T_mean))

#### Anderungen fuer die Setzwerte sollen langsam gemacht werden,
#### zB 20 Messungen mit 5 Sekunden Unterschied

#### spaeter soll zwischen Q_elec und dem Q_cal  gewichtet werden bis zum 6. Durchlauf,
#### sodass ein langsamer uebergang zur Regelung entsteht
    
    # if the fluid has tranistted 6 times, the control method is activated
    if 6 * getTransitTime() < (time.time() - start_clock) :
        V_SCR = 3
        V_pump = 1
    else:
#### hier muessen die ungeregelten werte hin        
        V_SCR = 0
        V_pump = 0
    
    print(6 * getTransitTime(), V_SCR, V_pump)




# TODO: - Zeitproblem loesen (Abfrage nach jeder sekunde notwendig?)
#       

######################### execution part ################################
command = "sudo ./master"
getParameters()

start_clock = time.time()

#tmpTime = -1
maxTimesteps = sys.maxsize

# initialize logfile
file_obj = open("logfile.txt","w")
file_obj.write("Date/Time\t\t\t Flowrate(obs)[L/s]\t Heatinlet(obs)[W]\t Temperature 1(obs)[K]\t Temperature 2(obs)[K]\n")
file_obj.close()

while True:
    end_clock = time.time()
    currentTimestep = round(end_clock - start_clock,0)
    #print(currentTimestep,tmpTime+1)
	
#    if currentTimestep == (tmpTime + 1):
#    if True:		
        #tmpTime = currentTimestep
       # print(currentTimestep, tmpTime +1)
    # run master script: read and set DA, read and write AD
    master_clock1 = time.time()
    os.system(command)
    master_clock2 = time.time()
    # print("time masterscript: %f\n" % (master_clock2-master_clock1))

    # read setpoints, if required
    if currentTimestep%set_time == 0:
        readSetpoints(currentTimestep)
        print("max Timesteps: %d, Q_set %f, flow_set: %f \n" % (maxTimesteps, Q_set, volumeFlowrate_set))
        
    # calculation of DA_values
    calculationOfDA(Q_set, volumeFlowrate_set)

    # write Logfile and DA values
    writeLogAndDAfile()    

    if(currentTimestep >= maxTimesteps):
	break

total_end = time.time()
print("End of TRT script, total time: %f" % (total_end-start_clock))
