import os
import time
import datetime
import math
import sys
from CoolProp.CoolProp import PropsSI

# set timezone
os.environ['TZ'] = 'America/Chicago'

#########################    functions   #################################
def getTransitTime():
    global volumeFlow

    t_transit = total_length * math.pi * pow(diameter,2)/ (4 * volumeFlow)

    ## testing
    t_transit = 2.361
    ##
    return t_transit

def writeLogfile():
    global volumeFlow, Q_tot_in, Q_elec, T_1, T_2

    # times 1000 for recalculation from m^3/s to L/s
    file_obj = open("logfile.txt","a")
    file_obj.write("%s\t %f\t\t %f\t\t %f\t\t %f\t\t %f\n" % (time.strftime("%x %I:%M:%S %p"), volumeFlow * 1000, Q_tot_in, Q_elec, T_1, T_2))
    file_obj.close()

def writeDAfile():
    global  V_pump, V_SCR

    DAfile_list = open("DA_vals.txt").readlines()     # necessary?? (for reading the reference voltage?)
    DAfile_obj = open("DA_vals.txt","w")
    DAfile_obj.write("%f\t %f\t %s\nDAC0\t\t DAC1\t\nV_pump\t\t V_SCR\t\t V_ref_set (5.0 or 3.3)" % (V_pump, V_SCR, DAfile_list[0].split()[2]))
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
   
    global V_SCR, V_pump, T_1, T_2, Q_tot_in, Q_elec, volumeFlow, counter, sum_error_Q
    R_T_ref = 1500
    totalPressure = 101325
    a_stein1 = 1.326e-3
    b_stein1 = 2.903e-4
    c_stein1 = -6.474e-6
    d_stein1 = 3.772e-7
    
    a_stein2 = 1.387e-3 
    b_stein2 = 2.670e-4
    c_stein2 = -3.599e-6
    d_stein2 = 2.619e-7
    
######## update correct values
    a_flow = 2          
    b_flow = 4.5
    
    AD_list = open("AD_vals.txt").readlines()
    V_S             = float(AD_list[0].split()[0])
    V_G             = float(AD_list[0].split()[1])
    V_R1            = float(AD_list[0].split()[2])
    V_R2            = float(AD_list[0].split()[3])
    V_flow          = float(AD_list[0].split()[4])
    V_deltaPressure = float(AD_list[0].split()[5])
    V_Q_elec        = float(AD_list[0].split()[6])
   
    # calculate resistance of thermistors
    R_T1 = (V_S - V_R1) / (V_R1 - V_G) * R_T_ref;
    R_T2 = (V_S - V_R2) / (V_R2 - V_G) * R_T_ref;    

    
    # calculate temperature with Steinhart-Hart equation
    T_1 = 1/(a_stein1 + b_stein1 * math.log(R_T1) + c_stein1 * pow(math.log(R_T1),2) + d_stein1 * pow(math.log(R_T1),3))
    T_2 = 1/(a_stein2 + b_stein2 * math.log(R_T2) + c_stein2 * pow(math.log(R_T2),2) + d_stein2 * pow(math.log(R_T2),3))


    volumeFlow = a_flow + b_flow * V_flow #(m^3/s)
    
    # testing
    volumeFlow = 0.3e-3
    T_2 = T_1 + 2
    delta_p = 200000  
    ########
    
#### add real equations for flow and pressure
              

    T_mean = (T_1 + T_2)/2

    c_p = PropsSI('Cpmass','T', T_mean, 'P', totalPressure,'Water')
    rho = PropsSI('Dmass','T', T_mean, 'P', totalPressure,'Water')

    Q_calc = math.fabs(rho * volumeFlow * c_p * (T_2 - T_1))
    Q_fric = delta_p * volumeFlow

    Q_tot_in = Q_calc + Q_fric

    error_Q = Q_set_ - Q_tot_in

    sum_error_Q = sum_error_Q + error_Q
    counter = counter + 1

# add correct calculation from V_Q_elec
    Q_elec = 3000

    print("counter: %d, sum_error_Q: %d" % (counter, sum_error_Q))
    error_Flow = volumeFlow - volumeFlowrate_set_
    
#### setzwerte aufsummieren und mittelwert bilden?
#### wie PID Regler verwirklichen
  
    
    print("Q_tot_in: %d, Q_set_: %d V_R1: %f, V_R2: %f, T_1: %f, T_2: %f, R_T1: %f, R_T2: %f, T_mean: %f" % (int(Q_tot_in), int(Q_set_), V_R1, V_R2, T_1-273.15, T_2-273.15, R_T1, R_T2, T_mean))

#### Anderungen fuer die Setzwerte sollen langsam gemacht werden,
#### zB 20 Messungen mit 5 Sekunden Unterschied

#### spaeter soll zwischen Q_elec und dem Q_cal  gewichtet werden bis zum 6. Durchlauf,
#### sodass ein langsamer uebergang zur Regelung entsteht
    
    # if the fluid has transitted 6 times, the control method is activated
    # every time the fluid has transitted one time new values are set and written into the DA file.
    currentTime = time.time() - start_clock
# 6*  in die if Abfrage!!!
    if (1 * getTransitTime() < currentTime) & (round(currentTime,0) % round(getTransitTime(),0) == 0):

        mean_error_Q = sum_error_Q / counter
        if abs(mean_error_Q) > 5:

            if mean_error_Q >= 0:
                
                Q_new = Q_calc + 3/4 * mean_error_Q

            else:
                
                Q_new = Q_calc - 3/4 * mean_error_Q


            SCR_supply  = 20000 # 0-20 kW
            plug_supply = 120 # volts
            wireLoops   = 3 # loops of wires around ...
            DArange     = 5 # 0-5 volts
            WTrange     = 10 # 0-10 V

            V_SCR = DArange * Q_new / (SCR_supply / wireLoops)
        else:

            DAfile_list = open("DA_vals.txt").readlines()
            V_SCR = DAfile_list[0].split()[1]                   # V_SCR remains uncontrolled
             
        counter = 0
        sum_error_Q = 0

## control pump

        V_pump = 0

        writeDAfile()
        print(round(currentTime,0),round(getTransitTime(),0))
        print("V_SCR: %f" % V_SCR)






########################################################################
        
# TODO: - Zeitproblem loesen (Abfrage nach jeder sekunde notwendig?)
#       
# EXPLANATION:

# thermistor 1: orange wire
# thermistor 2: purple wire

# analog input channels:
# AD 0:    V_S, source voltage (from PI)
# AD 1:    V_G, ground of PI (AGND)    
# AD 2:    V_R1, thermistor 1
# AD 3:    V_R2, thermistor 2
# AD 4:    V_flow, signal from omega meter (pulse counter)
# AD 5:    V_deltaPressure, signal from pressure sensor  
# AD 6:    V_Q_elec, signal from WT
# AD 7:    

# digital output channels:
# DAC0:     V_pump, signal to VFD
# DAC1:     V_SCR, signal to SCR for the heaters


######################### execution part ################################
command = "sudo ./master"
getParameters()
start_clock = time.time()

maxTimesteps = sys.maxsize
sum_error_Q = 0
counter = 0

# initialize logfile
file_obj = open("logfile.txt","w")
file_obj.write("Date/Time\t\t Flowrate(obs)[L/s]\t Heatinlet(obs)[W]\t Heatinlet(elec)[W]\t Temperature 1(obs)[K]\t Temperature 2(obs)[K]\n")
file_obj.close()

while True:

    currentTimestep = round(time.time() - start_clock,0)

    # run master script: read and set DA, read and write AD
    master_clock1 = time.time()
    os.system(command)
    master_clock2 = time.time()

    # read setpoints, if required
    if currentTimestep % set_time == 0:
        readSetpoints(currentTimestep)

#    AD_list = open("AD_vals.txt").readlines()
#    print("ch0: %f, ch1: %f, ch2: %f, ch3: %f, ch4: %f, ch4: %f, ch6: %f" % (float(AD_list[0].split()[0]),float(AD_list[0].split()[1]),float(AD_list[0].split()[2]),float(AD_list[0].split()[3]),float(AD_list[0].split()[4]),float(AD_list[0].split()[5]),float(AD_list[0].split()[6])))
        
    # calculation of DA_values
    calculationOfDA(Q_set, volumeFlowrate_set)

    # write Logfile
    writeLogfile()

    if(currentTimestep >= maxTimesteps):
	break

total_end = time.time()
print("End of TRT script, total time: %f" % (total_end-start_clock))
