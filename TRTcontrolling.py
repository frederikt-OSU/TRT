import os
import time
import datetime
import math
import sys
from collections import deque
from CoolProp.CoolProp import PropsSI

# set timezone
os.environ['TZ'] = 'America/Chicago'

#########################    functions   #################################
def getTransitTime():
    global avg_volumeFlow
    
    t_transit = total_length * math.pi * pow(diameter,2)/ (4 * avg_volumeFlow)

    ## testing
    print(t_transit)
    t_transit = 1.561
    ##
    return t_transit

def writeLogfile():

    global V_SCR, V_pump, V_S, V_G, V_R1, V_R2, V_flow, V_deltaPressure, V_Q_elec
    global avg_volumeFlow, Q_tot_in, Q_elec, Q_fric, T_1, T_2, delta_p

    # times 1000 for recalculation from m^3/s to L/s
    file_obj = open("logfile.txt","a")
    file_obj.write("%s\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t\t %f\n" % (time.strftime("%x %I:%M:%S %p"), avg_volumeFlow * 1000, Q_tot_in, Q_fric, Q_elec, T_1, T_2, delta_p, V_SCR, V_pump, V_S, V_G, V_R1, V_R2, V_flow, V_deltaPressure, V_Q_elec))
    #testing
    global expectedFlow
    file_obj.write("expected flow rate: %f\n" % (expectedFlow * 1000))
    ####
    file_obj.close()

def writeDAfile():
    global  V_pump, V_SCR, V_ref_set

    DAfile_obj = open("DA_vals.txt","w")
    DAfile_obj.write("%f\t %f\t %s\nDAC0\t\t DAC1\t\nV_pump\t\t V_SCR\t\t V_ref_set (5.0 or 3.3)" % (V_pump, V_SCR, V_ref_set))
    DAfile_obj.close()
    
def getParameters():
    global diameter, total_length, set_time, logfile_period, V_ref_set

    file_list = open("Parameters.txt").readlines()
    diameter = float(file_list[1].split()[0])
    total_length = float(file_list[1].split()[1]) + float(file_list[1].split()[2])
    set_time = int(file_list[1].split()[3])
    logfile_period = int(file_list[1].split()[4])
    V_ref_set = float(file_list[1].split()[5])
    
def readSetpoints(currentTs):
    global maxTimesteps, Q_set, volumeFlowrate_set

    file_list = open("Setpoints.txt").readlines() 
    maxTimesteps = int(file_list[len(file_list)-1].split()[0])
    
    for currentLine in range(1,len(file_list)):

        if currentTs < int(file_list[currentLine].split()[0]):

            Q_set = float(file_list[currentLine-1].split()[1])
            volumeFlowrate_set = float(file_list[currentLine-1].split()[2])
            break

def averageVariables():
    global avg_volumeFlow, sum_volumeFlow, counterAverageFunc
    global avg_Q_elec, sum_Q_elec, avg_T1, avg_T2, sum_T1, sum_T2, sum_delta_p, avg_delta_p
    
    avg_volumeFlow = sum_volumeFlow / counterAverageFunc
    avg_Q_elec = sum_Q_elec / counterAverageFunc
    avg_T1 = sum_T1 / counterAverageFunc
    avg_T2 = sum_T2 / counterAverageFunc
    avg_delta_p = sum_delta_p / counterAverageFunc
    
    # set the counter and the summation to zero in order to calculate the new average
    sum_volumeFlow = 0
    sum_Q_elec = 0
    sum_T1 = 0
    sum_T2 = 0
    sum_delta_p = 0 
    counterAverageFunc = 0

def calculationOfDA(Q_set_,volumeFlowrate_set_):
   
    global V_SCR, V_pump, V_S, V_G, V_R1, V_R2, V_flow, V_deltaPressure, V_Q_elec
    global T_1, T_2, Q_tot_in, Q_elec, Q_fric, avg_volumeFlow, avg_Q_elec, delta_p
    global counter, counterAverageFunc, sum_Q_elec, sum_error_Q, sum_volumeFlow, average_switch, secondPoint_bool, secondPoint_switch
    global sum_T1, sum_T2, avg_T1, avg_T2, sum_delta_p, avg_delta_p
    
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
            
    SCR_maxAmps  = 30               # max ampere from the SCR to the system
    SCR_voltSupply = 200            # volt supply from the SCR to the system
    DAlowSCR    = 1                 # figured out with trail
    DAhighSCR   = 4.6               # figured out with trail
    ADrange     = 5                 # 0-5 V
    pressureRange = 25 * 6894.757   # range of pressure sensor (0-25 psi) in pascal (1 psi = 6894.757 Pa)
    P_range     = 5                 # 0-5 V
    
    DArangeSCR  = DAhighSCR - DAlowSCR    
    Q_max = SCR_maxAmps * SCR_voltSupply  # max possible heat inlet of heaters
    
######## update correct values
    a_flow = 0          
    b_flow = 1.44
    
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

    #volumeFlow = (a_flow + b_flow * V_flow ) * 0.000063090197
    volumeFlow = (a_flow + b_flow * (V_flow - V_G)) * 0.000063090197 #(m^3/s) (1 gpm = 0.000063090197 m^3/s) 

    #delta_p = pressureRange * V_deltaPressure / P_range    
    delta_p = pressureRange * (V_deltaPressure - V_G) / P_range

    # calculate electrical heat inlet    
    Q_elec = V_Q_elec  * 1000                   # the output of the watt transducer is directly proportional to the electrical heat inlet
# test
    Q_elec = 2.94 * 1000 # measured with multimeter
    #Q_elec = (V_Q_elec - V_G) * 1000
    
    # in order to get a start value for the volume flow rate, the first calculated value is passed
    if average_switch:
        #time.sleep(2)                           # wait for the system to start flowing (only necessary, if system is switched directly before starting the program
        avg_volumeFlow = volumeFlow
        avg_Q_elec = Q_elec
        avg_T1 = T_1
        avg_T2 = T_2
        avg_delta_p = delta_p
        average_switch = 0

# 20 sec        
    # every 20 seconds the readed flow rate is averaged in order to counteract to fluctuations of the pulse meter    
    currentTime = time.time() - start_clock
    if round(currentTime,0) % 8 == 0:
        averageVariables()

    # calculate parameters cp and rho
    T_mean = (avg_T1 + avg_T2)/2
    c_p = PropsSI('Cpmass','T', T_mean, 'P', totalPressure,'Water')
    rho = PropsSI('Dmass','T', T_mean, 'P', totalPressure,'Water')     

        
    # testing
    global expectedFlow
    expectedFlow = avg_Q_elec / (rho * c_p * math.fabs(avg_T2 - avg_T1))
    ############
    
    # calculate calorimetric and frictional heat inlet
    Q_calc = math.fabs(rho * avg_volumeFlow * c_p * (T_2 - T_1))
    Q_fric = delta_p * avg_volumeFlow
    Q_tot_in = Q_calc + Q_fric
    
    error_Q = Q_set_ - Q_tot_in

    # for every run of the main loop, the volumeflow, Q_elec, temperatures, pressure and the error in Q are summed in order to
    # average them afterwards with the "counter"
    sum_Q_elec = sum_Q_elec + Q_elec
    sum_volumeFlow = sum_volumeFlow + volumeFlow
    sum_T1 = sum_T1 + T_1
    sum_T2 = sum_T2 + T_2
    sum_delta_p = sum_delta_p + delta_p
    sum_error_Q = sum_error_Q + error_Q
    counterAverageFunc = counterAverageFunc + 1         # counter for averaging of flow rate, temperatures and electrical heat inlet  
    counter = counter + 1                       # counter for averaging error in Q
    
    currentTime = time.time() - start_clock
    transitTime = getTransitTime()

    # every time the fluid has transitted one time, new values are set and written into the DA file.
    if (round(currentTime,0) % round(transitTime,0) == 0) or currentTime < 5:

        DAfile_list = open("DA_vals.txt").readlines()
        V_SCR_old = float(DAfile_list[0].split()[1])
        v_scr.append(V_SCR_old)                 # the current V_SCR results 
        f_v_scr.append(Q_tot_in)                # in a certain heat inlet
        
        # if the fluid has transitted 6 times, the control method is activated
# 6 * transitTime !!
        if 6 * transitTime < currentTime:
            avg_error_Q = sum_error_Q / counter

# 6.5 
            if 6.5 * transitTime < currentTime and secondPoint_switch == 1: # actual control
                # secant method, modifing to a "root finding problem"
                f_n_minus1 = float(f_v_scr[1]) - Q_set_
                f_n_minus2 = float(f_v_scr[0]) - Q_set_
                x_n_minus1 = float(v_scr[1])
                x_n_minus2 = float(v_scr[0])

                #print(v_scr,f_v_scr)
                #print(f_n_minus1, x_n_minus1, f_n_minus2, x_n_minus2)
                V_SCR =  x_n_minus1 - f_n_minus1 * (x_n_minus1 - x_n_minus2) / (f_n_minus1 - f_n_minus2)
# 6.5                    
            if 6.5 * transitTime < currentTime and secondPoint_bool == 1: # in order to get a second grid point (only needed one time)
                #print("fourth if statement")
                if avg_error_Q >= 0:                
                    Q_new = Q_calc + 3/4 * avg_error_Q
                else:                
                    Q_new = Q_calc - 3/4 * abs(avg_error_Q)
                V_SCR = DArangeSCR * Q_new / Q_max + DAlowSCR
                secondPoint_bool = 0                 # deactivate getting a second grid point
                secondPoint_switch = 1               # activate control

            counter = 0
            sum_error_Q = 0

        else:                                        # V_SCR in the first 6 transit loops          
            V_SCR = DArangeSCR * Q_set_ / Q_max + DAlowSCR
            
    ## control of the pump is made with VFD so far
        V_pump = 0 
        # testing
        V_SCR = 2.5 # 1.6 -> 1000 W
        
        writeDAfile()

########################################################################
        
# TODO: - Zeitproblem loesen (Abfrage nach jeder sekunde notwendig?)
# #### spaeter soll zwischen Q_elec und dem Q_cal  gewichtet werden bis zum 6. Durchlauf,
# #### sodass ein langsamer uebergang zur Regelung entsteht      
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
f_v_scr = deque(maxlen = 2)
v_scr   = deque(maxlen = 2)
secondPoint_bool = 1
secondPoint_switch = 0

maxTimesteps = sys.maxsize
sum_error_Q = 0
sum_volumeFlow = 0
sum_Q_elec = 0
sum_T1 = 0
sum_T2 = 0
sum_delta_p = 0
counter = 0
counterAverageFunc = 0
average_switch = 1

# initialize logfile
file_obj = open("logfile.txt","w")
file_obj.write("Date/Time\t\t Flow rate(obs)[L/s]\t Heat inlet(tot)[W]\t Heat inlet(fric)[W]\t Heat inlet(elec)[W]\t Temperature 1(obs)[K]\t Temperature 2(obs)[K]\t Pressure drop[Pa]\t V_SCR\t\t V_pump\t\t V_S\t\t V_G\t\t V_R1\t\t V_R2\t\t V_flow\t\t V_deltaPressure\t V_Q_elec\n")
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
        
    # calculation of DA_values
    calculationOfDA(Q_set, volumeFlowrate_set)

    AD_list = open("AD_vals.txt").readlines()
    DA_list = open("DA_vals.txt").readlines()
    print("V_S(ch0): %f, V_G(ch1): %f, V_R1(ch2): %f, V_R2(ch3): %f, V_flow(ch4): %f, V_pres(ch5): %f, noSig-V_Q_elec(ch6): %f, noSig: %f" % (float(AD_list[0].split()[0]),float(AD_list[0].split()[1]),float(AD_list[0].split()[2]),float(AD_list[0].split()[3]),float(AD_list[0].split()[4]),float(AD_list[0].split()[5]),float(AD_list[0].split()[6]),float(AD_list[0].split()[7])))
    print("DA0: %f, DA1: %f" % (float(DA_list[0].split()[0]), float(DA_list[0].split()[1])))
    
    # write Logfile
    if currentTimestep % logfile_period == 0:
        writeLogfile()

    if(currentTimestep >= maxTimesteps):
	break

total_end = time.time()
print("End of TRT script, total time: %f" % (total_end-start_clock))
