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

# Function calculates the transit time for the fluid.
# parameters: none
# return: calculated transit time.
def getTransitTime():
    global volumeFlow
    
    t_transit = total_length * math.pi * pow(diameter,2)/ (4 * volumeFlow)

    return t_transit

# Function adds the values of the listed global variables to the logfile. 
# parameters: none
def writeLogfile():
    global V_SCR, V_pump, avgV_S, avgV_G, avgV_R1, avgV_R2, avgV_flow, avgV_deltaPressure, avgV_Q_elec
    global volumeFlow, Q_tot_in, Q_elec, Q_fric, T_1, T_2, delta_p

    # times 1000 for recalculation from m^3/s to L/s
    file_obj = open("logfile.txt","a")
    file_obj.write("%s\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t\t %f\n" % (time.strftime("%x %I:%M:%S %p"), volumeFlow * 1000, Q_tot_in, Q_fric, Q_elec, T_1, T_2, delta_p, V_SCR, V_pump, avgV_S, avgV_G, avgV_R1, avgV_R2, avgV_flow, avgV_deltaPressure, avgV_Q_elec))
    file_obj.close()

# Function writes the current DA file with the listed global variables for the master program.  
# parameters: none
def writeDAfile():
    global  V_pump, V_SCR, V_ref_set

    DAfile_obj = open("DA_vals.txt","w")
    DAfile_obj.write("%f\t %f\t %s\nDAC0\t\t DAC1\t\nV_pump\t\t V_SCR\t\t V_ref_set (5.0 or 3.3)" % (V_pump, V_SCR, V_ref_set))
    DAfile_obj.close()

# Function takes the parameters provided from the Parameters.txt-file.  
# parameters: none
# set_time --> time period after which the set file is read repeatedly. 
def getParameters():
    global diameter, total_length, set_time, logfile_period, V_ref_set

    file_list = open("Parameters.txt").readlines()
    diameter = float(file_list[1].split()[0])
    total_length = float(file_list[1].split()[1]) + float(file_list[1].split()[2])
    set_time = int(file_list[1].split()[3]) * 60                        # times 60, because program is working in seconds
    logfile_period = int(file_list[1].split()[4])
    V_ref_set = float(file_list[1].split()[5])

# Function takes the updated setPoint of heat inlet and volume flowrate from the Setpoints.txt-file.
# parameters: currentTs --> setPointStep, timestep of the whole measurement
#             Q_setOld --> current set point of heat inlet
#             volumeFlowrate_setOld --> current set point of flowrate   
def readSetpoints(overallTs, Q_setOld, volumeFlowrate_setOld):
    global maxTimesteps, Q_set, volumeFlowrate_set, start_clock, secondPoint_bool, secondPoint_switch
    
    file_list = open("Setpoints.txt").readlines() 
    maxTimesteps = int(file_list[len(file_list)-1].split()[0]) * 60     # times 60, because program is working in seconds

    # Every line of the file is read and checked, if the line provides the setpoint for the current timestep. 
    for currentLine in range(1,len(file_list)):
        if overallTs < int(file_list[currentLine].split()[0]) * 60:     # times 60, because program is working in seconds

            Q_set = float(file_list[currentLine-1].split()[1])
            volumeFlowrate_set = float(file_list[currentLine-1].split()[2])
            
            # start_clock has to be started again in order to wait another 6 transit times before
            # switching to the controlling, if a new (different) Q_set or volumeFlowrate_set is set.
            # Furthermore, switch variables has to be reset in order to get another second gridpoint.            
            if Q_set != Q_setOld or volumeFlowrate_set != volumeFlowrate_setOld:
                print("--------------------- Q_set or flow_set switched -------------------")
                start_clock = time.time()
                secondPoint_bool = 1
                secondPoint_switch = 0
                
            break

# Function averages all seven measured voltage variables.
def averageVariables():
    global counterAverageFunc
    global avgV_S, avgV_G, avgV_R1, avgV_R2, avgV_flow, avgV_deltaPressure, avgV_Q_elec
    global sumV_S, sumV_G, sumV_R1, sumV_R2, sumV_flow, sumV_deltaPressure, sumV_Q_elec

    avgV_S = sumV_S / counterAverageFunc
    avgV_G = sumV_G / counterAverageFunc
    avgV_R1 = sumV_R1 / counterAverageFunc
    avgV_R2 = sumV_R2 / counterAverageFunc
    avgV_flow = sumV_flow / counterAverageFunc
    avgV_deltaPressure = sumV_deltaPressure / counterAverageFunc
    avgV_Q_elec = sumV_Q_elec / counterAverageFunc
    
    # set the counter and the summation to zero in order to calculate the new average    
    sumV_S = 0
    sumV_G = 0
    sumV_R1 = 0
    sumV_R2 = 0
    sumV_flow = 0
    sumV_deltaPressure = 0
    sumV_Q_elec = 0
    counterAverageFunc = 0

# Function calculates the temperatures, flow rate and heat inlets based on the measured voltages from the AD_vals.txt- file,
# which is provided by the master program. Afterwards the controlling method provides the DA values. 
# parameters: Q_set_ --> current setpoint for heat inlet
#             volumeFlowrate_set_ --> current setpoint for flow rate
def calculationOfDA(Q_set_,volumeFlowrate_set_):   
    global V_SCR, V_pump, V_S, V_G, V_R1, V_R2, V_flow, V_deltaPressure, V_Q_elec
    global T_1, T_2, Q_tot_in, Q_elec, Q_fric, volumeFlow, delta_p
    global counterAverageFunc, average_switch, secondPoint_bool, secondPoint_switch, currentTimestep

    global avgV_S, avgV_G, avgV_R1, avgV_R2, avgV_flow, avgV_deltaPressure, avgV_Q_elec
    global sumV_S, sumV_G, sumV_R1, sumV_R2, sumV_flow, sumV_deltaPressure, sumV_Q_elec
    
    R_T_ref = 1500                              # unit Ohm
    totalPressure = 101325                      # unit Pa

    # Steinhart-Hart coefficients provided by thermistor calibration 
    a_stein1 = 1.326e-3
    b_stein1 = 2.903e-4
    c_stein1 = -6.474e-6
    d_stein1 = 3.772e-7
    
    a_stein2 = 1.387e-3 
    b_stein2 = 2.670e-4
    c_stein2 = -3.599e-6
    d_stein2 = 2.619e-7

    # flow meter coefficients provided by flow meter calibration        
    b_flow = 0.0001004007
    
    SCR_maxAmps  = 30                           # max ampere from the SCR to the system
    SCR_voltSupply = 200                        # volt supply from the SCR to the system
    DAlowSCR    = 1                             # figured out with measurement
    DAhighSCR   = 4.6                           # figured out with measurement
    pressureRange = 25 * 6894.757               # range of pressure sensor (0-25 psi) in pascal (1 psi = 6894.757 Pa)
    P_range     = 5                             # range for voltage of pressure sensor: 0-5 V
    V_SCR_mod = 0.00001                         # value to manipulate (modify) the V_SCR signal in order to reactivate the secant method
     
    DArangeSCR  = DAhighSCR - DAlowSCR    
    Q_max = SCR_maxAmps * SCR_voltSupply        # max possible heat inlet of heaters

    # Read voltages from the AD_vales - file, which is provided by the master program.
    AD_list = open("AD_vals.txt").readlines()
    V_S             = float(AD_list[0].split()[0])
    V_G             = float(AD_list[0].split()[1])
    V_R1            = float(AD_list[0].split()[2])
    V_R2            = float(AD_list[0].split()[3])
    V_flow          = float(AD_list[0].split()[4])
    V_deltaPressure = float(AD_list[0].split()[5])
    V_Q_elec        = float(AD_list[0].split()[6])

    # For every run of the main loop, voltage values are summed in order to
    # average them afterwards with the counterAverageFunc.
    sumV_S = sumV_S + V_S
    sumV_G = sumV_G + V_G
    sumV_R1 = sumV_R1 + V_R1
    sumV_R2 = sumV_R2 + V_R2
    sumV_flow = sumV_flow + V_flow
    sumV_deltaPressure = sumV_deltaPressure + V_deltaPressure
    sumV_Q_elec = sumV_Q_elec + V_Q_elec
    counterAverageFunc = counterAverageFunc + 1  # counter for averaging of flow rate, temperatures and electrical heat inlet
    
    # In order to get a start value for the volume flow rate, the first calculated value is passed.
    # Afterwards this if statement is "switched off".
    if average_switch:
        #time.sleep(2)                           # wait for the system to start flowing (only necessary, if program starts before the pump is started.
        avgV_S = V_S
        avgV_G = V_G
        avgV_R1 = V_R1
        avgV_R2 = V_R2
        avgV_flow = V_flow
        avgV_deltaPressure = V_deltaPressure
        avgV_Q_elec = V_Q_elec
        average_switch = 0
##
    old_transitTime = sys.maxsize                # initialized value which is higher as every possible transitTimes
    if currentTimestep > 1:
        old_transitTime = getTransitTime()
##    
    # Before the logfile is written all voltage inputs are averaged.    
    if currentTimestep % logfile_period  == 0:
        print("------------------- average variables ----------------------")
        averageVariables()
   
    # calculate resistance of thermistors
    R_T1 = (avgV_S - avgV_R1) / (avgV_R1 - avgV_G) * R_T_ref;
    R_T2 = (avgV_S - avgV_R2) / (avgV_R2 - avgV_G) * R_T_ref;    
    
    # calculate temperature with Steinhart-Hart equation
    T_1 = 1/(a_stein1 + b_stein1 * math.log(R_T1) + c_stein1 * pow(math.log(R_T1),2) + d_stein1 * pow(math.log(R_T1),3))
    T_2 = 1/(a_stein2 + b_stein2 * math.log(R_T2) + c_stein2 * pow(math.log(R_T2),2) + d_stein2 * pow(math.log(R_T2),3))

    # calculate flow rate
    volumeFlow = b_flow * (avgV_flow - avgV_G)   # (m^3/s) 

    # calculate pressure drop over the borehole   
    delta_p = pressureRange * (avgV_deltaPressure - avgV_G) / P_range

    # calculate electrical heat inlet    
    Q_elec = avgV_Q_elec  * 1000                 # the output of the watt transducer is directly proportional to the electrical heat inlet

# has to be figured out
    Q_elec = 0
    #Q_elec = V_Q_elec * 1000

    # calculate parameters cp and rho
    T_mean = (T_1 + T_2)/2
    c_p = PropsSI('Cpmass','T', T_mean, 'P', totalPressure,'Water')
    rho = PropsSI('Dmass','T', T_mean, 'P', totalPressure,'Water')     
    
    # calculate calorimetric and frictional heat inlet
    Q_calc = math.fabs(rho * volumeFlow * c_p * (T_2 - T_1))
    Q_fric = delta_p * volumeFlow
    Q_tot_in = Q_calc + Q_fric    

    transitTime = getTransitTime()
#
    print("transitTime: %f, old transitTime: %f" % (transitTime, old_transitTime))
#
    # Every wait_factor * transitTime a voltage value for the SCR is calculated and set.
    # In case, that the transit time increases to much with the new averaged values, the if statement could be entered two times in a row.
    # To avoid this, the difference has to be smaller than one. Every time the fluid has transitted one time, new values are set and written into the DA file.
    if ((currentTimestep % round(wait_factor * transitTime,0) == 0) or currentTimestep < 2) and transitTime - old_transitTime < 1:
        print("------------------------------------------ In statement -------------------------------------------")
        DAfile_list = open("DA_vals.txt").readlines()
        V_SCR_old = float(DAfile_list[0].split()[1])
        v_scr.append(V_SCR_old)                 # The current V_SCR results
        f_v_scr.append(Q_tot_in)                # in a certain heat inlet.
        
        # If the fluid has transitted 6 times, the control method is activated.
        if 6 * transitTime < currentTimestep:
            
            # In order to get a second grid point, V_SCR is manipulated by a small amount. (only needed one time)
            if (6 + wait_factor) * transitTime < currentTimestep and secondPoint_bool == 1:
                print("-------------------- get second grid point --------------------")

                V_SCR = V_SCR - 0.05             # -0.05 is randomly picked. The only requirement is that V_SCR is changed.
                secondPoint_bool = 0             # deactivate getting a second grid point
                secondPoint_switch = 1           # activate control

            # actual control method   
            if (6 + wait_factor) * transitTime < currentTimestep and secondPoint_switch == 1:
                
                # secant method, modifing to a "root finding problem"
                f_n_minus1 = float(f_v_scr[1]) - Q_set_
                f_n_minus2 = float(f_v_scr[0]) - Q_set_
                x_n_minus1 = float(v_scr[1])
                x_n_minus2 = float(v_scr[0])
                list1 = [float(f_v_scr[1]), float(f_v_scr[0])]
                list2 = [float(v_scr[1]), float(v_scr[0])]
#                
                print(v_scr,f_v_scr)
                print(f_n_minus1, x_n_minus1, f_n_minus2, x_n_minus2)
#                
                # Accounts for "wrong points": If the order that a higher V_SCR gives a higher heat inlet is not provided,
                # no change in V_SCR is made.
                if list1.index(max(list1)) == list2.index(max(list2)):  
                    print("----------------------- correct points ------------------------")
                    
                    # Only if the total heat inlet of the both points is different, the secant step is applied.
                    # Otherwise a division by zero causes errors.
                    if f_n_minus1 != f_n_minus2:    
                        print("------------------------- secant step -------------------------")
                        V_SCR =  x_n_minus1 - f_n_minus1 * (x_n_minus1 - x_n_minus2) / (f_n_minus1 - f_n_minus2)
                        
                # If the V_SCR value remains the same, the control is not active anymore, because the secant method stays at the same point.
                # Therefore, the V_SCR value is manipulated by a small amount.
                if x_n_minus1 == x_n_minus2: 
                    print("----------------------- manipulate V_SCR ----------------------")
                    V_SCR = V_SCR + V_SCR_mod 
                   


        else:
            # V_SCR in the first 6 transit loops
            # Because the control of the SCR based on DArangeSCR is not reliable and provides to high values for V_SCR, 
            # the value for V_SCR is subtracted with 0.3. 
            V_SCR = DArangeSCR * Q_set_ / Q_max + DAlowSCR - 0.3
            
    # control of the pump is made with VFD so far. For the control of the pump with the PI, here is probably the best place.       
        V_pump = 0

    # If the V_SCR should be switched off manually, set V_SCR = 0 here and run the program until the red diode at the SCR is
    # not shining anymore (should occur immediately).
        V_SCR = 0
        
        writeDAfile()

########################################################################
            
# EXPLANATION:
# Because the thermistors have different coefficients for the Steinhart - Hart equation, it is important to know which thermistor, has
# which coefficients. Both thermistors are named thermistor 1 and thermistor 2 in the program.
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
setPoints_clock = time.time()
end_timer = 0
start_timer = 0
currentTimestep = 0
wait_factor = 3             # factor * transit time provides the time, after which a new V_SCR value is set. 
maxTimesteps = sys.maxsize  # maximal number of seconds to run the program (reading in via Setpoints.txt)

f_v_scr = deque(maxlen = 2)
v_scr   = deque(maxlen = 2)

sumV_S = 0
sumV_G = 0
sumV_R1 = 0
sumV_R2 = 0
sumV_flow = 0
sumV_deltaPressure = 0
sumV_Q_elec = 0

counterAverageFunc = 0      # counter for averaging all 7 voltage inputs.
average_switch = 1          # Switch for setting the initiale average values of the voltages.
secondPoint_bool = 1        # Variables, which serve as switches for getting
secondPoint_switch = 0      # the second grid point.

# initialize first set point
setpoints = open("Setpoints.txt").readlines() 
Q_set = float(setpoints[1].split()[1])
volumeFlowrate_set = float(setpoints[1].split()[2])

# initialize logfile
file_obj = open("logfile.txt","w")
file_obj.write("Date/Time\t\t Flow rate(obs)[L/s]\t Heat inlet(tot)[W]\t Heat inlet(fric)[W]\t Heat inlet(elec)[W]\t Temperature 1(obs)[K]\t Temperature 2(obs)[K]\t Pressure drop[Pa]\t V_SCR\t\t V_pump\t\t V_S\t\t V_G\t\t V_R1\t\t V_R2\t\t V_flow\t\t V_deltaPressure\t (no signal)V_Q_elec\n")
file_obj.close()

while True:
    # calculate the current and overall timestep, 
    currentTimestep = round(time.time() - start_clock)
    overallTimestep = round(time.time() - setPoints_clock)
    
    # sleep function in order to ensure that the program is executed one time every second
    time.sleep(math.fabs(1 - (end_timer - start_timer)))
    start_timer = time.time()
    
    # run master script (written in C): reads and sets DA, reads and writes AD
    os.system(command)

    # calculation of DA_values
    calculationOfDA(Q_set, volumeFlowrate_set)
    
    # read setpoints, if required
    if overallTimestep % set_time == 0:
        readSetpoints(overallTimestep, Q_set, volumeFlowrate_set)
            
    # write Logfile
    if currentTimestep % logfile_period == 0:
        print("------------------- write logfile ----------------------")
        writeLogfile()

    end_timer = time.time()
##
    AD_list = open("AD_vals.txt").readlines()
    DA_list = open("DA_vals.txt").readlines()
    print("V_S(ch0): %f, V_G(ch1): %f, V_R1(ch2): %f, V_R2(ch3): %f, V_flow(ch4): %f, V_pres(ch5): %f, noSig-V_Q_elec(ch6): %f, noSig: %f" % (float(AD_list[0].split()[0]),float(AD_list[0].split()[1]),float(AD_list[0].split()[2]),float(AD_list[0].split()[3]),float(AD_list[0].split()[4]),float(AD_list[0].split()[5]),float(AD_list[0].split()[6]),float(AD_list[0].split()[7])))
    print("currentTimestep: %d, DA1: %f" % (currentTimestep, float(DA_list[0].split()[1])))
##    
    if(overallTimestep >= maxTimesteps):
	break

total_end = time.time()
print("End of TRT script, total time: %f" % (total_end-start_clock))
