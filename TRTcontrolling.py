import os
import time
from threading import Timer

 
#########################    functions   #################################





######################### execution part ################################
command = "sudo ./master"

start_clock = time.time()

tmpTime = -1
maxTimesteps = 3

while True:
    end_clock = time.time()
    currentTimestep = round(end_clock - start_clock,0)
#    print(currentTimestep,tmpTime+1)
	
    if currentTimestep == (tmpTime + 1):
        
	tmpTime = currentTimestep
#	print("currentTime: %d, tmpTime: %d "% (currentTimestep, tmpTime))
        # run master script: read and set DA, read  and write AD
	master_clock1 = time.time()
        os.system(command)
	master_clock2 = time.time()
	print("time masterscript: %f\n" % (master_clock2-master_clock1))

	# read setpoints, if required


	# calculation


	# write Logfile and DA values



    if(currentTimestep >= maxTimesteps):
        break



total_end = time.time()
print("End of TRT script, total time: %f" % (total_end-start_clock))
