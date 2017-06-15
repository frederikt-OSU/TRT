import os
import time
import datetime
import threading
#import _thread
from threading import Timer

 
######################## functions ##########################



def counter():
    currentTimestep += 1
    #myTimer.start()





############################# main function ###############################



def main():
    #myTimer.start()
    tmpTime = 0
    currentTimestep = 1
    while True:
        if currentTimestep == tmpTime + 1:

	
            print("currentTime: %d, tmpTime: %d "% (currentTimestep, tmpTime))

		
            
            maxTimesteps = 6
            os.system(command)
            tmpTime += 1

            if(currentTimestep >= maxTimesteps):
		break





######################### execution part ################################
command = "sudo ./master"

#myTimer = Timer(1.0,counter)


main()



    
