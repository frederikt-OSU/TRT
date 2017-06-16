#include "bcm2835.h"
#include "ads1256_input.h"
#include "dac8532_output.h"
#include "master.h"
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*	Explanation:
 * get_input(float [], int) : is the main function from ads1256_test --> floating values with volts as unit
 * set_output(float[],float): is the main function from dac8532_test --> floating values with volts as unit
 * input_voltage			: 0 -->  V_T_ref1, 3 --> V_T2
 * 							  1 -->  V_T1	 , 4 --> V_flow
 * 							  2 -->  V_Tref_2, 5 --> V_pressure 
 * V_ref_set 				: reference voltage for setting of voltage output (5.0 or 3.3 V)
*/
/*
void read_setvalues(float *Q, float *volumeflowrate, int length_firstline_Setfile, int *maxTime, int currentTime){
	
	printf("======================== read setfile ===========================\n");
	
	int num_lines = -1;	// first line is not counted, because of text
	int counter = 0;
	char line[200];
	char buf1[15],buf2[15],buf3[15];
	char * pbuf1 = buf1;
	char * pbuf2 = buf2;
	char * pbuf3 = buf3;
	FILE * setfile_pointer;
	
	setfile_pointer = fopen("Setfile","r+");
	fseek(setfile_pointer,length_firstline_Setfile,SEEK_SET);
	
	// count number input lines
	while(fgets(line,sizeof(line),setfile_pointer) != 0){		
			num_lines++;
	}

	float set_information[num_lines][3];	// 3 because of settime, Q and volumeflowrate	
	fseek(setfile_pointer,length_firstline_Setfile,SEEK_SET);
	
	while(counter < num_lines){
		fscanf(setfile_pointer, "%s\t\t%s\t\t%s",buf1, buf2, buf3);
		
		set_information[counter][0] = atof(pbuf1);	
		set_information[counter][1] = atof(pbuf2);
		set_information[counter][2] = atof(pbuf3);
		
		// read the last time step
		if(counter == (num_lines-1)) 
			*maxTime = (int) set_information[counter][0];
		
		// set heatinlet and flowrate
		if(currentTime < (int) atof(pbuf1)){
			*Q 			    = set_information[counter-1][1];
			*volumeflowrate = set_information[counter-1][2];			
			break;
		}			
		counter++;
	}
		
	printf("==================== read setfile SUCCESSFUL ====================\n");	
	fclose(setfile_pointer);

} */

//main function for controlling the pi input and output
int main(){
  
	int num_input_ch = 6;				// up to 8 possible	
	float DA_values[2];	
	float input_voltage[num_input_ch];
	float V_ref_set;	
	FILE * setfile_pointer;
	FILE * ADfile_pointer;
	char buf1[15],buf2[15],buf3[15];
	char * pbuf1 = buf1;
	char * pbuf2 = buf2;
	char * pbuf3 = buf3;		    
	
	// Read and set DA_values
	setfile_pointer = fopen("DA_vals.txt","r+");
	fscanf(setfile_pointer, "%s\t%s\t%s",buf1, buf2, buf3);
	fclose(setfile_pointer);
			
	DA_values[0] = atof(pbuf1);	
	DA_values[1] = atof(pbuf2);
	V_ref_set	 = atof(pbuf3);	

	set_output_voltage(DA_values, V_ref_set);	  
	
	// Read AD values and write file
	get_input(input_voltage, num_input_ch);	
	
	ADfile_pointer = fopen("AD_vals.txt","w+");
	fprintf(ADfile_pointer,"%f\t%f\t%f\t%f\t%f\t%f", input_voltage[0], input_voltage[1], input_voltage[2], input_voltage[3], input_voltage[4], input_voltage[5]);
	fclose(ADfile_pointer);
			
	//printf("AD values: %f\t%f\t%f\t%f\t%f\t%f\n", input_voltage[0], input_voltage[1], input_voltage[2], input_voltage[3], input_voltage[4], input_voltage[5]);
	//printf("DA0 values: %f, DA1 values: %f\n", DA_values[0], DA_values[1]);

	return 0;  
}  
 
/*	Do loop and time handeling which should be implemented in python
 * 	///////////////////////////// initialization ///////////////////////
	int num_input_ch = 6;				// up to 8 possible
	int i,p, currentTime;
	int tmpTime = -1;
	int maxTime = INT_MAX;
	int set_time = 3;					// timeperiod to check the setfile in seconds 
	int length_firstline_Setfile = 27;	// length of the very first line in the setfile
	
	float output_voltage[2];	
	float input_voltage[num_input_ch];
	float V_ref_set = 5.0;				// 5.0 or 3.3	
	float Q_set, volumeflowrate_set;
	FILE * fPointer;
	time_t date_and_time;
	clock_t start_clock,end_clock;
	
	// start of taking time
	start_clock = clock();
	fPointer = fopen("InAndOutputData.txt","w");
	fprintf(fPointer, "Input(ch0) \t Input(ch1) \t Input(ch2) \t Input(ch3) \t Input(ch4) \t Input(ch5) \t Output(ch0) \t Output(ch1) \t Date/Time \n");
	fclose(fPointer);		
	////////////////////////////////////////////////////////////////////
	* 
	* 
	do{		
		end_clock = clock();
		currentTime = (( end_clock - start_clock)/ CLOCKS_PER_SEC) /1;	// Todo: durch 3600 teilen um auf stunden zu kommen oder ueber set_time regeln??
		
		// check if one second is passed
		if (currentTime == tmpTime +1){		

			//check if it is necessary to check and set the heatinlet and volumeflowrate
			if((currentTime % set_time) == 0){	
				read_setvalues(&Q_set, &volumeflowrate_set, length_firstline_Setfile, &maxTime, currentTime);
				printf("Setvalues: currentTime = %d, Q = %f, flow = %f\n", currentTime, Q_set, volumeflowrate_set);
			}

			get_input(input_voltage, num_input_ch);			
			calculation_of_output(input_voltage, output_voltage, Q_set, volumeflowrate_set);
			
			set_output_voltage(output_voltage, V_ref_set);
			//getchar();
				
			// write data to output file
			date_and_time = time(0);			
			fPointer = fopen("InAndOutputData.txt","a");
			fprintf(fPointer, "%f, \t %f, \t %f, \t %f, \t %f, \t %f, \t %f, \t %f, \t %s", input_voltage[0], input_voltage[1], input_voltage[2], input_voltage[3], input_voltage[4], input_voltage[5], output_voltage[0], output_voltage[1], ctime(&date_and_time));
			fclose(fPointer);			
			tmpTime = currentTime;
		}
	} while(currentTime < maxTime);
*/



