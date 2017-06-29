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
 * V_ref_set 				: reference voltage for setting of voltage output (5.0 or 3.3 V)
*/

//main function for getting and setting input and output
int main(){
  
	int num_input_ch = 8;				// up to 8 possible	
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
	fprintf(ADfile_pointer,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", input_voltage[0], input_voltage[1], input_voltage[2], input_voltage[3], input_voltage[4], input_voltage[5], input_voltage[6], input_voltage[7]);
	fclose(ADfile_pointer);
			
	//printf("AD values: %f\t%f\t%f\t%f\t%f\t%f\n", input_voltage[0], input_voltage[1], input_voltage[2], input_voltage[3], input_voltage[4], input_voltage[5]);
	//printf("DA0 values: %f, DA1 values: %f\n", DA_values[0], DA_values[1]);

	return 0;  
}  
