/*
 * DAC8532_test.c:
 *	Very simple program to test the serial port. Expects
 *	the port to be looped back to itself
 *
 */
 
/*
             define from bcm2835.h                       define from Board 
DVK511
                 3.3V | | 5V               ->                 3.3V | | 5V
    RPI_V2_GPIO_P1_03 | | 5V               ->                  SDA | | 5V 
    RPI_V2_GPIO_P1_05 | | GND              ->                  SCL | | GND
       RPI_GPIO_P1_07 | | RPI_GPIO_P1_08   ->                  IO7 | | TX
                  GND | | RPI_GPIO_P1_10   ->                  GND | | RX
       RPI_GPIO_P1_11 | | RPI_GPIO_P1_12   ->                  IO0 | | IO1
    RPI_V2_GPIO_P1_13 | | GND              ->                  IO2 | | GND
       RPI_GPIO_P1_15 | | RPI_GPIO_P1_16   ->                  IO3 | | IO4
                  VCC | | RPI_GPIO_P1_18   ->                  VCC | | IO5
       RPI_GPIO_P1_19 | | GND              ->                 MOSI | | GND
       RPI_GPIO_P1_21 | | RPI_GPIO_P1_22   ->                 MISO | | IO6
       RPI_GPIO_P1_23 | | RPI_GPIO_P1_24   ->                  SCK | | CE0
                  GND | | RPI_GPIO_P1_26   ->                  GND | | CE1

::if your raspberry Pi is version 1 or rev 1 or rev A
RPI_V2_GPIO_P1_03->RPI_GPIO_P1_03
RPI_V2_GPIO_P1_05->RPI_GPIO_P1_05
RPI_V2_GPIO_P1_13->RPI_GPIO_P1_13
::
*/

#include <bcm2835.h>  
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <errno.h>

//CS      -----   SPICS_  
//DIN     -----   MOSI
//DOUT  -----   MISO
//SCLK   -----   SCLK
//DRDY  -----   ctl_IO     data  starting
//RST     -----   ctl_IO     reset



#define	SPICS_	RPI_GPIO_P1_16	//P4 	//modified SPICS to SPICS_

#define CS_1_() bcm2835_gpio_write(SPICS_,HIGH) // modified CS_1() to CS_1()
#define CS_0_() bcm2835_gpio_write(SPICS_,LOW) // modified CS_0() to CS_0()


/* Unsigned integer types  */
#define uint8_t unsigned char
#define uint16_t unsigned short    

#define channel_A   0x30
#define channel_B   0x34
	  	
//void  bsp_DelayUS(uint64_t micros);
void Write_DAC8532(uint8_t channel, uint16_t Data);
uint16_t Voltage_Convert_(float Vref, float voltage);


void  bsp_DelayUS_(uint64_t micros)
{
		bcm2835_delayMicroseconds (micros);
}


void Write_DAC8532(uint8_t channel, uint16_t Data)
{
	uint8_t i;

	 CS_1_() ;
	 CS_0_() ;
      bcm2835_spi_transfer(channel);
      bcm2835_spi_transfer((Data>>8));
      bcm2835_spi_transfer((Data&0xff));  
      CS_1_() ;
}

uint16_t Voltage_Convert(float Vref, float voltage)
{
	uint16_t _D_;
	_D_ = (uint16_t)(65536 * voltage / Vref);
    
	return _D_;
}


int  set_output_voltage(float output_values[], float V_ref_set)
{

   uint16_t   i,tmp;
   
	if (!bcm2835_init())
        return 1;

    bcm2835_spi_begin();
    bcm2835_spi_setBitOrder(BCM2835_SPI_BIT_ORDER_LSBFIRST );      // The default
    bcm2835_spi_setDataMode(BCM2835_SPI_MODE1);                   // The default;
    bcm2835_spi_setClockDivider(BCM2835_SPI_CLOCK_DIVIDER_1024); // The default
    bcm2835_gpio_fsel(SPICS_, BCM2835_GPIO_FSEL_OUTP);//
    bcm2835_gpio_write(SPICS_, HIGH); 	

	Write_DAC8532(0x30, Voltage_Convert(V_ref_set,output_values[0]));    //Write channel A buffer (0x30) --> DAC0
    Write_DAC8532(0x34, Voltage_Convert(V_ref_set,output_values[1]));    //Write channel B buffer (0x34) --> DAC1

    bcm2835_spi_end();
    bcm2835_close();
	
    return 0;
}
