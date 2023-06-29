# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:18:04 2023

@author: Paula Sánchez Almagro
"""

#%% Importing libraries, files and setting general parameters
import numpy as np
import matplotlib.pyplot as plt
import write_memfile as w8 
import math

ref_clk=18e9;
if(ref_clk>20e9):
    raise ValueError('Reference clock exceed limit of 20 GHz');
mem_end=32768;
row_size=32;

fs=ref_clk * 2;
freq_1row =fs/32;
T_1row=1/freq_1row;
T_1symb=T_1row/32;
delta_t=1/fs;

#%% Segment signals

class segment:
    def  __init__ (self, start=0, length=0, repeat_msb=0, repeat_lsb=0, wait_msb=0, wait_mid=0, wait_lsb=0, next_segment=0):
        self.start = start;
        self.length =length;
        self.repeat_msb = repeat_msb;
        self.repeat_lsb =repeat_lsb;
        self.wait_msb =wait_msb;
        self.wait_mid =wait_mid;
        self.wait_lsb =wait_lsb;
        self.next_segment = next_segment;
        
        self.symb= 0;
        self.yi = np.zeros(mem_end);   
        
    def lin_test(self):
        yi=np.zeros(mem_end);

        for i in range (0,256):
            j=i*128;
            for k in range (j, j+128):
                yi[k]=i;    
            
        filename = "mem_lin";
        header = "# Linearity test ramp 0 to 255";
        
        self.yi=yi;
        self.length=len(self.yi)/32;


    def ramp_up (self, full_memory, lowest_val_ramp, highest_val_ramp, Time_ramp):
        self.full_memory = full_memory;
        self.lowest_val_ramp = lowest_val_ramp;
        self.highest_val_ramp =highest_val_ramp;
        self.Time_ramp = Time_ramp;
        
        yi=np.zeros(mem_end);

        if (self.highest_val_ramp < self.lowest_val_ramp):
            raise ValueError('Parameters highest and lowest value for the ramp not correct');
            
        N_symb = int (self.Time_ramp/T_1symb);
            
        for i in range (0,mem_end-N_symb):
            yi[i]=self.lowest_val_ramp;
        for i in range (mem_end-N_symb, mem_end):
            yi[i]=(i*(self.highest_val_ramp-self.lowest_val_ramp)+self.lowest_val_ramp*(mem_end-1)-self.highest_val_ramp*(mem_end-N_symb))*1/(N_symb-1);
        
        self.symb=N_symb; 
        self.length=math.ceil(N_symb/(row_size))-1;
            
        if (full_memory):
            filename = "mem_ramp_up_fm";
            header = "#Ramp up test full mem";
            self.yi = yi;
            self.length=len(self.yi)/32;
        else:   
            filename = "mem_ramp_up";
            header = "#Ramp up test";
            self.yi = yi[mem_end-((self.length+1)*32):mem_end];
            self.length=len(self.yi)/32;

            
    def ramp_down (self, full_memory, lowest_val_ramp, highest_val_ramp, Time_ramp):
        self.lowest_val_ramp = lowest_val_ramp;
        self.highest_val_ramp =highest_val_ramp;
        self.Time_ramp = Time_ramp;
        
        yi=np.zeros(mem_end);

        if (self.highest_val_ramp < self.lowest_val_ramp):
            raise ValueError('Parameters highest and lowest value for the ramp not correct');
            
        N_symb = int (self.Time_ramp/T_1symb);
        
        for i in range (0,N_symb-1):
            yi[i]=(i*(self.lowest_val_ramp-self.highest_val_ramp)+self.highest_val_ramp*(N_symb-1))*1/(N_symb-1);
        for i in range (N_symb-1, mem_end):
            yi[i]=self.lowest_val_ramp;       
        
        self.symb=N_symb; 
        self.length=math.ceil(N_symb/(row_size))-1;
        
        if (full_memory):
            filename = "mem_ramp_down_fm";
            header = "#Ramp down test full mem";
            self.yi = yi;
            self.length=len(self.yi)/32;
        else:   
            filename = "mem_ramp_down";
            header = "#Ramp down test";
            self.yi = yi[0:((self.length+1)*32)];
            self.length=len(self.yi)/32;
            
            
    def sine (self, full_memory, frequency, offset, amp):
        self.frequency =frequency;
        self.offset = offset;
        self.amp=amp;
        
        numb_freqs=len(self.frequency);
        
        flag_several_rows = True;
        
        if(numb_freqs==0):
            raise ValueError('Frequency not defined');
        elif (numb_freqs==1):
            N = freq_1row/self.frequency[0];
            if(N>=1):
                N=np.floor(N);
                freq_1row_new = N * self.frequency[0];
                flag_several_rows = True;                
            else:
                N2 = np.ceil(self.frequency[0]/freq_1row);
                freq_1row_new = self.frequency[0]/N2;
                flag_several_rows = False;

            ref_clk = freq_1row_new * 16;
            if(ref_clk>10e9):
                raise ValueError('Reference clock exceed limit of 10e9 Hz');
            print('New ref_clk: ' + str(ref_clk));
            
            fs=ref_clk * 2;
            delta_t=1/fs; 
            time=(np.arange(0, mem_end, 1))*delta_t;
            yi=np.zeros(mem_end);
            yi= self.amp * np.sin(2 * np.pi * self.frequency[0] * time);            
            for k in range (0, len(time)):
                yi[k]=yi[k] + self.offset;
            if (full_memory):
                self.yi = yi;
                self.length=len(self.yi)/32;
            else:
                if(N>1):
                    self.yi = yi[0:int(32*N)];
                    self.length=len(self.yi)/32;
                else:
                    self.yi = yi[0:32];
                    self.length=len(self.yi)/32;
        else:
            time=(np.arange(0, mem_end, 1))*delta_t;
            yi=np.zeros(mem_end);
            for f in self.frequency:
               yi += self.amp/numb_freqs * np.sin(2 * np.pi * int(f) * time);
            for k in range (0, len(time)):
                yi[k]=yi[k] + self.offset;
            self.yi = yi;
            self.length=len(self.yi)/32;
                
        if (full_memory):
            filename = "mem_sine_fm";
            header = "#Sine test full mem";
        else:  
            filename = "mem_sine";
            header = "#Sine test";
        
                        

    def rcos (self, frequency, offset, amp, T, beta, span=mem_end):
        self.frequency =frequency;
        self.offset = offset;
        self.amp = amp;
        self.T = T;
        self.beta = beta;
        self.span = span;
        
        numb_freqs=len(self.frequency);
        
        T_delta_filter = 1/float(fs);
        sample_num_filter = np.arange(self.span);
        h_rc = np.zeros(self.span, dtype=float);
            
        for x in sample_num_filter:
            t_filter = (x-self.span/2)*T_delta_filter;
            if t_filter == 0.0:
                h_rc[x] = 1.0;
            elif self.beta != 0 and t_filter == self.T/(2*self.beta):
                h_rc[x] = (np.pi/4)*(np.sin(np.pi*t_filter/self.T)/(np.pi*t_filter/self.T)); 
            elif self.beta != 0 and t_filter == -self.T/(2*self.beta):
                h_rc[x] = (np.pi/4)*(np.sin(np.pi*t_filter/self.T)/(np.pi*t_filter/self.T));
            else:
                h_rc[x] = (np.sin(np.pi*t_filter/self.T)/(np.pi*t_filter/self.T))* \
                        (np.cos(np.pi*self.beta*t_filter/self.T)/(1-(((2*self.beta*t_filter)/self.T)*((2*self.beta*t_filter)/self.T))));
        
        if(self.beta<0)or(self.beta>1):
            raise ValueError('Beta parameter must be in [0,1]');
        
        time=(np.arange(0, mem_end, 1))*delta_t;
        if(numb_freqs==0):
            raise ValueError('Frequency not defined');
        elif (numb_freqs==1):
            yi = h_rc *self.amp * np.sin(2 * np.pi * self.frequency[0] * time);      
            
        else:
            yi=np.zeros(mem_end);
            for i in range (0, numb_freqs):
                yi = yi + h_rc *self.amp * np.sin(2 * np.pi * self.frequency[i] * time);
            
        yi=yi+self.offset;
        self.yi=yi;
        self.length=len(self.yi)/32;
        
        filename = "mem_rcos_fm";
        header = "#Rcos test full mem";
        


    def gauss (self, full_memory, frequency, offset, amp, sigma, length=mem_end):
        self.frequency =frequency;
        self.offset = offset;
        self.amp = amp;
        self.sigma = sigma;
        self.length = length;        
        
        numb_freqs=len(self.frequency);
        
        
        time=(np.arange(0, mem_end, 1))*delta_t;
        x = np.linspace(-self.length/2, self.length/2, self.length);
        window = np.exp(-(x**2) / (2 * self.sigma**2));
        
        window_length = len([a for a in window if a > 0]);
        if(window_length%32 !=0):
            window_length = np.ceil(window_length/32) *32;
            
        
        if(numb_freqs==0):
            raise ValueError('Frequency not defined');
        elif (numb_freqs==1):
            sine=self.amp * np.sin(2 * np.pi * self.frequency[0] * ((time+(1/(2*self.frequency[0])))));           
        else:
            sine=np.zeros(mem_end);
            for f in self.frequency:
               sine += self.amp * np.sin(2 * np.pi * int(f) * time);
            
        yi=window*sine;   
        yi=yi+self.offset;
        
        if (full_memory):
            filename = "mem_gauss_fm";
            header = "#Gauss test full mem";
            self.yi = yi;
            self.length=len(self.yi)/32;
        else:   
            filename = "mem_gauss";
            header = "#Gauss test";
            self.yi = yi[int((mem_end/2)-(window_length/2)):int((mem_end/2)+(window_length/2))]; 
            self.length=len(self.yi)/32;
        
        

    def clk (self, frequency, low_value, high_value, periods=0): 
        self.frequency = frequency;  
        self.low_value = low_value;
        self.high_value = high_value;
        self.periods = periods;
        
        yi=np.zeros(mem_end);

        if (self.high_value < self.low_value):
            raise ValueError('Parameters high and lowe value for the clk not correct');
            
        period = 1/self.frequency;
        half_wave_off = np.floor(period/delta_t/2);
        half_wave_on = np.ceil(period/delta_t/2);
        count=0;
        flag=True;
        
        for i in range (0, mem_end):
            count = count + 1;
            if (flag):
                if (count > half_wave_on):
                    count = 0;
                    flag = not (flag);
                yi[i]=self.high_value;
            else:
                if (count > half_wave_off):
                    count = 0;
                    flag = not (flag);
                yi[i]=self.low_value;    
                    
        filename = "mem_clk_fm";
        header = "#Clk test full mem";
        self.yi = yi;
        self.length=len(self.yi)/32;

    
     
#%% regpgenQC

class regpgenQC:
    # def __init__ (self, s0=None, s1=None, s2=None, s3=None, s4=None, s5=None, s6=None, s7=None):
    def __init__ (self, s0, s1, s2, s3, s4, s5, s6, s7):
        self.s0 = s0;
        self.s1 = s1;   
        self.s2 = s2;
        self.s3 = s3;
        self.s4 = s4;
        self.s5 = s5;
        self.s6 = s6;
        self.s7 = s7;
        self.segments=[self.s0, self.s1, self.s2, self.s3, self.s4, self.s5, self.s6, self.s7];
    
    def generate_mem_file(self):
        y=np.zeros(0);
        for i in range (0, 8):
            if(self.segments[i]!=None):
                y=np.concatenate((y, self.segments[i].yi));
        
        plt.title('Signal shape')
        plt.xlabel('Data memory')
        plt.ylabel('Number of output levels')
        plt.plot(y);
        
        filename = "mem.sequence";
        header = "#mem sequence";
        w8.write_memfile(filename, header, y, mem_end, row_size);      

        
    def generate_rg_file(self):
        f=open("reg_pgenQC_pyfile", "w");
        f.write("# GUI_ORDER = 20\n## FOREGROUND_COLOR = #111111\n## BACKGROUND_COLOR = light blue\n\n#\n# comments and blank lines allowed\n# only first two columns (addr/data) are parsed/used\n# supported data formats (dec,hex,oct,bin), e.g.:  25  0x19  031  0b11001\n# underscores can be used to group digits, e.g.: 10_000  0x3a_e7  0b10_111_01_1000\n#\n# addr	data	name\n\n");
        
        if(self.s0 != None):
            f.write("0x3a0	" + str(int(self.s0.start)) +"	3/s0_start		# 2:11	10 bit SRAM address\n");
            f.write("0x3a1	" + str(int(self.s0.length-1)) +"	3/s0_length		# 2:11	10 bit, length of sequence (in SRAM words)\n");
            f.write("0x3a2	" + str(int(self.s0.repeat_msb)) +"	3/s0_repeat_msb		# 8:11	16 bit, number of repeats (of entire sequence)\n");
            f.write("0x3a3	" + str(int(self.s0.repeat_lsb)) + "	3/s0_repeat_lsb		# 0:11\n");
            f.write("0x3a4	" + str(int(self.s0.wait_msb)) + "	3/s0_wait_msb		# 4:11	32 bit, wait after sequence (in clock cycles)\n");
            f.write("0x3a5	" + str(int(self.s0.wait_mid)) + "	3/s0_wait_mid		# 0:11\n");
            f.write("0x3a6	" + str(int(self.s0.wait_lsb)) + "	3/s0_wait_lsb		# 0:11\n");
            f.write("0x3a7	" + str(int(self.s0.next_segment)) + "	3/s0_next		# 8:11	4 bit, next seq. number. Special value 0xf can be used to end pattern and go back to idle\n\n");
        else:
            f.write("0x3a0	0	3/s0_start		# 2:11	10 bit SRAM address\n");     
            f.write("0x3a1	0	3/s0_length		# 2:11	10 bit, length of sequence (in SRAM words)\n");
            f.write("0x3a2	0	3/s0_repeat_msb		# 8:11	16 bit, number of repeats (of entire sequence)\n");
            f.write("0x3a3	0	3/s0_repeat_lsb		# 0:11\n");
            f.write("0x3a4	0	3/s0_wait_msb		# 4:11	32 bit, wait after sequence (in clock cycles)\n");
            f.write("0x3a5	0	3/s0_wait_mid		# 0:11\n");
            f.write("0x3a6	0	3/s0_wait_lsb		# 0:11\n");
            f.write("0x3a7	0	3/s0_next		# 8:11	4 bit, next seq. number. Special value 0xf can be used to end pattern and go back to idle\n\n");
        
        
        if(self.s1 != None):
            f.write("0x3a8	" + str(int(self.s1.start+self.s0.length)) +"	3/s1_start\n");
            f.write("0x3a9	" + str(int(self.s1.length-1)) +"	3/s1_length\n");
            f.write("0x3aa	" + str(int(self.s1.repeat_msb)) +"	3/s1_repeat_msb\n");
            f.write("0x3ab	" + str(int(self.s1.repeat_lsb)) + "	3/s1_repeat_lsb\n");
            f.write("0x3ac	" + str(int(self.s1.wait_msb)) + "	3/s1_wait_msb\n");
            f.write("0x3ad	" + str(int(self.s1.wait_mid)) + "	3/s1_wait_mid\n");
            f.write("0x3ae	" + str(int(self.s1.wait_lsb)) + "	3/s1_wait_lsb\n");
            f.write("0x3af	" + str(int(self.s1.next_segment)) + "	3/s1_next\n\n");
        else:
            f.write("0x3a8	0	3/s1_start\n");
            f.write("0x3a9	0	3/s1_length\n");
            f.write("0x3aa	0	3/s1_repeat_msb\n");
            f.write("0x3ab	0	3/s1_repeat_lsb\n");
            f.write("0x3ac	0	3/s1_wait_msb\n");
            f.write("0x3ad	0	3/s1_wait_mid\n");
            f.write("0x3ae	0	3/s1_wait_lsb\n");
            f.write("0x3af	0	3/s1_next\n\n");
        
        if(self.s2 != None):
            f.write("0x3b0	" + str(int(self.s2.start+self.s0.length+self.s1.length)) +"	3/s2_start\n");
            f.write("0x3b1	" + str(int(self.s2.length-1)) +"	3/s2_length\n");
            f.write("0x3b2	" + str(int(self.s2.repeat_msb)) +"	3/s2_repeat_msb\n");
            f.write("0x3b3	" + str(int(self.s2.repeat_lsb)) + "	3/s2_repeat_lsb\n");
            f.write("0x3b4	" + str(int(self.s2.wait_msb)) + "	3/s2_wait_msb\n");
            f.write("0x3b5	" + str(int(self.s2.wait_mid)) + "	3/s2_wait_mid\n");
            f.write("0x3b6	" + str(int(self.s2.wait_lsb)) + "	3/s2_wait_lsb\n");
            f.write("0x3b7	" + str(int(self.s2.next_segment)) + "	3/s2_next\n\n");
        else:
            f.write("0x3b0	0	3/s2_start\n");
            f.write("0x3b1	0	3/s2_length\n");
            f.write("0x3b2	0	3/s2_repeat_msb\n");
            f.write("0x3b3	0	3/s2_repeat_lsb\n");
            f.write("0x3b4	0	3/s2_wait_msb\n");
            f.write("0x3b5	0	3/s2_wait_mid\n");
            f.write("0x3b6	0	3/s2_wait_lsb\n");
            f.write("0x3b7	0	3/s2_next\n\n");
        
        if(self.s3 != None):
            f.write("0x3b8	" + str(int(self.s3.start+self.s0.length+self.s1.length+self.s2.length)) +"	3/s3_start\n");
            f.write("0x3b9	" + str(int(self.s3.length-1)) +"	3/s3_length\n");
            f.write("0x3ba	" + str(int(self.s3.repeat_msb)) +"	3/s3_repeat_msb\n");
            f.write("0x3bb	" + str(int(self.s3.repeat_lsb)) + "	3/s3_repeat_lsb\n");
            f.write("0x3bc	" + str(int(self.s3.wait_msb)) + "	3/s3_wait_msb\n");
            f.write("0x3bd	" + str(int(self.s3.wait_mid)) + "	3/s3_wait_mid\n");
            f.write("0x3be	" + str(int(self.s3.wait_lsb)) + "	3/s3_wait_lsb\n");
            f.write("0x3bf	" + str(int(self.s3.next_segment)) + "	3/s3_next\n\n");
        else:
            f.write("0x3b8	0	3/s3_start\n");
            f.write("0x3b9	0	3/s3_length\n");
            f.write("0x3ba	0	3/s3_repeat_msb\n");
            f.write("0x3bb	0	3/s3_repeat_lsb\n");
            f.write("0x3bc	0	3/s3_wait_msb\n");
            f.write("0x3bd	0	3/s3_wait_mid\n");
            f.write("0x3be	0	3/s3_wait_lsb\n");
            f.write("0x3bf	0	3/s3_next\n\n");
        
        if(self.s4 != None):
            f.write("0x3c0	" + str(int(self.s4.start+self.s0.length+self.s1.length+self.s2.length+self.s3.length)) +"	3/s4_start\n");
            f.write("0x3c1	" + str(int(self.s4.length-1)) +"	3/s4_length\n");
            f.write("0x3c2	" + str(int(self.s4.repeat_msb)) +"	3/s4_repeat_msb\n");
            f.write("0x3c3	" + str(int(self.s4.repeat_lsb)) + "	3/s4_repeat_lsb\n");
            f.write("0x3c4	" + str(int(self.s4.wait_msb)) + "	3/s4_wait_msb\n");
            f.write("0x3c5	" + str(int(self.s4.wait_mid)) + "	3/s4_wait_mid\n");
            f.write("0x3c6	" + str(int(self.s4.wait_lsb)) + "	3/s4_wait_lsb\n");
            f.write("0x3c7	" + str(int(self.s4.next_segment)) + "	3/s4_next\n\n");
        else:
            f.write("0x3c0	0	3/s4_start\n");
            f.write("0x3c1	0	3/s4_length\n");
            f.write("0x3c2	0	3/s4_repeat_msb\n");
            f.write("0x3c3	0	3/s4_repeat_lsb\n");
            f.write("0x3c4	0	3/s4_wait_msb\n");
            f.write("0x3c5	0	3/s4_wait_mid\n");
            f.write("0x3c6	0	3/s4_wait_lsb\n");
            f.write("0x3c7	0	3/s4_next\n\n");
        
        if(self.s5 != None):
            f.write("0x3c8	" + str(int(self.s5.start+self.s0.length+self.s1.length+self.s2.length+self.s3.length+self.s4.length)) +"	3/s5_start\n");
            f.write("0x3c9	" + str(int(self.s5.length-1)) +"	3/s5_length\n");
            f.write("0x3ca	" + str(int(self.s5.repeat_msb)) +"	3/s5_repeat_msb\n");
            f.write("0x3cb	" + str(int(self.s5.repeat_lsb)) + "	3/s5_repeat_lsb\n");
            f.write("0x3cc	" + str(int(self.s5.wait_msb)) + "	3/s5_wait_msb\n");
            f.write("0x3cd	" + str(int(self.s5.wait_mid)) + "	3/s5_wait_mid\n");
            f.write("0x3ce	" + str(int(self.s5.wait_lsb)) + "	3/s5_wait_lsb\n");
            f.write("0x3cf	" + str(int(self.s5.next_segment)) + "	3/s5_next\n\n");
        else:
            f.write("0x3c8	0	3/s5_start\n");
            f.write("0x3c9	0	3/s5_length\n");
            f.write("0x3ca	0	3/s5_repeat_msb\n");
            f.write("0x3cb	0	3/s5_repeat_lsb\n");
            f.write("0x3cc	0	3/s5_wait_msb\n");
            f.write("0x3cd	0	3/s5_wait_mid\n");
            f.write("0x3ce	0	3/s5_wait_lsb\n");
            f.write("0x3cf	0	3/s5_next\n\n");
        
        if(self.s6 != None):
            f.write("0x3d0	" + str(int(self.s6.start+self.s0.length+self.s1.length+self.s2.length+self.s3.length+self.s4.length+self.s5.length)) +"	3/s6_start\n");
            f.write("0x3d1	" + str(int(self.s6.length-1)) +"	3/s6_length\n");
            f.write("0x3d2	" + str(int(self.s6.repeat_msb)) +"	3/s6_repeat_msb\n");
            f.write("0x3d3	" + str(int(self.s6.repeat_lsb)) + "	3/s6_repeat_lsb\n");
            f.write("0x3d4	" + str(int(self.s6.wait_msb)) + "	3/s6_wait_msb\n");
            f.write("0x3d5	" + str(int(self.s6.wait_mid)) + "	3/s6_wait_mid\n");
            f.write("0x3d6	" + str(int(self.s6.wait_lsb)) + "	3/s6_wait_lsb\n");
            f.write("0x3d7	" + str(int(self.s6.next_segment)) + "	3/s6_next\n\n");
        else:
            f.write("0x3d0	0	3/s6_start\n");
            f.write("0x3d1	0	3/s6_length\n");
            f.write("0x3d2	0	3/s6_repeat_msb\n");
            f.write("0x3d3	0	3/s6_repeat_lsb\n");
            f.write("0x3d4	0	3/s6_wait_msb\n");
            f.write("0x3d5	0	3/s6_wait_mid\n");
            f.write("0x3d6	0	3/s6_wait_lsb\n");
            f.write("0x3d7	0	3/s6_next\n\n");
        
        if(self.s7 != None):
            f.write("0x3d8	" + str(int(self.s7.start+self.s0.length+self.s1.length+self.s2.length+self.s3.length+self.s4.length+self.s5.length+self.s6.length)) +"	3/s7_start\n");
            f.write("0x3d9	" + str(int(self.s7.length-1)) +"	3/s7_length\n");
            f.write("0x3da	" + str(int(self.s7.repeat_msb)) +"	3/s7_repeat_msb\n");
            f.write("0x3db	" + str(int(self.s7.repeat_lsb)) + "	3/s7_repeat_lsb\n");
            f.write("0x3dc	" + str(int(self.s7.wait_msb)) + "	3/s7_wait_msb\n");
            f.write("0x3dd	" + str(int(self.s7.wait_mid)) + "	3/s7_wait_mid\n");
            f.write("0x3de	" + str(int(self.s7.wait_lsb)) + "	3/s7_wait_lsb\n");
            f.write("0x3df	" + str(int(self.s7.next_segment)) + "	3/s7_next\n\n");
        else:
            f.write("0x3d8	0	3/s7_start\n");
            f.write("0x3d9	0	3/s7_length\n");
            f.write("0x3da	0	3/s7_repeat_msb\n");
            f.write("0x3db	0	3/s7_repeat_lsb\n");
            f.write("0x3dc	0	3/s7_wait_msb\n");
            f.write("0x3dd	0	3/s7_wait_mid\n");
            f.write("0x3de	0	3/s7_wait_lsb\n");
            f.write("0x3df	0	3/s7_next\n\n");
        
        f.close();
        
        
#%%  ____________________________________________ EXAMPLE  ____________________________________________

#%%
#(full_memory, lowest_val_ramp, highest_val_ramp, Time_ramp):
premanip=segment()
premanip.ramp_up(False, 0, 127, 40e-9)

#%%
#(full_memory, lowest_val_ramp, highest_val_ramp, Time_ramp):
postmanip=segment()
postmanip.ramp_down(False, 0, 127, 40e-9)#%%

#%%(full_memory, lowest_val_ramp, highest_val_ramp, Time_ramp):
lineseg=segment()
lineseg.ramp_down(False, 127, 127, 20e-9)

#%%
#(full_memory, frequency, offset, amp):
sineseg=segment()
sineseg.sine(False, [10e6], 127, 64)

#%%
mem_seq4 = regpgenQC(None, premanip, lineseg, sineseg, postmanip, None, None, None)
mem_seq4.generate_mem_file()
mem_seq4.generate_rg_file()
