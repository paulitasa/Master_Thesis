# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 14:09:21 2023

@author: Paula Sánchez Almagro
"""

import numpy as np

def write_memfile(filename, header, y, mem_end, row_size):
    z=len(y)

    for i in range (0, z):
        if (y[i]<0):
            y[i]=0
    
    y=np.asarray(y, dtype='int')
    
    if(z<(row_size-1)):
        raise ValueError('data length < 63 Bytes');
    elif z>mem_end :
        raise ValueError('data length > 32 kBytes');

    fid=open(filename, "w");
    fid.write(header + "\n");
    
    
    ind =1;
    count=0;
    fid.write("   "+str(0)+" ")
    for i in range (0, z):
        if(count>31):
            fid.write("\n")
            if(ind<10):
                fid.write("   ")
            elif(ind<100):
                fid.write("  ")
            elif(ind<1000):    
                fid.write(" ")
            fid.write(str(ind) + " ")
            count=1;
            ind=ind+1;
        else:    
            count=count+1;
            
        if(y[i]<16):
            fid.write("0"+str(hex(y[i]).removeprefix('0x'))) 
        else:
            fid.write(str(hex(y[i]).removeprefix('0x'))) 

    
    fid.close();
    
    
def write_memfile_part(filename, header, y, mem_end, row_size, length):
    y=y[0:((length+1)*32)]
    z=len(y)

    for i in range (0, z):
        if (y[i]<0):
            y[i]=0
    
    y=np.asarray(y, dtype='int')
    
    if(z<(row_size-1)):
        raise ValueError('data length < 63 Bytes');
    elif z>mem_end :
        raise ValueError('data length > 32 kBytes');
        

    fid=open(filename, "w");
    fid.write(header + "\n");
    
    
    ind =1;
    count=0;
    fid.write("   "+str(0)+" ")
    for i in range (0, z):
        if(count>31):
            fid.write("\n")
            if(ind<10):
                fid.write("   ")
            elif(ind<100):
                fid.write("  ")
            elif(ind<1000):    
                fid.write(" ")
            fid.write(str(ind) + " ")
            count=1;
            ind=ind+1;
        else:    
            count=count+1;
            
        if(y[i]<16):
            fid.write("0"+str(hex(y[i]).removeprefix('0x'))) 
        else:
            fid.write(str(hex(y[i]).removeprefix('0x'))) 

    
    fid.close();


   
    
