import astropy.units as u
import astropy.constants as const
import numpy as np
import matplotlib.pyplot as plt
import sys



def calc_density(pressure, stellar_type):
    
    type_id = id(stellar_type)
 
    if type_id == id("NS") or type_id == id('ns'):
        density = (pressure/(5.4*10**9))**(3/5)
    
    elif type_id == id('WD') or type_id == id('wd'):
        density = (pressure/1e13)**(3/5) * 2
        
    else: 
        sys.exit("This function only defines the density of white dwarfs (WD) \
                 and neutron stars (NS).")
        
    return density
        

def euler_method(func, y0, x, dx):    
    y = np.zeros(len(x))
    y[0] = y0
    
    for i in range(1, len(x)):
        y[i] = y[i-1] + dx * func(x[i-1], y[i-1])
            
    return y
        
def heuns_method(func, y0, x, dx):
    y = np.zeros(len(x))
    y[0] = y0
    
    for i in range(1, len(x)):
        func_val = func(x[i-1], y[i-1])
        y[i] = y[i-1] + dx/2 * (func_val + func(x[i], y[i] + dx*func_val))
            
    return y
        

def rk4_method(func, y0, x, dx):
    y = np.zeros(len(x))
    y[0] = y0
    
    for i in range(1, len(x)):
        k1 = dx * func(x[i-1], y[i-1])
        k2 = dx * func(x[i-1] + dx/2, y[i-1] + k1/2)
        k3 = dx * func(x[i-1] + dx/2, y[i-1] + k2/2)
        k4 = dx * func(x[i-1] + dx, y[i-1] + k3)
        y[i] = y[i-1] + 1/6*(k1 + 2*k2 + 2*k3 + k4)
       
    return y

    

def dm_dr(r, pressure, stellar_type):
    pi = np.pi 
    return 4 * pi * r**2 * calc_density(pressure, stellar_type)   

def dP_dr(r, pressure, M_enc, stellar_type):
    G = const.G.cgs    
    return (-G * M_enc * calc_density(pressure, stellar_type))/r**2 
    
    