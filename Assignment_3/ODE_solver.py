import astropy.units as u
import astropy.constants as const
import numpy as np
import matplotlib.pyplot as plt
import sys
        

def euler_method(func, y0, x, dx):  
    '''
    Summary:
    Numerically integrates a given function using Euler's method.  
    
    Parameters
    ----------
    func :  function, a function with 1 independent variable, x, 
            and n-number of dependent variables y[i].
    y0 :    int, float, or list, the initial values of the depedent 
            variables, y[i].
    x :     list, a list of values (independent variables) that
            that we will solve the integral at.
    dx :    int or float, the step size used in numerical integration.
        
    '''
    
    ## Check if y0 is a number (int or float) and not a list.
    ## This tells us if we are solving an uncoupled system or
    ## a system with y[i] coupled variables.
    if isinstance(y0, int) or isinstance(y0, float):
        
        ## If it is a number (int or float) we convert it to a
        ## 1D list so the solver inteprets it correctly. 
        y0 = [y0]
        
    ## Initialize our solution array to have as many rows as y0 
    ## and as many columns as x.
    y = np.zeros((len(x), len(y0[:])))
        
    ## Set the initial values of the solution matrix.
    for i in range(len(y0)):
        y[0,i] = y0[i]
    
    ## Perform Euler's method to numerically integrate a function, func.
    ## This is modelled after Eq. 4 of Michael Lam's ODE Notes.
    for i in range(1, len(x)):
        y[i,:] = y[i-1,:] + np.multiply(dx, func(x[i-1], y[i-1,:]))
         
    return y
        
def heuns_method(func, y0, x, dx):
    '''
    Summary:
    Numerically integrates a given function using Heuns's method.  
    
    Parameters
    ----------
    func :  function, a function with 1 independent variable, x, 
            and n-number of dependent variables y[i].
    y0 :    int, float, or list, the initial values of the depedent 
            variables, y[i].
    x :     list, a list of values (independent variables) that
            that we will solve the integral at.
    dx :    int or float, the step size used in numerical integration.
    '''
    
    ## Check if y0 is a number (int or float) and not a list.
    ## This tells us if we are solving an uncoupled system or
    ## a system with y[i] coupled variables.
    if isinstance(y0, int) or isinstance(y0, float):
        
        ## If it is a number (int or float) we convert it to a
        ## 1D list so the solver inteprets it correctly. 
        y0 = [y0]
        
    ## Initialize our solution array to have as many rows as y0 
    ## and as many columns as x.
    y = np.zeros((len(x), len(y0[:])))
    
    ## Set the initial values of the solution matrix.
    for i in range(len(y0)):
        y[0,i] = y0[i]
    
    ## Perform Heuns's method to numerically integrate a function, func.
    ## This is modelled after Eq. 8 of Michael Lam's ODE Notes.
    for i in range(1, len(x)):
        
        ## Compute value of our function at a given point.
        func_val = func(x[i-1], y[i-1,:])
        
        ## Compute the bracketed term in Eq. 8 of Michael Lam's ODE notes.
        y_curl = np.add(func_val, func(x[i], np.add(y[i-1,:], np.multiply(dx, func_val))))
        
        ## Compute y[i] from y_curl.
        y[i,:] = y[i-1,:] + np.multiply(dx/2, y_curl)
        
        ## Notes, this could all be done in one line as it is in
        ## Eq. 8. But, that would look ugly because we have to 
        ## do matrix operations with the use of numpy. Breaking
        ## up the operations is much cleaner.
            
    return y
        

def rk4_method(func, y0, x, dx):
    '''
    Summary:
    Numerically integrates a given function using RK4.  
    
    Parameters
    ----------
    func :  function, a function with 1 independent variable, x, 
            and n-number of dependent variables y[i].
    y0 :    int, float, or list, the initial values of the depedent 
            variables, y[i].
    x :     list, a list of values (independent variables) that
            that we will solve the integral at.
    dx :    int or float, the step size used in numerical integration.
    '''  
    
    ## Check if y0 is a number (int or float) and not a list.
    ## This tells us if we are solving an uncoupled system or
    ## a system with y[i] coupled variables.
    if isinstance(y0, int) or isinstance(y0, float):
        
        ## If it is a number (int or float) we convert it to a
        ## 1D list so the solver inteprets it correctly. 
        y0 = [y0]
        
    ## Initialize our solution array to have as many rows as y0 
    ## and as many columns as x.
    y = np.zeros((len(x), len(y0[:])))
        
    ## Set the initial values of the solution matrix.
    for i in range(len(y0)):
        y[0,i] = y0[i]
            
    ## Perform RK4 to numerically integrate a function, func.
    ## This is modelled after the set of Eq's at the bottom
    ## of page 4 in Michael Lam's ODE notes.
    for i in range(1, len(x)):
        
        ## Compute k1, k2, k3, k4.
        k1 = np.multiply(dx, func(x[i-1], y[i-1,:]))
        k2 = np.multiply(dx, func(x[i-1] + dx/2, y[i-1,:] + k1/2))
        k3 = np.multiply(dx, func(x[i-1] + dx/2, y[i-1,:] + k2/2))
        k4 = np.multiply(dx, func(x[i-1] + dx, y[i-1, :] + k3))
        
        ## Calculate y[i].
        y[i, :] = y[i-1, :] + 1/6*(k1 + 2*k2 + 2*k3 + k4)

    return y

    