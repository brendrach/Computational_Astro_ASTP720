import numpy as np


def symmetric_derivative(function, starting_x, ending_x, num_points):
    '''
    Summary:
    Finds the symmetric derivative of a function.
    
    Parameters
    ----------
    function : a function you want to find the derivative of.
    starting_x : int OR float, The beginning x-value in your dataset.
    ending_x : int OR float, The last x-value in your dataset.
    num_points : int, the number of points you want to discretize the x-range into. 
    '''
    
    ## Calculate the list of x datapoints.
    x = np.linspace(starting_x, ending_x, num_points)
    
    ## Calculate the spacial step size. 
    delta_x = (ending_x - starting_x)/num_points
    
    ## Initialize the array that will hold the values of the derivative.
    y_prime = np.zeros(len(x)-1)
    
    ## Calculate the value of the derivative at each point using the midpoint formula.
    ## Found in the class notes titled "Numerical Calculus" -- Eq. 6.
    for i in range(1, len(x)-1):
        y_prime[i] = (function(x[i+1]) - function(x[i-1]))/(2*delta_x)
        
    ## Return the values of the derivative and the data points we calculated them at. 
    return x[:-1], y_prime
    

def midpoint_integration(function, starting_x, ending_x, num_points):
    '''
    Summary:
    Finds the integral of the function using the midpoint rule.
    
    Parameters
    ----------
    function : a function you want to find the derivative of.
    starting_x : int OR float, The beginning x-value in your dataset.
    ending_x : int OR float, The last x-value in your dataset.
    num_points : int, the number of points you want to discretize the x-range into. 
    '''
    
    ## Calculate the list of x datapoints.
    x = np.linspace(starting_x, ending_x, num_points)
    
    ## Calculate the spacial step size. 
    delta_x = (ending_x - starting_x)/num_points
    
    ## Initialize the array that will hold the values of the integral.
    y = np.zeros(len(x)-1)
    
    ## Calculate the value of the integral at each point using the midpoint formula.
    ## Found in the class notes titled "Numerical Calculus" -- Eq. 6.
    for i in range(1, len(x)-1):
        x_bar = (x[i-1] + x[i+1])/2
        y[i] = (function(x_bar))*delta_x
      
    y = np.sum(y)
    ## Return the values of the integral and the data points we calculated them at. 
    return y
    
    

def trapezoid_integration(function, starting_x, ending_x, num_points):
    '''
    Summary:
    Finds the integral of the function using the trapezoid rule.
    
    Parameters
    ----------
    function : a function you want to find the derivative of.
    starting_x : int OR float, The beginning x-value in your dataset.
    ending_x : int OR float, The last x-value in your dataset.
    num_points : int, the number of points you want to discretize the x-range into. 
    '''
    
    ## Calculate the list of x datapoints.
    x = np.linspace(starting_x, ending_x, num_points)
    
    ## Calculate the spacial step size. 
    delta_x = (ending_x - starting_x)/num_points
    
    ## Initialize the array that will hold the values of the integral.
    y = np.zeros(len(x)-1)
    
    ## Calculate the value of the integral at each point using the trapezoid formula.
    ## Found in the class notes titled "Numerical Calculus" -- Eq. 35.
    for i in range(0, len(x)-1):
        y[i] = (function(x[i]) + function(x[i+1]))*delta_x/2
      
    y = np.sum(y)
    ## Return the values of the integral and the data points we calculated them at. 
    return y
    
 
def simpsons_integration(function, starting_x, ending_x, num_points):
    '''
    Summary:
    Finds the integral of the function using the simpsons rule.
    
    Parameters
    ----------
    function : a function you want to find the derivative of.
    starting_x : int OR float, The beginning x-value in your dataset.
    ending_x : int OR float, The last x-value in your dataset.
    num_points : int, the number of points you want to discretize the x-range into. 
    '''
    
    ## Calculate the list of x datapoints.
    x = np.linspace(starting_x, ending_x, num_points)
    
    ## Calculate the spacial step size. 
    delta_x = (ending_x - starting_x)/num_points
    
    ## Initialize the array that will hold the values of the integral.
    y = np.zeros(len(x)-2)
    
    ## Calculate the value of the integral at each point using the trapezoid formula.
    ## Found in the class notes titled "Numerical Calculus" -- Eq. 35.
    for i in range(0, int(len(x)/2)-1):
        print(i)
        y[i] = (function(x[2*i]) + 4*function(x[2*i+1])+ function(x[2*i+2]))*delta_x/3
      
    y = np.sum(y)
    ## Return the values of the integral and the data points we calculated them at. 
    return y
    
