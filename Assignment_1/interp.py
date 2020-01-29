import numpy as np

def piecewise_linear_interp(x, y):
    '''
    Summary:
    Given a set of x and y data points, this function returns a function that will interpolate
    to any given x data point. It is built to handle either a single x value as input
    or a set of them (an array).
    
    Parameters
    ----------
    x : a list of x coordinate data
    y : a list of y coordinate data
    '''
    
    ## The slope of the line connecting neighboring data points
    ## in the input dataset.
    dy_dx = np.zeros(len(x))
    
    for i in range(1, len(x)):
        dy_dx[i] = (y[i] - y[i-1])/(x[i] - x[i-1])
        
        
    
    def interpolated_func(xp):
        '''
        Summary:
        Interpolates to an input x-value, xp, and returns the value. 
        Works with single data points or arrays.
        
        Parameters
        ----------
        xp : int or list, datapoint(s) that we will interpolate to. 
        '''


        ## Check if the input, xp, is a single data point. 
        if type(xp) == int:
            ## Make sure xp is in the valid domain of interpolation
            if xp > x[-1] or xp < x[0]:
                print("Your value, " + str(xp) + "is outside the range of x.")
                exit()
            
            ## Compute the line connecting the neighboring points and calculate
            ## y(xp) along that line.
            for i in range(len(x)):
                if xp == x[i]:
                    y_tmp = y[i]
                if xp > x[i] and xp < x[i+1]:
                    y_tmp = y[i] + dy_dx[i] * (xp - x[i])
                        
            return y_tmp
        
        ## Handle xp differently if the input is an array.
        else:
            y_tmp = []
            for j in range(0, len(xp)):
                ## Make sure xp is in the valid domain of interpolation
                if xp[j] > x[-1] or xp[j] < x[0]:
                    print("Your value, " + str(xp[j]) + "is outside the range of x.")
                    exit()
            
                ## Compute the line connecting the neighboring points and calculate
                ## y(xp) along that line.
                for i in range(len(x)):
                    if xp[j] == x[i]:
                        y_tmp.append(y[i])
                    if xp[j] > x[i] and xp[j] < x[i+1]:
                        y_tmp.append(y[i] + dy_dx[i] * (xp[j] - x[i]))
                        
            return y_tmp
                
    return interpolated_func
        