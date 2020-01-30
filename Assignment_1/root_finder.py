def function_quadratic(x):
    '''
    Summary:
    Returns the value of a given quadratic function at a specific point, x.
    
    Parameters
    ----------
    x : an integer, x value. f(x) value will be returned. 
    '''
    
    ## Define the function you want to find the root of.
    func = x**2.0 - 4 
    
    ## Return the value of the function at point, x. 
    return func

def function_quadratic_derivative(x):
    '''
    Summary:
    Returns the value of the derivative of a given quadratic function at a specific point, x.
    
    Parameters
    ----------
    x : an integer, x value. f'(x) value will be returned. 
    '''
    
    ## Define the function you want to find the root of.
    func = 2*x 
    
    ## Return the value of the function at point, x. 
    return func

def bisection(function, endpoint_a, endpoint_b, tolerance, maximum_iterations):
    '''
    Summary:
    Finds the root of a function using the bisection method.
    
    Parameters
    ----------
    function : a function defined elsewhere.
    endpoint_a : float, The starting point from one side of the y-axis.
    endpoint_b : float, The starting point from the other side of the y-axis.
    tolerance : float, How close to the root you want to converge to before stopping.
    maximum_iterations : int, How long we will try to converge before exiting. 
    '''
    
    ## Check to see if f(a) and f(b) are both positive or both negative. 
    ## If so, our root cannot be between the values a and b. Therefore,
    ## we will exit and pick a new point!
    if function(endpoint_a)*function(endpoint_b) > 0:
        print("You have not chosen points on both sides on the x-axis. \
              convergence will be impossible")
        exit()
        
    ## Loop over the maximum number of iterations.
    for i in range(0, maximum_iterations):
        
        ## Set the initial endpoints.
        if i == 0:
            a = endpoint_a
            b = endpoint_b
        
        ## Calculate the midpoint.
        c = (a+b)/2
        
        ## If f(a) * f(c) is less than zero, that means that f(c) lies
        ## on the same side of the y-axis as f(b).
        ## Therefore, c becomes the new b. 
        if function(a) * function(c) < 0:
            b = c
            
        ## Identical reason as above. 
        if function(b) * function(c) < 0:
            a = c
            
        ## If the value of f(c) is below our tolerance, we have found our root.
        if abs(function(c)) < tolerance:
            print("We've reached our tolerance threshold of "+str(tolerance)+" in \
                  "+str(i)+" iterations.")
            print("The value of the root is " + str(c) + ".")
            return(c, i)
        
        ## If we reach our maximum number of iterations without converging, 
        ## we must've chosen a bad starting point.
        if i == maximum_iterations:
            print("We cannot find a root! Try different starting \
                  points.")
        
        

def newton_method(function, function_derivative, starting_point, tolerance, maximum_iterations):
    '''
    Summary:
    Finds the root of a function using Newton's Method
    
    Parameters
    ----------
    function : a function defined elsewhere.
    function_derivative : the derivative of function
    starting_point : float, The value where we will begin our convergence.
    tolerance : float, How close to the root you want to converge to before stopping.
    maximum_iterations : int, How long we will try to converge before exiting. 
    '''
    
    ## Initialize a temporary array to store our values as we converge.
    ## Append our starting point.
    temp_x = []
    temp_x.append(starting_point)
    
    ## Loop over the maximum number of iterations.
    for i in range(0, maximum_iterations):
        
        ## Compute x_(i+1) using Newton's Method
        x_i_1 = temp_x[i] - function(temp_x[i])/function_derivative(temp_x[i])
        
        ## Append the value to the temp_x array.
        temp_x.append(x_i_1)
        
        ## If the difference between successive values of x are below the tolerance 
        ## than we have likely reached our root. Convergence will continue but 
        ## at a precision beyond what we need.
        if abs(temp_x[i+1]-temp_x[i]) < tolerance:
            print("We've reached our tolerance threshold of "+str(tolerance)+" in \
                  "+str(i)+" iterations.")
            print("The value of the root is " + str(temp_x[i+1]) +".")
            return(temp_x[i+1], i)
        
        ## If we reach our maximum number of iterations without converging, 
        ## we must've chosen a bad starting point or something else must be wrong!
        if i == maximum_iterations:
            print("We cannot find a root! Try different starting \
                  points.")
    
    
    
def secant_method(function, starting_point, secondary_point, tolerance, maximum_iterations):
    '''
    Summary:
    Finds the root of a function using the secant method.
    
    Parameters
    ----------
    function : a function defined elsewhere.
    starting_point : float, The value where we will begin our convergence.
    tolerance : float, How close to the root you want to converge to before stopping.
    maximum_iterations : int, How long we will try to converge before exiting. 
    '''
    
    ## Initialize a temporary array to store our values as we converge.
    ## Append our starting point.
    temp_x = []
    temp_x.append(starting_point)
    temp_x.append(secondary_point)
    
    ## Loop over the maximum number of iterations.
    for i in range(1, maximum_iterations):
        
        
        ## Compute x_(i+1) using the secant method
        ## First, I will define the denominator as it is written in Eq. 22 of 
        ## Michael Lam's Root Finding Notes - RIT MyCourses access required.
        ## https://mycourses.rit.edu/d2l/le/content/819279/viewContent/6251151/View
        denom_x_i_1 = (function(temp_x[i]) - function(temp_x[i-1]))/(temp_x[i] - temp_x[i-1])
        x_i_1 = temp_x[i] - function(temp_x[i])/denom_x_i_1
        
        ## Append the value to the temp_x array.
        temp_x.append(x_i_1)
        
        ## If the difference between successive values of x are below the tolerance 
        ## than we have likely reached our root. Convergence will continue but 
        ## at a precision beyond what we need.
        if abs(temp_x[i+1]-temp_x[i]) < tolerance:
            print("We've reached our tolerance threshold of "+str(tolerance)+" in \
                  "+str(i)+" iterations.")
            print("The value of the root is " + str(temp_x[i+1]) +".")
            return(temp_x[i+1], i)
        
        ## If we reach our maximum number of iterations without converging, 
        ## we must've chosen a bad starting point or something else must be wrong!
        if i == maximum_iterations:
            print("We cannot find a root! Try different starting \
                  points.")
    
  
    
    