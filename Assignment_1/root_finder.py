def function(x):
    # Takes an argument, x, and returns the value of the function
    # we are trying to find the root of. 
    
    # This will only work for 1-D functions. 
    return x**2.0 - 4

def bisection(function, endpoint_a, endpoint_b, tolerance, maximum_iterations):
    
    if function(endpoint_a)*function(endpoint_b) > 0:
        print("You have not chosen points on both sides on the x-axis. \
              convergence will be impossible")
        exit()
        
    for i in range(0, maximum_iterations):
        if i == 0:
            a = endpoint_a
            b = endpoint_b
        
        c = (a+b)/2
        
            
        if function(a) * function(c) < 0:
            b = c
        if function(b) * function(c) < 0:
            a = c
            
        if abs(function(c)) < tolerance:
            print("We've reached our tolerance threshold of "+str(tolerance)+" in \
                  "+str(i)+" iterations")
            return(function(c),c)
        
        if i == maximum_iterations:
            print("We cannot find a root! Try different starting \
                  points.")
        
        
    
    
    
''' 
def newton_method(f,g):
    
    
def secant_method(f,g):
'''   
    
    