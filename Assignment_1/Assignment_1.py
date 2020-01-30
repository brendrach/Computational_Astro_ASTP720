import numpy as np
import interp as pol
import root_finder as rf
import matplotlib.pyplot as plt
import astropy.units as u
from astropy import constants as const



def problem_3():
    '''
    Summary:
    Solves problem 3 on Assignment 1. I numerically calculate the
    FWHM for the pseudo isothermal sphere model. I produce a plot 
    showing the number of iterations as a function of the threshold. 
    '''
    
    ## Equation 14 in the assignment set equal to 0 for root finding.
    def Ne_Half(x_p):
        return (1 + x_p**2.0)**(-1/2) - 1/2

    ## Equation 15 in the assignment set equal to 0 for root finding.
    def Ne_Half_Derivative(x_p):
        return -(1 + x_p**2.0)**(-3/2) * x_p
    
    ## Define the number of iterations before we stop trying to converge.
    max_iters = 4000
    
    ## Initialize our tolerance array. We will index this array to 
    ## test different tolerances. 
    tol = [1]
    for i in range(0, 16):
        tol.append(tol[i]*0.1)
        
    ## Initialize arrays that will hold the interations-until-convergence values.
    con_bi = []
    con_ne = []
    con_se = []
    
    ## Run the root finding algorithms for every tolerance in the tol array.
    ## For documentation on the functions used below, see root_finder.py.
    for i in range(0, len(tol)):
    
        bi_zero, bi_iter = rf.bisection(Ne_Half, -5, 0, tol[i],  max_iters)
        ne_zero, ne_iter = rf.newton_method(Ne_Half, Ne_Half_Derivative, 5, tol[i], max_iters)
        se_zero, se_iter = rf.secant_method(Ne_Half, -5, 0, tol[i], max_iters)
        
        con_bi.append(bi_iter)
        con_ne.append(ne_iter)
        con_se.append(se_iter)
        
    
    ## Plot the data.    
    plt.semilogx(tol, con_bi, 'r', label = 'Bisection Convergence')
    plt.semilogx(tol, con_ne, 'b', label = 'Newtons Method Convergence')
    plt.semilogx(tol, con_se, 'g', label = 'Secant Method Convergence')
    plt.legend()
    plt.xlabel("log(threshold)")
    plt.ylabel("Iterations")
    plt.savefig("threshold_vs_iterations.png")
    plt.show()
    

def problem_4():
    '''
    Summary:
    Solves problem 4 on Assignment 1. I numerically raytrace the particles
    around a Gaussian Lens.
    '''
    
    ## Define x_prime as it is found in Eq. 11. 
    ## I believe there is a typo in this equation.
    ## I suspect it should read lambda^2 * r_e instead of 
    ## lambda * r_e^2. 
    def x_prime(lam, r_e, N0, D, a, x):
        xp = x * (1 + lam**2.0 * r_e * N0 * D / (np.pi * a**2.0) * np.exp(-(x/a)**2.0))
        return xp
    
    ## Define constants using Astropy
    x = np.linspace(-10,10,500) * u.au
    lam = 21 * u.cm
    r_e = 2.818e-15 * u.m
    N0 = 0.01 * u.pc / u.cm ** 3
    D = 1 * u.kpc
    a = 1 * u.au
    
    ## Calculate x_prime for each value of x. Plot each line on the
    ## same plot to illustrate the ray-tracing.
    for i in range(0, len(x)):
        xp = x_prime(lam, r_e, N0, D, a, x[i])
        plt.plot([x[i].value , xp.value], [1,0] , color = 'r' , linewidth = .3)
        
    plt.ylabel("Distance (Kpc)")
    plt.xlabel("x (AU)")
    plt.savefig("gaussian_lense.png")
    plt.show()
    
    
def problem_5():
    '''
    Summary:
    Solves problem 5 on Assignment 1. I numerically raytrace the particles
    around a pseudo-isothermal sphere.
    '''
    
    ## Define x_prime by combining Eq. 15, 6, and 7 in that order. 
    ## I believe there is a typo in Eq. 15. There should be an N0 in 
    ## the numerator. I've included that here and think I got the right answer. 
    ## Also, the units won't work out without including it.
    def x_prime(lam, r_c, r_e, N0, D, a, x):
        dN_dx = -x * N0 / (r_c**2.0 * (1 + (x/r_c)**2.0)**(3/2))
        theta_r = lam**2.0 * r_e / (2 * np.pi) * dN_dx
        xp = x - theta_r * D
        return xp
    
    ## Define constants using Astropy
    x = np.linspace(-5,5,500) * u.au
    lam = 21 * u.cm
    r_e = 2.818e-15 * u.m
    N0 = 0.01 * u.pc / u.cm ** 3
    D = 1 * u.kpc
    a = 1 * u.au
    r_c = 1 * u.au
    
    ## Calculate x_prime for each value of x. Plot each line on the
    ## same plot to illustrate the ray-tracing.
    for i in range(0, len(x)):
        xp = x_prime(lam, r_c, r_e, N0, D, a, x[i])
        plt.plot([x[i].value , xp.value], [1,0] , color = 'g' , linewidth = .3)
    
    plt.ylabel("Distance (Kpc)")
    plt.xlabel("x (AU)")
    plt.savefig("pseudo_isothermal_sphere.png")
    plt.show()


        

    
    
    
    
    
problem_3()
problem_4()
problem_5()

