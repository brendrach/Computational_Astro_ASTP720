import numpy as np
import emcee

def least_squares(x, y):
    '''
    Summary:
    Solves the equation X.T * Y = X.T * X * theta

    Parameters
    ----------
    x : design matrix. 
    y : observation matrix.
    ''' 
    ## Add a row of ones to account for
    ## the unweighted term.
    x = np.vstack([x, np.ones(len(x))]).T
    
    ## Compute the RHS
    XT = np.transpose(x)
    XTX = np.matmul(XT, x)
    
    ## Compute the LHS
    XTY = np.matmul(XT, y)
    
    ## Solve for theta.
    coefs = np.linalg.solve(XTX,XTY)
     
    return coefs

def lnprior(theta):
    '''
    Summary:
    Define the prior distribution.

    Parameters
    ----------
    theta : Contains the two fitting parameters - slope and y_int
    ''' 
    
    ## The parameters are stored as a vector of values, so unpack them
    slope, y_int = theta
    
    ## Keep the MCMC from wandering outside of this range.
    if slope > 150 or slope < 0:
        return -np.inf
    if y_int > 10000 or y_int < 0:
        return -np.inf
    return 0


def lnlike(theta, x, y, yerr):
    '''
    Summary:
    Define the log likelihood.

    Parameters
    ----------
    theta : Contains the two fitting parameters - slope and y_int
    x : the x-data that we will use to test the goodness of fit.
    y : the y-data that we will use to test the goodness of fit.
    yerr : the error on the distance measurements.
    ''' 
    
    ## Unpack the Parameters
    slope, y_int = theta
    
    ## Define the linear model.
    model = slope*x + y_int
    
    ## Define a "denominator" term for ease.
    denom = np.sqrt(2*np.pi*yerr)
    
    ## compute the log likelihood. 
    lp = np.sum((-(y - model)**2.0/(2 * yerr)) - np.log(denom))
    return lp

def lnprob(theta, x, y, yerr):
    '''
    Summary:
    Define the log likelihood.

    Parameters
    ----------
    theta : Contains the two fitting parameters - slope and y_int
    x : the x-data that we will use to test the goodness of fit.
    y : the y-data that we will use to test the goodness of fit.
    yerr : the error on the distance measurements.
    ''' 
    
    ## Keep us bounded by our prior range.
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    
    ## Compute the log probability.
    return lp + lnlike(theta, x, y, yerr)


def MCMC(Nwalker, Ndim, iterations, initial_guess, x, y, yerr):
    '''
    Summary:
    Define the log likelihood.

    Parameters
    ----------
    Nwalker : The number of walkers used by the MCMC to explore the
              paramater space.
    Ndim : The number of parameters we are fitting for. This could 
           probably have been hardcoded to be 2 for this specific 
           use case but I've kept it tweakable.
    iterations : The number of steps to take in our MCMC.
    initial_guess : An array of starting values for each parameter.
    theta : Contains the two fitting parameters - slope and y_int
    x : the x-data that we will use to test the goodness of fit.
    y : the y-data that we will use to test the goodness of fit.
    yerr : the error on the distance measurements.
    ''' 
    ## Initialize a walker at a random state that is slightly
    ## perturbed from our initial guess.
    p0 = [initial_guess+1.10*np.random.randn(Ndim) for i in range(Nwalker)]
    
    ## Initialize the sampler.
    sampler = emcee.EnsembleSampler(Nwalker,Ndim,lnprob, args=(x, y, yerr))
    
    ## Perform the MCMC.
    pos,prob,state = sampler.run_mcmc(p0, iterations, progress = True)
    
    return sampler, pos, prob, state