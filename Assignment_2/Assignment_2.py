import astropy.units as u
import astropy.constants as const
import numpy as np
import matplotlib.pyplot as plt
import numerical_calculus



def v_c(r, r_200, v_200, c):
    '''
    Summary: 
    Calculates the circular velocity at a given radius.
    Parameters
    ----------
    r : float, the radius in kiloparsecs.
    r_200 : float, virial radius.
    v_200 : float, circular velocity at r_200.
    c : dimensionless concentration factor. 
    '''
    x = r / (v_200.value * u.kpc)
    numer = np.log(1 + c * x) - (c * x)/(1 + c * x)
    denom = x * (np.log(1 + c) - (c)/(1 + c))
    
    return np.sqrt(v_200**2 * numer/denom)


def M_enc(r, r_200, v_200, c): 
    '''
    Summary: 
    Calculates the mass interior to a radius in solar masses.
    Parameters
    ----------
    r : float, the radius in kiloparsecs.
    r_200 : float, virial radius.
    v_200 : float, circular velocity at r_200.
    c : dimensionless concentration factor. 
    '''
    return (r * v_c(r, r_200, v_200, c)**2 / const.G).to(u.solMass)


def density(r , r_200 , v_200 , c):
    '''
    Summary: 
    Calculates the density as a function of radius, r.
    Parameters
    ----------
    r : float, the radius in kiloparsecs.
    r_200 : float, virial radius.
    v_200 : float, circular velocity at r_200.
    c : dimensionless concentration factor. 
    '''
    
    r_s = r_200 / c
    return 1 / ((r / r_s) * (1 + r / r_s) ** 2)


def M(r , r_200 , v_200 , c):
    '''
    Summary: 
    Calculates the mass as a function of radius, r.
    Parameters
    ----------
    r : float, the radius in kiloparsecs.
    r_200 : float, virial radius.
    v_200 : float, circular velocity at r_200.
    c : dimensionless concentration factor. 
    '''
    
    return 4 * np.pi * r**2 * density(r , r_200 , v_200 , c)


def problem_2():
    '''
    Summary: 
    Produces plots for problem 2.
    '''
    
    r_200 = 230 * u.kpc
    v_200 = 250 * u.km / (u.s)
    c = 8
    
    m_r = []
    m_enc = []
    
    def M_func(i):
        M_func = M(i * u.kpc, r_200 , v_200 , c).value
        return M_func
    
    x_grad, m_grad = numerical_calculus.symmetric_derivative(M_func, 0, 300, 301)
    
    for i in np.linspace(0, 300, 301):
            m_r.append(M(i * u.kpc , r_200 , v_200 , c).value)
            m_enc.append(M_enc(i * u.kpc , r_200 , v_200 , c).value)
        
    plt.plot(np.linspace(0, 300, 301) , np.log10(m_r))
    plt.title(r"$log(M(r)) \ vs \ r, \ c \ = \ $" + str(c) + r"$, \ v_{200} \ = \ $" + str(v_200))
    plt.xlabel(r"$r \ (kpc)$")
    plt.ylabel(r"$\log(M(r)) \ M_\odot$")
    plt.savefig("v200_"+str(v_200.value)+"_c_"+str(c)+"_m.png")
    plt.show()
    
    plt.plot(x_grad , m_grad)
    plt.title(r"$log(\frac{d M(r)}{dr}) \ vs \ r, \ c \ = \ $" + str(c) + r"$, \ v_{200} \ = \ $" + str(v_200))
    plt.xlabel(r"$r \ (kpc)$")
    plt.ylabel(r"$log(\frac{d M(r)}{dr}) \ M_\odot$")
    plt.savefig("v200_"+str(v_200.value)+"_c_"+str(c)+"_dmdr.png")
    plt.show()
    
    plt.plot(np.linspace(0, 300, 301) , np.log10(m_enc))
    plt.title(r"$log(M_{enc}(r)) \ vs \ r, \ c=$" + str(c) + r"$ \ , \ v_{200} \ = \ $" + str(v_200))
    plt.xlabel(r"$r \ (kpc)$")
    plt.ylabel(r"$\log(M_{enc}(r)) \ M_\odot$")
    plt.savefig("v200_"+str(v_200.value)+"_c_"+str(c)+"_menc.png")
    plt.show()
    
    M_galaxy = M_enc(15 * u.kpc , r_200 , v_200 , c)
    M_DM = M_enc(100 * u.kpc, r_200, v_200, c)
    print (r"Mass of Dark Matter Halo = %.4g M_sun" % (M_DM - M_galaxy).value )
    
problem_2()


def freq(upper , lower):
    '''
    Summary: 
    Computes the frequency of a photon emitted when moving
    from energy state upper to energy state lower.
    Parameters
    ----------
    upper : int, energy state.
    lower : int, energy state.
    '''
    
    ## Mimics an in text Equation directly under Eq. 8 in the assignment.
    h = const.h
    freq = abs(( -13.6 * u.eV * (upper)**(-2) ) - ( -13.6 * u.eV * (lower)**(-2)))/h
    return freq


def J(temp , upper , lower):
    '''
    Summary: 
    Calculates the spherically averaged mean intensity from a photon
    emitted from energy state upper and falls to energy state lower.
    Parameters
    ----------
    temp : float, the temperature of the gas.
    upper : int, energy state.
    lower : int, energy state.
    '''
    
    h = const.h
    c = const.c
    k = const.k_B
    
    ## Calculate the frequency
    freq_value = freq(upper , lower).to(1 / u.s)
    
    ## Calculate the averaged mean intensity.
    ## Modelled after Eq. 8 in the assignment.
    temp_const = 2 * h * freq_value**3 * c**(-2)
    temp_exp = (h * freq_value / (k * temp))
    temp_exp = (h * freq_value / (k * temp))
    return temp_const / (np.exp(temp_exp) - 1)
    
def read_file():
    '''
    Summary: 
    Reads in the file containing the Einstein Coeffs.
    Returns the upper, lower, and A matrices.
    '''
    
    file = np.genfromtxt('A_coefficients.dat', delimiter=',')
    
    lower_coeffs = file[:,0]
    upper_coeffs = file[:,0]
    A_coeffs = file[:,2]
    
    A = matrix.Matrix([8,8])
    lower = matrix.Matrix([8,8])
    upper = matrix.Matrix([8,8])
    
    ## Loop through the matrix and store the data triangularly.
    k = 0
    for i in range(0, 8):
        for j in range(i, 8):
            A.elems[i][j] = A_coeffs[k]
            lower.elems[i][j] = lower_coeffs[k]
            upper.elems[i][j] = upper_coeffs[k]
            k = k+1
            
    return lower.transpose(), upper, A


def get_lower(lower, upper):
    '''
    Summary: 
    Returns the value of the lower matrix at a given set of indices.
    ----------
    upper : int, energy state.
    lower : int, energy state.
    '''
    
    lower = read_file()[0]
    
    return lower.elems[lower][upper]


def get_upper(lower, upper):
    '''
    Summary: 
    Returns the value of the upper matrix at a given set of indices.
    ----------
    upper : int, energy state.
    lower : int, energy state.
    '''
    
    upper = read_file()[1]
    
    return upper.elems[lower][upper]


def get_Aval(lower, upper):
    '''
    Summary: 
    Returns the value of the Einstein coeffs at a given set of indices.
    ----------
    upper : int, energy state.
    lower : int, energy state.
    '''
    
    Aval = read_file()[2]
    
    return Aval.elems[lower][upper], Aval.elems[upper][lower]
    

def get_Bval(lower, upper):
    '''
    Summary: 
    Computes the value of the Einstein Coeff, B.
    ----------
    upper : int, energy state.
    lower : int, energy state.
    '''
    
    h = const.h
    c = const.c
    
    A_temp = get_Aval(lower, upper)[0]
    
    if A_temp == 0:
        A_temp = get_Aval(lower, upper)[1] * ( ((2 * lower ** 2) / (2 * upper ** 2)) ** 1)
    
    a = (c ** 2) / (2 * h * freq(lower , upper) ** 3)
    temp_B = a * A_temp
    return temp_B


