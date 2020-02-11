import astropy.units as u
import astropy.constants as const
import numpy as np
import matplotlib.pyplot as plt
import numerical_calculus



def v_c(r, r_200, v_200, c):
    x = r / (v_200.value * u.kpc)
    numer = np.log(1 + c * x) - (c * x)/(1 + c * x)
    denom = x * (np.log(1 + c) - (c)/(1 + c))
    
    return np.sqrt(v_200**2 * numer/denom)


def M_enc(r, r_200, v_200, c): 
    return (r * v_c(r, r_200, v_200, c)**2 / const.G).to(u.solMass)


def density(r , r_200 , v_200 , c):
	
    r_s = r_200 / c
    return 1 / ((r / r_s) * (1 + r / r_s) ** 2)


def M(r , r_200 , v_200 , c):
	return 4 * np.pi * r**2 * density(r , r_200 , v_200 , c)


def problem_2():
    r_200 = 230 * u.kpc
    v_200 = 200 * u.km / (u.s)
    c = 15
    
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
    plt.xlabel("r (kpc)")
    plt.ylabel("M(r) M_sun")
    plt.show()
    
    plt.plot(x_grad , m_grad)
    plt.xlabel("r (kpc)")
    plt.ylabel("M'(r) M_sun")
    plt.show()
    
    plt.plot(np.linspace(0, 300, 301) , np.log10(m_enc))
    plt.xlabel("r (kpc)")
    plt.ylabel("M_enc (r) M_sun")
    plt.show()
    
    Mtot = M_enc(1e8 * u.kpc , r_200 , v_200 , c)
    print (np.log10(Mtot.value))
    
    print (M(5 * u.kpc , r_200 , v_200 , c))


problem_2()