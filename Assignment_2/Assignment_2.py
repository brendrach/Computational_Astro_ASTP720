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


def freq(transition_start , transition_end):
    h = const.h
    freq = abs(( -13.6 * u.eV * transition_end**(-2) ) - ( -13.6 * u.eV * transition_start**(-2)))/h
    return freq


def J(temp , transition_start , transition_end):
    h = const.h
    c = const.c
    k = const.k_B
    freq_value = freq(transition_start , transition_end)
    freq_value = freq_value.to(1 / u.s)
    
    temp_const = 2 * h * freq_value**3 * c**(-2)
    temp_exp = h * freq_value / (k * temp)
    return temp_const / (np.exp(temp_exp) - 1)
    
def read_file():
    
    file = np.genfromtxt('A_coefficients.dat', delimiter=',')
    
    lower_coeffs = file[:,0]
    upper_coeffs = file[:,0]
    A_coeffs = file[:,2]
    
    A = matrix.Matrix([8,8])
    lower = matrix.Matrix([8,8])
    upper = matrix.Matrix([8,8])
    k = 0
    for i in range(0, 8):
        for j in range(i, 8):
            A.elems[i][j] = A_coeffs[k]
            lower.elems[i][j] = lower_coeffs[k]
            upper.elems[i][j] = upper_coeffs[k]
            k = k+1
            
    return lower.transpose(), upper, A

def get_lower(lower, upper):
    lower = read_file()[0]
    
    return lower.elems[lower][upper]

def get_upper(lower, upper):
    upper = read_file()[1]
    
    return upper.elems[lower][upper]

def get_Aval(lower, upper):
    Aval = read_file()[2]
    
    
    return Aval.elems[lower][upper], Aval.elems[upper][lower]
    
def calc_B(lower, upper):
    h = const.h
    c = const.c
    
    A_temp = get_Aval(lower, upper)[0]
    
    if A_temp == 0:
        A_temp = get_Aval(lower, upper)[1] * ( ((2 * lower ** 2) / (2 * upper ** 2)) ** 1)
    
    a = (c ** 2) / (2 * h * freq(lower , upper) ** 3)
    temp_B = a * A_temp
    return temp_B


def get_freq(l , up):
    '''
    finds the frequency for a particular line transition of hydrogen
    takes in two energy levels (integers)
    returns a frequency
    '''
    E = abs(( -13.6 * u.eV / (up ** 2) ) - ( -13.6 * u.eV / (l ** 2) ))
    return E / const.h
    
def densities(T):

    '''
    This function takes in a temperature T
    computes all of our number densities
    returns a Matrix object containing all of the number densities
    '''
    
    ###(3x3)
    b = np.matrix((0,0,1))
    A = np.zeros((3,3))
    
    A.elements[0][0] = ((calc_B(1 , 2) * J(T , 1 , 2) + calc_B(1 , 3) * J(T , 1 , 3))).value
    A.elements[0][1] = -1 * (get_Aval[0](2 , 1) + calc_B(2 , 1) * J(T , 2 , 1)).value
    A.elements[0][2] = -1 * (get_Aval[0](3 , 1) + calc_B(3 , 1) * J(T , 3 , 1)).value
    A.elements[1][1] = (calc_B(2 , 1) * J(T , 2 , 1) + get_Aval(2 , 1) + calc_B(2 , 3) * J(T , 2 , 3)).value
    A.elements[1][0] = -1 * (calc_B(1 , 2) * J(T  , 1 , 2)).value
    A.elements[1][2] = -1 * (calc_B(3 , 2) * J(T , 3  ,2) + calc_Aval(3 , 2)).value
    '''
    A.elements[2][0] = -1 * (get_B(1 , 3) * J(T , 1 , 3)).value
    A.elements[2][1] = -1 * (get_B(2 , 3) * J(T , 2 , 3)).value
    A.elements[2][2] = (get_B(3 , 1) * J(T , 3 , 1) + get_A(3 , 1) + get_A(3, 2) + get_B(3 , 2) * J(T , 3 , 2)).value
    '''
    A.elements[2][0] = 1
    A.elements[2][1] = 1
    A.elements[2][2] = 1
    x = matrix.solve_eq(A , b)
    return x


def nden(T):
    '''
    This function takes in a temperature T
    computes all of our number densities
    returns a Matrix object containing all of the number densities
    '''
    
    n_levels = 9
    
    A = matrix.Matrix((n_levels , n_levels))
    b = matrix.Matrix((n_levels , 1))
    b.elements[n_levels - 1][0] = 1
    
    for i in range(n_levels):
        ###These for loops set up the matrix A, which contains many einstein coefficieints
        ### A contains all of the coefficients for a system of equations for the various number densities
        
        for j in range(n_levels):
            if i == n_levels - 1:
                A.elements[i][j] = 1
                continue
            elm = 0
            if i == j:
                for k in range(n_levels):
                    if k == i: ###No transitions from one state to itself
                        continue
                    elm += (get_A(i + 1 , k + 1).value)
                    elm += (get_B(i + 1 , k + 1) * J(T , i + 1 , k + 1)).value
                A.elements[i][j] = elm
                
            else:
                #for k in range(n_levels):

                elm = (get_B(j + 1 , i + 1) * J(T , j + 1, i + 1) + get_A(j + 1, i + 1)).value
                A.elements[i][j] = -1 * elm
        

    ###Now we solve our matrix problem
    
    x = matrix.solve_eq(A , b)

    return x

###Runs code to produce all of our plots


