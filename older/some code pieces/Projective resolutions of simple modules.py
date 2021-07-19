######################################################################
# Projective resolutions of simple modules


def tilde(w):
    return w.inverse()*w0


def hat(w):
    return w0*w.inverse()


def proj_resolution_param_L(w):
    '''Returns the parameters of the projective resolution of L_w,
    as dictionary where keys are homological degrees, and values are again dictionaries,
    with keys elements in W and values are multiplicities.
    As the resolution is linear, the shifts are not recorded.'''
    
    P = char_P( tilde(w) )
    P_res = {}
    for i in P.component:
        P_res[i] = {}
        for x in P.component[i]:
            P_res[i][hat(x)] = P.component[i][x]
    return P_res



def euler_char(seq):
    '''Returns the Euler characteristics of the sequence, i.e. the alternating sum of terms.'''
    
    suma = char_0()
    for i in range(len(seq)):
        suma += (-1)^i * seq[i]
    return suma


def Pv(w):
    """Returns the KL polynomial p_{e,w}, but in the 'v'-normalization. See [KMM: Bigrassmannians and Vermas].
    Not a very efficient method, though."""
    
    dictM = dict_mult(char_Delta(e),w)
    #dictionary with items (k:m), where m is the multiplicity of w in k-th graded piece of char.
    
    pol = 0
    for k in dictM:
        pol += dictM[k]*q^k
    
    return pol