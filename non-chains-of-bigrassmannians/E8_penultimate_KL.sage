from sage.interfaces.gap import  get_gap_memory_pool_size, set_gap_memory_pool_size
set_gap_memory_pool_size(24364842180)



W = WeylGroup("E8", prefix="s")
[s1,s2,s3,s4,s5,s6,s7,s8] = W.simple_reflections()

##################################################################################

n = rank(W)
w0 = W.long_element()
e = W(1)


####### Kazhdan-Lusztig polynomials ##########################

# A faster implementation of KL-polynomials (using the optional package Coxeter 3) is given by this
# Fokko Ducloux’s Coxeter3 C++ library.

# Had to install it: I just typed "sage -i coxeter3" in the terminal.

# It seems that one can direcly coerce from WeylGroup to CoxeterGroup and vice versa.
# I will therefore use CoxeterGroup to calculate KL-polynomials, but for all other Bruhat business I will use WeylGroup.

##################################################################################

n = rank(W)
w0 = W.long_element()
e = W(1)


####### Kazhdan-Lusztig polynomials ##########################

# A faster implementation of KL-polynomials (using the optional package Coxeter 3) is given by this
# Fokko Ducloux’s Coxeter3 C++ library.

# Had to install it: I just typed "sage -i coxeter3" in the terminal.

# It seems that one can direcly coerce from WeylGroup to CoxeterGroup and vice versa.
# I will therefore use CoxeterGroup to calculate KL-polynomials, but for all other Bruhat business I will use WeylGroup.

R.<q> = LaurentPolynomialRing(QQ)

KL = KazhdanLusztigPolynomial(W,q)  # KL-polynomials implemented in standard Sage way
# http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/kazhdan_lusztig.html


CoxeterPackage = CoxeterGroup(W, implementation="coxeter3")

def KLP(x,y):
    '''Returns the KL-polynomial, implemented in "Coxeter3" package by Fokko du Cloux.
    http://math.univ-lyon1.fr/~ducloux/coxeter/coxeter3/english/coxeter3_e.html'''
    
    if x not in W:
        x = convert_from_123(x)
    if y not in W:
        y = convert_from_123(y)
    
    return CoxeterPackage.kazhdan_lusztig_polynomial(CoxeterPackage(x), CoxeterPackage(y))
    # If "coxeter3" is not installed, remove the line 'CoxeterPackage = CoxeterGroup(W, implementation="coxeter3")'
    # and in this function return KL.P(x,y)

#Point:
#    - standard Sage way: KL.P(x,y)
#    - faster way: KLP(x,y) 


def mu(w,x):
    '''Returns the KL mu-function with arguments w,x.
    By Humphrey's BGG book p. 175 and p. 169, for w<x we have:
    mu(x,w) = mu(w0*w,w0*x),
    mu(w,x) = dim Exit^1 (L_w,L_x) = dim Exit^1(L_x,L_w) = dim Exit^1(Delta_x,L_w).'''

    if w not in W:
        w = convert_from_123(w)
    if x not in W:
        x = convert_from_123(x)
        
    if w.bruhat_le(x):
        poly_dict = KLP(w,x).dict()       
        j = (x.length()-w.length()-1)/2 
        if j not in poly_dict.keys():
            return 0
        return poly_dict[j]

    return 0


def convert_to_123(w):
    '''Converts an element from W to the "123" string notation.
    Does not work with coefficients, as "convert_to_123_long".'''
    
    if w == W(1):
        return "e"
    
    return str(w).replace("s","").replace("*","")


def convert_from_123(string):
    '''Converts one element from W in the "123" string notation to the usual "s1*s2*s3" notation.'''
    
    if type(string)== Integer:
        string = str(string)
        
    if string == "e":
        return W(1)
    
    string = "*".join([char for char in string])
    
    for i in range(1,n+1):
        string = string.replace(str(i),"s%s"%i)
    
    return eval(string)


def DR(w):
    '''Returns the set of simple right descents of w.'''
    
    return {W.simple_reflections()[i] for i in w.descents()}


def DL(w):
    '''Returns the set of simple left   descents of w.'''
        
    return DR(w.inverse())


def AL(w):
    '''Returns the set of simple left ascends of w.'''
    
    DescLe = list(DL(w))
    AscLe = [s for s in W.simple_reflections() if s not in DescLe]
    return set(AscLe)


def AR(w):
    '''Returns the set of simple right ascends of w.'''
    
    DescRi = list(DR(w))
    AscRi = [s for s in W.simple_reflections() if s not in DescRi]
    return set(AscRi)



####################################

def coeff_one(i):
    if i != 1:
        return str(i)
    else:
        return ""
    
def KLPv_e(w):
    p = KLP(e,w)
    l = w.length()
    pv = sum(p.coefficients()[i] * q^(l - 2*p.exponents()[i]) for i in range(len(p.exponents())))
    
    # Format pv nicely for latex:
    
    return " + ".join([ coeff_one(pv.coefficients()[i])+"v^{%s}"%pv.exponents()[i] for i in range(len(pv.exponents())) ])




# This is only for E8

small = [ convert_from_123(str(x)) for x in [1, 13, 134, 1342, 1345, 13456, 134567, 1345678,
 2, 24, 243, 245, 2456, 24567, 245678,
 3, 34, 345, 3456, 34567, 345678,
 4, 45, 456, 4567, 45678,
 5, 56, 567, 5678,
 6, 67, 678,
 7, 78,
 8]] + [ convert_from_123(str(x)).inverse() for x in [13, 134, 1342, 1345, 13456, 134567, 1345678,
 24, 243, 245, 2456, 24567, 245678,
 34, 345, 3456, 34567, 345678,
 45, 456, 4567, 45678,
 56, 567, 5678,
 67, 678,
 78]]

penultimate = [ w0*x for x in small ]






penultimate_table = [ [e for i in range(n)] for i in range(n) ]

for x in penultimate:
    i = eval(convert_to_123(list(AL(x))[0]))
    j = eval(convert_to_123(list(AR(x))[0]))
    penultimate_table[i-1][j-1]=x



def penultimate_two_cell_E_latex():
    print("\\begin{tabular}{|l||%s|}"%("|").join(["l"]*n) + " \\hline")
    
    print("$%s_{%s}$ & "%(CartanType(W)[0],str(CartanType(W)[1])) + " & ".join([ "$%d$"%j for j in range(1,n+1) ]) + " \\"+"\\"+" \\hline\\hline")
    
    for i in range(n):
        
        print( "$%d$ & "%(i+1) + " & ".join( [ "$"+ str(KLPv_e(penultimate_table[i][j])) +"$" for j in range(n)   ] ) + " \\"+"\\"+" \\hline" )

    print("\\end{tabular}\n")



# penultimate_two_cell_E_latex()



w = penultimate[0]
print("w = %s"%convert_to_123(w))

l = w.length()
print("length = %d"%l)

print("\n****\n")

for i in range(1,l):
    x = convert_from_123((convert_to_123(w)[:i]))
    print("length = %d"%i)
    print("x = %s"%convert_to_123(x))
    print("KLP(e,x) = %s\n"%KLP(e,x) )
	
    











