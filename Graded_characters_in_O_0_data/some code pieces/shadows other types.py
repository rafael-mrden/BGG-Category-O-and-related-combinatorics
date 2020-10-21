W = WeylGroup("D4", prefix="s")
[s1,s2,s3,s4] = W.simple_reflections()

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




############### BIGRASSMANNIANS ###############

bigrassmannians = [x for x in W if len(DR(x))==1 and len(DL(x))==1]
print("#bigrassmannians = %d"%len(bigrassmannians))
print("bigrassmannians = %s"%bigrassmannians)
print()


############### JOIN-IRREDUCIBLES ###############

i=0
join_irreducibles = []
bigrassmannians_copy = [x for x in bigrassmannians]

for y in W:
    inter = W.bruhat_interval(e,y)
    complement = [z for z in W if z not in inter]
    
    remove_b = []
    
    for b in bigrassmannians_copy: # Is b minimum in complement?
        if b in complement:
            is_min = 1
            for z in complement:
                if z.bruhat_le(b) and z != b:
                    is_min = 0
                    break
            if is_min == 1:
                print("%s is join-irreducible!"%b)
                join_irreducibles.append(b)
                if b not in remove_b:
                    remove_b.append(b)
                    
    for b in remove_b:
        bigrassmannians_copy.remove(b)
  
    i += 1
    print("%d/%d"%(i,len(W)))
    

print("***************************************************\n")
print("#join_irreducibles = %d"%len(join_irreducibles))
print("join_irreducibles = %s"%join_irreducibles)
print()

############### DISSECTORS ###############

def is_dissector(x):
    A = W.bruhat_interval(x,w0)
    B = [y for y in W if y not in A]
    for y in B:
        if set(B) == set(W.bruhat_interval(e,y)):
            return True
    return False

i=0
dissectors = []
for x in bigrassmannians:
    if is_dissector(x):
        print("%s is a dissector!"%x)
        dissectors.append(x)
    i += 1
    print("%d/%d"%(i,len(bigrassmannians)))

print("***************************************************\n")

print("#dissectors = %d"%len(dissectors))
print("dissectors = %s"%dissectors)
print()



############### PENULTIMATE ###############

penultimate = two_cell(s1*w0)

# Option 2, only for type D (see my notes 13.7.2020. for type B).
def path(i,j):
    if j<i:
        return path(j,i)[::-1] # after assumed i<j 
    if (i,j)==(n-1,n):
        return "%d%d%d"%(n-1,n-2,n)  
    res = ""
    for k in range(i,j+1):
        res += str(k)
    return res.replace("%d%d"%(n-1,n),str(n))

small_cell = []
for i in range(1,1+n):
    for j in range(1,n+1):
        small_cell.append(convert_from_123(path(i,j)))
penultimate = [x*w0 for x in small_cell]

# This prints the sum of coeff. of all KL polynomails P_e,w for w in penultimate.
pol = 0
for w in penultimate:
    pol += KLP(e,w)
print(pol(q=1))



############### PLOTTING NON-CHAINS OF BIGRASSMANNIANS ###############
W_poset = W.bruhat_poset()

default_vertex_color = "#fec7b8"
for i in range(1,n+1):
    for j in range(1,n+1):
        bigrassmannians_fixed = [x for x in bigrassmannians if DR(x)=={W.simple_reflections()[i]} and DL(x)=={W.simple_reflections()[j]}]
        if bigrassmannians_fixed != []:
            P = W_poset.subposet(bigrassmannians_fixed)
            if not P.is_chain():
                Q = P.relabel(lambda x: convert_to_123(x))
                colors =  { "yellow" : [ Q(convert_to_123(x)) for x in bigrassmannians_fixed if x not in dissectors],
                           "blue" : [ Q(convert_to_123(x)) for x in bigrassmannians_fixed if x not in join_irreducibles],
                          default_vertex_color : [Q(convert_to_123(x)) for x in bigrassmannians_fixed if x in dissectors]}
                
                grap = Q.plot(figsize = 8.3,cover_color="gray",element_colors = colors,element_shape = "s" )
                show(grap)
 #               grap.save('plots/%d%d.png'%(i,j))      
                print("*************************")
                
                



# Later version

badpairs = []
B = join_irreducibles
for i in range(1,n+1):
    for j in range(1,n+1):
        B_ij = [x for x in join_irreducibles if DL(x)=={W.simple_reflections()[i]} and DR(x)=={W.simple_reflections()[j]}]
        
        for x in B_ij:
            for y in B_ij:
                inter = W.bruhat_interval(x,y)
                if len([z for z in inter if z in B_ij])==2:
                    
                    exists = 0
                    
                    for w in W:
                        if ( W.simple_reflections()[i]  not in DL(w) or W.simple_reflections()[j] not in DR(w)) and x.bruhat_le(w) and w.bruhat_le(y):
                            exists = 1
                            break
                    
#                    if exists == 1:
#                        print(convert_to_123(x),convert_to_123(y),"ok",convert_to_123(w))
                    if exists == 0:
                        badpairs.append((convert_to_123(x)+"(%d)"%x.length(),convert_to_123(y)+"(%d)"%y.length()))
                        print("%s-%s"%(convert_to_123(x),convert_to_123(y)))
             
            
default_vertex_color = "#fec7b8"
for i in range(1,n+1):
    for j in range(1,n+1):
        join_irreducibles_fixed = [x for x in join_irreducibles if DR(x)=={W.simple_reflections()[i]} and DL(x)=={W.simple_reflections()[j]}]
        if join_irreducibles_fixed != []:
            P = W_poset.subposet(join_irreducibles_fixed)
            if True: #not P.is_chain():
                Q = P.relabel(lambda x: convert_to_123(x)+"(%d)"%x.length())
                #colors =  { #"yellow" : [ Q(convert_to_123(x)) for x in bigrassmannians_fixed if x not in dissectors],
                           #"blue" : [ Q(convert_to_123(x)) for x in bigrassmannians_fixed if x not in join_irreducibles],
                          #default_vertex_color : [Q(convert_to_123(x)) for x in bigrassmannians_fixed if x in dissectors]}
                
                cov_lb = []
                for cov in Q.cover_relations():
                    a,b = cov
                    if (str(a),str(b)) in badpairs:
                        cov_lb.append(cov+["N"])
                    else:
                        cov_lb.append(cov+[""])
                
                cov_red = []
                cov_gray = []
                for cov in Q.cover_relations():
                    a,b = cov
                    if (str(a),str(b)) in badpairs or (str(b),str(a)) in badpairs:
                        cov_red.append(cov)
                    else:
                        cov_gray.append(cov)
                
                grap = Q.plot(figsize = 8.3,cover_color="gray",element_shape = "s", cover_labels = cov_lb,
                              cover_colors = {"red" : cov_red, "gray" : cov_gray}) #element_colors = colors, 
                show(grap)
                
                if not os.path.isdir("posets"):  
                    os.mkdir("posets") 
                path = "posets/" + CartanType(W)[0]+str(CartanType(W)[1])
                if not os.path.isdir(path): 
                    os.mkdir(path)   
                    
                grap.save(path+'/J_%d%d.png'%(i,j))      
                print("*************************")                   
                    
        
              
############### JOINS ###############              
                
W_poset = W.bruhat_poset()
def join(S):
    SS = [convert_from_123(a) for a in S if a not in W] + [a for a in S if a in W]
    
    U = set(W.bruhat_interval(SS[0],w0))
    for a in SS[1:]:
        U = U.intersection(set(W.bruhat_interval(a,w0)))
        
    minU = (W_poset.subposet(list(U))).minimal_elements()
    
    if len(minU)==1:
        j = minU[0]
        return eval(convert_to_123(str(j)))
    else:
        return minU


def JM(w):
    if w not in W:
        w = convert_from_123(w)
    result = W_poset.subposet([x for x in join_irreducibles if x.bruhat_le(w) ]).maximal_elements()
    return [convert_from_123(convert_to_123(x)) for x in result]

def BM(w):
    if w not in W:
        w = convert_from_123(w)
    return W_poset.subposet([x for x in bigrassmannians if x.bruhat_le(w) ]).maximal_elements()

def JM2(w):
    if w not in W:
        w = convert_from_123(w)
    DLw = DL(w)
    DRw = DR(w)
    result = []
    for x in JM(w):
        x = convert_from_123(convert_to_123(x))
        if list(DL(x))[0] in DLw and list(DR(x))[0] in DRw:
            result.append(x)
    return result
