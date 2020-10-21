penultimate = two_cell( s1*w0  )   # Works only in type A


def all_reduced(w):
    if w not in W:
        w = convert_from_123(w)
    if w == e:
        return []
    return [eval("".join([str(d) for d in word])) for word in w.reduced_words()]


# 3D pictures:

from sage.plot.plot3d.shapes import Text
from sage.plot.plot3d.shapes2 import frame3d


ptsize = 25
ptsize_small = 10

def coords(w): # For w in the penultimate cell
    if w not in W:
        w = convert_from_123(w)
    return ( eval(convert_to_123(list(AL(w))[0])), eval(convert_to_123(list(AR(w))[0])),  w0.length()-w.length())


def picture_coshadow_Delta(t):
    
    tetra = point3d([])
    
    appears = {} # Here we collect those w from the penultimate that will appear in the picture, and in which heights.
    
    for w in penultimate:
        
        p_ew = KLP(e,w)
        p_tw = KLP(t,w)
        
        appears[w] = [m.degree() for m in (p_ew-p_tw).monomials()]
        
        x,y,z = coords(w)

        for i in [m.degree() for m in p_ew.monomials()]:
            if i in appears[w]:
                tetra += point3d((x,y, z + 2*i  ), size=ptsize)
                
            else:
                tetra += point3d((x,y, z + 2*i ), size=ptsize_small,color='red')

            tetra += Text(convert_to_123(w)).translate(x,y, z + 2*i -0.15)
            tetra += Text(  (coords(w)[0],coords(w)[1],coords(w)[2]+2*i)   ).translate(x,y, z + 2*i +0.15)
    
    remove_keys_with_value(appears,[])
    
    socle = deepcopy(appears)
    
    for u in appears:
        for v in appears:
            if mu(u,v) != 0:
                
                ux,uy,uz = coords(u)
                vx,vy,vz = coords(v)
                
                for i in appears[u]:
                    for skip in [1,-1]:
                        j = (uz - vz - skip)/2 + i
                        if j in appears[v]:
                            tetra += line3d([(ux,uy,uz + 2*i),(vx,vy,vz + 2*j)])
                            if skip==1 and i in socle[u]:
                                socle[u].remove(i)
                            if skip==-1 and j in socle[v]:
                                socle[v].remove(j)
     
    remove_keys_with_value(socle,[])  
    return (tetra,socle)


def tetra(w):
    if w not in W:
        w = convert_from_123(w)
    
    print( ", ".join( str(t) for t in all_reduced(w))   )
    picture = picture_coshadow_Delta(w)[0]
    socle = picture_coshadow_Delta(w)[1]
    
    print("Socle:")
    
    socle_heights = reversed(sorted(list(set([coords(w)[2]+2*socle[w][0] for w in socle ]))))
    for h in socle_heights:
        print(", ".join( sorted([str((coords(w)[0], coords(w)[1], coords(w)[2]+2*socle[w][0])) for w in socle if  coords(w)[2]+2*socle[w][0]==h]) ) )
    
    show(picture,frame=False)

    

#  Shadow and coshadow as graded characters:

def shadow(X):
    Y = char_0()
    Y.name=""
    for i in components(X):
        Y.component[i] = {}
        for w in X.component[i]:
            
            if w in penultimate:
                Y.component[i][w] = X.component[i][w]
    Y.remove_zeros()
    return Y

def shadow_Delta(x):
    return shadow(shift(char_Delta(x),-x.length()))

def coshadow_Delta(x):
    return shadow( char_Delta(e) - shift(char_Delta(x),-x.length()) )


# Matrices:

def matrix_coshadow_Delta(x):
    picture = [[0 for i in range(n)] for j in range(n)]
    
    for w in penultimate:
        i,j = [eval(convert_to_123(list(AL(w))[0])),eval(convert_to_123(list(AR(w))[0]))]
        
        pol = KLP(e,w)-KLP(x,w)
        if pol != 0:
            picture[i-1][j-1] = pol
    
    return matrix(picture)


# Graphs:

def picture_coshadow(x):
    picture = [["0" for i in range(n)] for j in range(n)]
    
    for w in penultimate:
        i,j = [eval(convert_to_123(list(AL(w))[0])),eval(convert_to_123(list(AR(w))[0]))]
        
        pol = KLP(e,w)-KLP(x,w)
        if pol != 0:
            picture[i-1][j-1] = str(pol).replace(" ","")
    
    return picture


#P = W.bruhat_poset().dual()
#Nodes = dict( (w, (convert_to_123(w) +"\n"+ str(table(picture_coshadow(w), align = "center",frame=True)))) for w in W)
#P.plot(element_labels=Nodes, figsize = 70, element_shape='s', element_size=14000)



shift = (n+1)/2

def coords_shifted(w): # For w in the penultimate cell
    if w not in W:
        w = convert_from_123(w)
    return ( -shift + eval(convert_to_123(list(AL(w))[0])), -shift + eval(convert_to_123(list(AR(w))[0])), -1+ w0.length()-w.length())

def beaut_picture(t, *args):
    
    tetra = point3d([])
    
    appears = {} # Here we collect those w from the penultimate that will appear in the picture, and in which heights.
    
    for w in penultimate:
        
        p_ew = KLP(e,w)
        p_tw = KLP(t,w)
        
        appears[w] = [m.degree() for m in (p_ew-p_tw).monomials()]
        
        x,y,z = coords_shifted(w)

        for i in [m.degree() for m in p_ew.monomials()]:
            if i in appears[w]:
                tetra += point3d((x,y, z + 2*i  ), size=ptsize, color='blue')
                
            else:
                tetra += point3d((x,y, z + 2*i ), size=ptsize_small, color='red')

#            tetra += Text(convert_to_123(w)).translate(x,y, z + 2*i -0.10)
            if args[0]=="labels":
                tetra += Text(  (shift+ coords_shifted(w)[0],shift+ coords_shifted(w)[1],w0.length()-coords_shifted(w)[2]-2*i)   ).translate(x,y, z + 2*i +0.10)
    
    remove_keys_with_value(appears,[])
    
    for u in appears:
        for v in appears:
            if mu(u,v) != 0:
                
                ux,uy,uz = coords_shifted(u)
                vx,vy,vz = coords_shifted(v)
                
                for i in appears[u]:
                    for skip in [1,-1]:
                        j = (uz - vz - skip)/2 + i
                        if j in appears[v]:
                            tetra += line3d([(ux,uy,uz + 2*i),(vx,vy,vz + 2*j)], color = "blue", thickness = arrow_thickness)

    tetra +=  frame3d([1-shift,1-shift,0],vector([n-shift,n-shift,n-1]), opacity=.5, color="black",linestyle="dashed")
    
    return tetra
	
	
	
	
path = "Graded_characters_in_O_0_data/graphics/"+ CartanType(W)[0]+str(CartanType(W)[1])+"/"
ptsize = 15
arrow_thickness = 3

CT = CartanType(W)[0] +str(CartanType(W)[1])
path = "Graded_characters_in_O_0_data/graphics/" #+ CT +"/"

if not os.path.isdir(path): 
    os.mkdir(path)  
            
for w in [w0]:
    td = beaut_picture(w,"")
    
#    folder = path+convert_to_123(w)
#    if not os.path.isdir(folder): 
#        os.mkdir(folder)
    
    for i in [2,20]: # range(40)
        td.rotate((0,0,1),i*pi/40).save_image(path+"/"+CT+"_angle%d.png"%i, frame=False, figsize=16)	
	
	
	
	
############################### Checking if posets of bigrassmannians with the same L/R descends are totally ordered.
W_poset = W.bruhat_poset()
bigrassmannians = [x for x in W if len(DR(x))==1 and len(DL(x))==1]


for i in range(n):
    for j in range(n):
        bigrassmannians_fixed = [x for x in bigrassmannians if DR(x)=={W.simple_reflections()[i]} and DL(x)=={W.simple_reflections()[j]}]
        if bigrassmannians_fixed != []:
            P = W_poset.subposet(bigrassmannians_fixed)
            if not P.is_chain():
                show(P.plot())
                print("*************************")
				
# Some variations:
for i in range(n):
    for j in range(n):
        bigrassmannians_fixed = [x for x in bigrassmannians if DR(x)=={W.simple_reflections()[i]} and DL(x)=={W.simple_reflections()[j]}]
        if bigrassmannians_fixed != []:
            P = W_poset.subposet(bigrassmannians_fixed)
            if not P.is_chain():
                
                Q = P.relabel(lambda x: convert_to_123(x))
                colors =  { "blue" : [ Q(convert_to_123(x)) for x in bigrassmannians_fixed if x not in join_irreducibles]}
                
                show(Q.plot(figsize = 5,cover_color="gray",element_colors = colors,element_shape = "s" )) #
                
                print("*************************")

				


# Some combinatorics:

N = n+1

def to_perm(w): # Different from the one in the main file.
    '''Converts an element from W to a permutation of type list.
    Mind that it is inverted in the proces.'''
    
    check_if_type_A()
    
    if w==1: w=W(1)
        
    return list(w.to_permutation())


def to_W(t): # Different from the one in the main file.
    '''The inverse of "to_perm".'''
    
    check_if_type_A()
    
    t = Permutation(list(t))
    red = prod(W.simple_reflections()[i] for i in t.reduced_word())
    
    if red==1: red=W(1)
    return red


def eval_perm(w,i):
    return to_perm(w)[i-1]


def Ess(w):
    if w not in W:
        w = convert_from_123(w)
    ess = set()
    for i in range(1,N):
        for j in range(1,N):
            if i < eval_perm(w.inverse(),j) and j < eval_perm(w,i) and eval_perm(w,i+1)<= j and eval_perm(w.inverse(),j+1)<=i:
                ess.add((i,j))
    return ess

def diagram(w):
    if w not in W:
        w = convert_from_123(w)
    M = [[0 for i in range(N)] for j in range(N)]
    for i in range(N):
        M[i][eval_perm(w,i+1)-1] = 1
    return matrix(M)
    

def b(i,j):
    if i <= j:
        return W.from_reduced_word(range(i,j+1))
    if i>j:
        return W.from_reduced_word(range(j,i+1))
        

def soc(w):
    dict_soc = picture_coshadow_Delta(w)[1]
    return { coords(w)[:2] for w in dict_soc.keys()  }


#def rank_fun(w):
#    return matrix( [[len([q for q in range(1,i+1) if eval_perm(w,q)<=j ]) for j in range(1,N+1)] for i in range(1,N+1)  ] )

def rank_fun(w,i,j):
    return len([q for q in range(1,i+1) if eval_perm(w,q)<=j ])
	
    
def corank_fun(w):
    return matrix( [[min(i,j)-len([q for q in range(1,i+1) if eval_perm(w,q)<=j ]) for j in range(1,N+1)] for i in range(1,N+1)  ] )


def Ess_RWY(w):
    not_below = [x for x in W if x not in W.bruhat_interval(e,w)]
    return (W_poset.subposet(not_below)).minimal_elements()


def soc2(w):
    if w not in W:
        w = convert_from_123(w)
    socle_ungr = [tuple(reversed(list(t))) for t in Ess(w)]
    socle = []
    for t in socle_ungr:
        i,j = t
        

        k = min(i-1,j-1,N-i-1,N-j-1) - min(i,j) + rank_fun(w,j,i)
        
        for wji in penultimate:
            if coords(wji)[:2]==(j,i):
                break
        
        socle.append( (i,j,wji.length()-2*(k+1))  )
    
    return socle


BI = [] 
for x in W:
    for y in Ess_RWY(x):
        y = eval(str(y))
        if y not in BI:
            BI.append(y)




def b(i,j,k):
    string = ""
    if i<=j:
        for s in reversed(range(i-k,i+1)):
            for t in range(s,s+k+j-i+1):
                string += str(t)
    else:
        for s in range(j-k,j+1):
            for t in reversed(range(s,s+k+i-j+1)):
                string += str(t)

    return convert_from_123(string)

def b2(i,j,k):
    r = min(i-1,j-1)-k  
    return to_W(list(range(1,r+1)) + list(range(i+1,i+j-r+1)) + list(range(r+1,i+1)) + list(range(i+j-r+1,n+2)))
    
