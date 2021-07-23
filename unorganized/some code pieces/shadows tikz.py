N = n+1

penultimate = two_cell(s1*w0)

def coords(w): # For w in the penultimate cell
    if w not in W:
        w = convert_from_123(w)
    return ( eval(convert_to_123(list(AL(w))[0])), eval(convert_to_123(list(AR(w))[0])))


def transform(p):
    x,y,z = p
    return tuple((x,y,-z))


def label(p):
#    return "\\tiny $%s$"%(p,)   # Comment this line if labels are not wanted
    return ""


def tetrahedron_tikz():
    points = []
    for i in range(1,N):
        for j in range(1,N):
            for k in range(0, 1+min(i-1,j-1,N-1-i,N-1-j) ):
                points.append( (i,j,N*(N-1)/2 - abs(N-i-j)-1-2*k))
    
    lines = []
    for p in points:
        for q in points:
            if ((p[0] == q[0] and abs(p[1]-q[1]) == 1) or (p[1] == q[1] and abs(p[0]-q[0]) == 1)) and p[2] == q[2]-1:
                lines.append((p,q))
    
    print("%% Tetrahedron, n=%d"%N)
    print("\\tdplotsetmaincoords{110}{130} \n\\begin{tikzpicture}[tdplot_main_coords, scale=1.3]\n")  

    # Print points
    for p in points:
        print("\\filldraw[black] %s circle (1.5pt) node[anchor=west] {%s};"%(transform(p),label(p)))
        
    print()

    # Print lines
    for line in lines:
        p,q = line

        print("\\draw %s -- %s;"%(transform(p),transform(q)))
    
    print()
    
    # Print cube
    top  = N*(N-1)/2 - N + 1
    bottom = N*(N-1)/2 - 1
    A = [(1,1,border) for border in [bottom,top]]
    B = [(1,N-1,border) for border in [bottom,top]]
    C = [(N-1,N-1,border) for border in [bottom,top]]
    D = [(N-1,1,border) for border in [bottom,top]]
    
    for i in range(2):
        print("\\draw[gray, dashed] %s -- %s -- %s -- %s -- cycle;" %( transform(A[i]), transform(B[i]), transform(C[i]), transform(D[i])))
    for p in [A,B,C,D]:
        print("\\draw[gray, dashed] %s -- %s;" %(transform(p[0]),transform(p[1])) )
    print("\\draw[gray, dashed] %s -- %s;" %(transform(B[0]),transform(D[0])) )
    print("\\draw[gray, dashed] %s -- %s;" %(transform(A[1]),transform(C[1])) )    
                
    print("\n\\end{tikzpicture}\n")        
#    return (points,lines)




def socle_tikz(x):
    if x not in W:
        x = convert_from_123(x)
        
    def label(p):
    #    return "\\tiny $%s$"%(p,)   # Comment this line if labels are not wanted
        return ""

    points = [] 
    shadow = []
    
    for w in penultimate:
        i,j = coords(w)
        
        p_ew = KLP(e,w)
        p_xw = KLP(x,w)
        
        in_shadow = [m.degree() for m in (p_ew-p_xw).monomials()]
        
        for k in range(0, 1+min(i-1,j-1,N-1-i,N-1-j) ):
            points.append( (i,j,N*(N-1)/2 - abs(N-i-j)-1-2*k) )
            if k in in_shadow:
                shadow.append( (i,j,N*(N-1)/2 - abs(N-i-j)-1-2*k) )
               
    lines = []
    for p in points:
        for q in points:
            if ((p[0] == q[0] and abs(p[1]-q[1]) == 1) or (p[1] == q[1] and abs(p[0]-q[0]) == 1)) and p[2] == q[2]-1:
                lines.append((p,q))

    socle = []
    for p in shadow:
        good = True
        for line in lines:
            if line[0] == p and line[1] in shadow:
                good = False
                break
        if good:
            socle.append(p)

    
    print("%% coshadow(%s), n=%d"%(convert_to_123(x),N))
    print("\\tdplotsetmaincoords{110}{130} \n\\begin{tikzpicture}[tdplot_main_coords, scale=1.3]\n")  

    # Print cube
    top  = N*(N-1)/2 - N + 1
    bottom = N*(N-1)/2 - 1
    A = [(1,1,border) for border in [bottom,top]]
    B = [(1,N-1,border) for border in [bottom,top]]
    C = [(N-1,N-1,border) for border in [bottom,top]]
    D = [(N-1,1,border) for border in [bottom,top]]
    
    for i in range(2):
        print("\\draw[gray, dashed] %s -- %s -- %s -- %s -- cycle;" %( transform(A[i]), transform(B[i]), transform(C[i]), transform(D[i])))
    for p in [A,B,C,D]:
        print("\\draw[gray, dashed] %s -- %s;" %(transform(p[0]),transform(p[1])) )
    print("\\draw[gray, dashed] %s -- %s;" %(transform(B[0]),transform(D[0])) )
    print("\\draw[gray, dashed] %s -- %s;" %(transform(A[1]),transform(C[1])) )  
    
    print()

    # Print lines
    for line in lines:
        p,q = line
        if p in shadow and q in shadow:
            print("\\draw %s -- %s;"%(transform(p),transform(q)))
        else:
            print("\\draw[gray, dashed] %s -- %s;"%(transform(p),transform(q)))
    
    print()
    
    # Print points
    for p in points:
        if p in socle:
            print("\\draw[fill=red] %s circle (2.0pt) node[anchor=west] {%s};"%(transform(p),label(p)))
        else:
            if p in shadow:
                print("\\filldraw[black] %s circle (1.5pt) node[anchor=west] {%s};"%(transform(p),label(p)))
            else:
                print("\\draw[fill=white] %s circle (1.5pt) node[anchor=west] {%s};"%(transform(p),label(p)))
                
    print("\n\\end{tikzpicture}")        
#    return (socle,shadow,points,lines)











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

def Di(w):
    if w not in W:
        w = convert_from_123(w)
    ess = set()
    for i in range(1,N):
        for j in range(1,N):
            if i < eval_perm(w.inverse(),j) and j < eval_perm(w,i):
                ess.add((i,j))
    return ess

def diagram(w):
    if w not in W:
        w = convert_from_123(w)
    M = [[0 for i in range(N)] for j in range(N)]
    for i in range(N):
        M[i][eval_perm(w,i+1)-1] = 1
    return matrix(M)
    

    

def soc(w):
    dict_soc = picture_coshadow_Delta(w)[1]
    return { coords(w)[:2] for w in dict_soc.keys()  }


def rank_fun(w):
    return matrix( [[len([q for q in range(1,i+1) if eval_perm(w,q)<=j ]) for j in range(1,N+1)] for i in range(1,N+1)  ] )

#def rank_fun(w,i,j):
#    return len([q for q in range(1,i+1) if eval_perm(w,q)<=j ])

def corank_fun(w):
    return matrix( [[min(i,j) - len([q for q in range(1,i+1) if eval_perm(w,q)<=j ]) for j in range(1,N+1)] for i in range(1,N+1)  ] )


def latex_Ess(w):
    if w not in W:
        w = convert_from_123(w)
        
    print("%% Ess(%s)"%convert_to_123(w))
        
    table = diagram(w)
    D = Di(w)
    ess = Ess(w)
    
    print("\\begin{tikzpicture}[scale=0.8, yscale=-1]")
    
    print("\n% Axes:")
    print("\\draw (0.5,-0.25)--(0.5,%s);"%str(N+0.25)[:4])
    print("\\draw (-0.25,0.5)--(%s,0.5);"%str(N+0.25)[:4])
    for i in range(1,N+1):
        print("\\node[] at (0,%s) {\\small{$%s$}};" %(i,i))
        print("\\node[] at (%s,0) {\\small{$%s$}};" %(i,i))

    print("\n%  Points:")
    for i in range(1,N+1):
        for j in range(1,N+1):

            if eval_perm(w,i) == j:
                print("\\draw (%d,%d)--(%d,%d);"%(j,i,j,N))                
                print("\\draw (%d,%d)--(%d,%d);"%(j,i,N,i))          
                print("\\draw[fill=white] (%d,%d) circle (2pt);"%(j,i))
                
            if (i,j) in D:
                if (i,j) in ess:
                    print("\\draw (%d,%d) node[cross] {};"%(j,i))
                    print("\\node at (%d,%d) [rectangle,draw] {};"%(j,i)) 
                else:
                    print("\\draw (%d,%d) node[cross] {};"%(j,i))  
 
    print("\\end{tikzpicture}")
