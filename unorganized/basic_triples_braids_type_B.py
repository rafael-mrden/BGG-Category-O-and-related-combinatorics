def b(k,p,q):
    
    if p<=0 or q==0 or k<=0 or k<=1-p-q or (p==1 and q<0):
        raise ValueError("Not a basic triple!")
    
    if q>=p:
        result = tuple( list(range(1,p)) + list(range(-q-k+1,-q+1)) + list(range(p,q)) )
    
    elif p>q and q>0:
        result = tuple( list(range(1,q)) + list(range(q+k,p+k)) + list(range(-q-k+1,-q+1)) )
    
    elif k>-q and -q>0:
        result = tuple( list(range(k+1,p+k)) + list(range(-k,q)) + list(range(1,-q+1)) )
    
    elif -q>=k:
        result = tuple( list(range(1,-q-k+1)) + list(range(-q+1,p+k)) + list(range(-q-k+1,-q+1)) )
    
    else:
        raise ValueError("Not a basic triple!")        

    return result


def notation(i):
    if i >= 0:
        return str(i)
    if i < 0:
        return "\\overline{%s}"% str(-i)


def braid(arg, n):
    
    if n < len(arg):
        raise ValueError("Size too small!")
    
    height = 2
    xfactor = 0.35

    print("\\scalebox{.9}{")
    print("\\begin{tikzpicture}[scale=1, baseline=-3mm]")
    
    print("\\draw[red,fill=white] (0,0) -- (0,%f);"%height)
    print("\\node[below] at (0,0) {\\tiny $0$};")
#    print("\\draw[red,fill=white] (0,0) circle (1pt) node[below] {};")
#    print("\\draw[red,fill=white] (0,%f) circle (1pt) node[below] {};" %height)
    
    for i in range(1,n+1):
  
        if i<=len(arg):
            value = arg[i-1]
        else:
            value = i

        print("\\draw (%f,0) -- (%f,%f);" %(i*xfactor, value*xfactor, height))
        print("\\draw (%f,0) -- (%f,%f);" %(-i*xfactor, -(value*xfactor), height))

        print("\\node[below] at (%f,0) {\\tiny $%s$};"%(i*xfactor,notation(i)))
        print("\\node[below=-0.45mm] at (%f,0) {\\tiny $%s$};"%(-i*xfactor,notation(-i)))
        
#        print("\\draw[fill=white] (%f,0) circle (1pt) node[below] {};"%(i*xfactor))
#        print("\\draw[fill=white] (%f,%f) circle (1pt) node[below] {};"% (i*xfactor, height ))
#        print("\\draw[fill=white] (%f,0) circle (1pt) node[below] {};"%(-i*xfactor))
#        print("\\draw[fill=white] (%f,%f) circle (1pt) node[below] {};"%(-i*xfactor, height))
        
    print("\\end{tikzpicture}}")    

    
def all_basic_triples(i,j,n):
    result = []
    
    p = j+1
    
    q = i+1
    for k in range(1+max(0,1-p-q), n+2-max(p,q)):
        result.append( (k,p,q) )
    
    if p>1 and i!=0:
        q = -i
        for k in range(1+max(0,1-p-q), n+2-max(p,q)):
            result.append( (k,p,q) )
    
    return result


def length_b(triple):
    k,p,q = triple
    
    if q>0:
        return (p+q-1)*k + k*(k-1)/2
    
    elif k> -q and -q>0:
        return p*k + k*(k-1)/2 - (-q+1)*(-q)/2
    
    elif -q >= k:
        return (p+q+k-1)*k
    
    else:
        raise ValueError("Not a basic triple!")
    


def braid_modified(arg, n):
    '''Has dots on every 2(mod3).'''
    
    if n < len(arg):
        raise ValueError("Size too small!")
    
    height = 3
    xfactor = 0.6

    print("\\scalebox{1.0}{")
    print("\\medmuskip=-2mu")
    print("\\begin{tikzpicture}[scale=1, baseline=-3mm]")
    
    print("\\draw[red,fill=white] (0,0) -- (0,%f);"%height)
    print("\\node[below] at (0,0) {\\tiny $0$};")
    print("\\node[above] at (0,%f) {\\tiny $0$};"%height)
#    print("\\draw[red,fill=white] (0,0) circle (1pt) node[below] {};")
#    print("\\draw[red,fill=white] (0,%f) circle (1pt) node[below] {};" %height)
    
    for i in range(1,n+1):
  
        if i<=len(arg):
            value = arg[i-1]
        else:
            value = i
        
        if i%3 != 2:
            print("\\draw (%f,0) -- (%f,%f);" %(i*xfactor, value*xfactor, height))
            print("\\draw (%f,0) -- (%f,%f);" %(-i*xfactor, -(value*xfactor), height))

            print("\\node[below] at (%f,0) {\\tiny $%s$};"%(i*xfactor,notation(i)))
            print("\\node[below=-0.45mm] at (%f,0) {\\tiny $%s$};"%(-i*xfactor,notation(-i)))

            print("\\node[above] at (%f,%f) {\\tiny $%s$};"%(i*xfactor, height, notation(i)))
            print("\\node[above] at (%f,%f) {\\tiny $%s$};"%(-i*xfactor, height, notation(-i)))
        
#        print("\\draw[fill=white] (%f,0) circle (1pt) node[below] {};"%(i*xfactor))
#        print("\\draw[fill=white] (%f,%f) circle (1pt) node[below] {};"% (i*xfactor, height ))
#        print("\\draw[fill=white] (%f,0) circle (1pt) node[below] {};"%(-i*xfactor))
#        print("\\draw[fill=white] (%f,%f) circle (1pt) node[below] {};"%(-i*xfactor, height))
        
        else:
            print("\\node[below] at (%f,0) {\\tiny $\\cdots$};"%(i*xfactor))
            print("\\node[below] at (%f,0) {\\tiny $\\cdots$};"%(-i*xfactor))

            print("\\node[above] at (%f,%f) {\\tiny $\\cdots$};"%(i*xfactor, height))
            print("\\node[above] at (%f,%f) {\\tiny $\\cdots$};"%(-i*xfactor, height))            
        
    print("\\end{tikzpicture}}")    

    

 
    