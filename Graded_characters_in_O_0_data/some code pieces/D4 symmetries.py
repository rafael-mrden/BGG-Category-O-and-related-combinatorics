def symm(w,u,v):
    '''Symmetry corresponding to replacing two nodes u and v in the Dynkin diagram.'''
    
    if w == e:
        return e
    
    u = str(u)
    v = str(v)
    if u[0] != "s":
        u = "s"+u
    if v[0] != "s":
        v = "s"+v
        
    return eval(str(w).replace(u,"0").replace(v,u).replace("0",v))

def symm1(w):
    return symm(w,1,3)

def symm2(w):
    return symm(w,3,4)

def orbit(w):
    return [w, symm1(w), symm2(w), symm1(symm2(w)), symm2(symm1(w)), symm1(symm2(symm1(w)))]
