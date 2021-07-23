def br_le(x,y):
    return x.bruhat_le(y)  


def W_subposet(L):
    '''This is equivalent to W.bruhat_poset().subposet(),
    but without calculating the full W.bruhat_poset().'''
      
    return Poset((L,br_le))


def join_old(S):
    SS = [convert_from_123(a) for a in S if a not in W] + [a for a in S if a in W]
    
    U = set(W.bruhat_interval(SS[0],w0))
    for a in SS[1:]:
        U = U.intersection(set(W.bruhat_interval(a,w0)))
        
    minU = (W_subposet(list(U))).minimal_elements()
    
    return minU


def join(S):
    S = [convert_from_123(a) for a in S if a not in W] + [a for a in S if a in W]   

    ld = set()
    rd = set()
    for s in S:
        ld = ld.union(DL(s))
        rd = rd.union(DR(s))

    U = []
    for y in [x for x in W if DL(x).issubset(ld) and  DR(x).issubset(rd)]:
        is_upper = True
        for s in S:
            if not s.bruhat_le(y):
                is_upper = False
        if is_upper == True:
            U.append(y)
    
    return (W_subposet(U)).minimal_elements()      
    
    
def JM_(w):
    if w not in W:
        w = convert_from_123(w)
    result = W_subposet([x for x in join_irreducibles if x.bruhat_le(w) ]).maximal_elements()
    return [convert_from_123(convert_to_123(x)) for x in result]

def BM(w):
    if w not in W:
        w = convert_from_123(w)
    return W_subposet([x for x in bigrassmannians if x.bruhat_le(w) ]).maximal_elements()

def JM__(w):
    if w not in W:
        w = convert_from_123(w)
    DLw = DL(w)
    DRw = DR(w)
    result = []
    for x in JM_(w):
        x = convert_from_123(convert_to_123(x))
        if list(DL(x))[0] in DLw and list(DR(x))[0] in DRw:
            result.append(x)
    return result

def JM(w):
    if w not in W:
        w = convert_from_123(w)
    DLw = DL(w)
    DRw = DR(w)
    result = W_subposet([x for x in join_irreducibles if x.bruhat_le(w) and (list(DL(x))[0] in DLw) and (list(DR(x))[0] in DRw) ]).maximal_elements()
    return [convert_from_123(convert_to_123(x)) for x in result]   



BG = {}
for x in bigrassmannians:
    i = eval(convert_to_123(list(DL(x))[0]))
    j = eval(convert_to_123(list(DR(x))[0]))
    if (i,j) not in BG:
        BG[(i,j)] = []  
    BG[(i,j)].append(x)
print("BG done!")

JI = {}
for x in join_irreducibles:
    i = eval(convert_to_123(list(DL(x))[0]))
    j = eval(convert_to_123(list(DR(x))[0]))
    if (i,j) not in JI:
        JI[(i,j)] = []  
    JI[(i,j)].append(x)
print("JI done!")

#load("JM_data/{}/JM_{}".format(CartanType(W)[0]+str(CartanType(W)[1]), convert_to_123(w0)))