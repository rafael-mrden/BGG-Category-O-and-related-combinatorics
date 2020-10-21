######################################################################
# Combinatorial projective resolutions of for general graded characters


def proj_resolution_param(X,t):
    '''Returns the parameters of the combinatorial projective resolution of a graded character X.
    Here t is a bound on the length of the resolution, to avoid to long computations.
    Output is a list, and the elements are
    dictionaries with keys parameters (w,shift) and values are multiplicities.
    Algorithm is given by Hankyung.
    TODO: Make t optional variable.'''
   
    resolution_param = []

    while X != char_0():

        P_param = {}   # to be the set of parameters in the projective resolution
        
        X1 = char_0()  # to be the next X
        C = copy(X)

        while C != char_0():
            
            top_P = char_0()
            top_deg = min(components(C))
            for w in C.component[top_deg]:
                
                top_P += C.component[top_deg][w] * shift(char_P(w),-top_deg)   
                
                if (w,-top_deg) in P_param:
                    P_param[(w,-top_deg)] += C.component[top_deg][w]
                else:
                    P_param[(w,-top_deg)] = C.component[top_deg][w]
                

            K = combinatorial_kernel(top_P,C)    
            C = combinatorial_cokernel(top_P,C)

            X1 += K
        
        resolution_param.append(P_param)
        
        if len(resolution_param)>t+1: # Here exit if too big.
            return [False]*1000
        
        X = X1
    
    return resolution_param


def proj_resolution(X,t):
    '''Returns the combinatorial projective resolution of a graded character X, as a list of graded characters.
    Here t is a bound on the length of the resolution, to avoid to long computations.
    TODO: Make t optional variable.'''
    
    resolution = []
    resolution_param = proj_resolution_param(X,t)

    for i in range(len(resolution_param)):
        P_i = char_0()
        for p in resolution_param[i]:
            P_i += shift(char_P(p[0]),p[1])*resolution_param[i][p]
            
        resolution.append(P_i)
    
    return resolution


def proj_dim(X):
    '''Returns the combinatorial projective dimension of X
    (assuming it exists and is less than 1000.'''
    
    return len(proj_resolution_param(X,1000))-1
    

    







# Sorting W by length and lexicographically and printing comb. proj. dim. of theta(Delta)'s and theta(Nabla)'s.

ListW_lex = [convert_to_123(w) for w in W if w != e]
ListW_lex.sort()
ListW_lex = ['e']+ListW_lex
def sort_lex(w):
    return ListW_lex.index(convert_to_123(w))
ListW_len_lex = [convert_to_123(w) for w in W if w != e]
ListW_len_lex.sort()
ListW_len_lex.sort(key=len)
ListW_len_lex = ['e']+ListW_len_lex
def sort_len_lex(w):
    return ListW_len_lex.index(convert_to_123(w))
def l(w):
    return w.length()
W_list = [convert_from_123(x) for x in ListW_len_lex]



print("\\begin{tabular}{|l||%s|}"%("|").join(["l"]*(len(W))) + " \\hline")
print("$y \\setminus x$ & " + " & ".join(["$"+convert_to_123(x)+"$" for x in W_list]) + " \\"+"\\"+" \\hline \\hline" )

for y in W_list:
    print( "$%s$ & "%convert_to_123(y) + " & ".join(["$"+str(proj_dim(theta(y,char_Delta(x))))+"$" for x in W_list]) + " \\"+"\\"+" \\hline" )

print("\\end{tabular}")

print("\n\n\n")

print("\\begin{tabular}{|l||%s|}"%("|").join(["l"]*(len(W))) + " \\hline")
print("$y \\setminus x$ & " + " & ".join(["$"+convert_to_123(x)+"$" for x in W_list]) + " \\"+"\\"+" \\hline \\hline" )

for y in W_list:
    print( "$%s$ & "%convert_to_123(y) + " & ".join(["$"+str(proj_dim(theta(y,char_Nabla(x))))+"$" for x in W_list]) + " \\"+"\\"+" \\hline" )

print("\\end{tabular}")