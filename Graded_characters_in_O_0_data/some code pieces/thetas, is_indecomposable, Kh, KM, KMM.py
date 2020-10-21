thetas = []
for aa in W:
    for b in W:
        if (aa,b.inverse()) in R:
            X = M(aa,b)
            if X not in thetas:
                thetas.append(X)

def is_decomposable(X):
	
    for b in thetas:
        if X-b in thetas:
                return True
    return False
	
	
	

	
def is_decomposable_M(x,y):
    
    X = M(x,y)
    thetas = []
    
    for u in L_cell(x):
        for v in R_cell(y):
            if R_smaller(u,v.inverse()):
                Y = M(u,v)
                if Y not in thetas:
                    thetas.append(Y)
    
    for i in range(len(thetas)):
        X1 = thetas[i]
        if X1.is_true_character():
            if X-X1 in thetas:
                return True
        print("%d/%d"%(i,len(thetas)))
    return False
	
	
	

def Kh(y):
    for z in W:
        if (z,y.inverse()) in R:
            X = M(z,y)
            for t in W:
                if t != z and (t,y.inverse()) in R:
                    Y = M(t,y)
                    if X == Y:
                        return False
    return True
	
	
	
	
def KMM(d,x):
    if not R_smaller(d,x.inverse()):
        return 0
    
    X = M(d,x)
    m = mult(X,a(d),x)
    
    return(m)
		
		
########################################################



file = open("Graded_characters_in_O_0_data/some data/B4 elementary isomorphisms ala Johan.txt", "r")
file_string = file.read()
file.close()
equal_products = [ line for line in file_string.split("\n") if line != "" and line[0]=="{"]

def check_notKM(y):
    if y in W:
        y = convert_to_123(y)
    y = str(y)
    for ep in equal_products:
        if " "+y+" *" in ep:
            if "decomposable" in ep:
                print(y + ": " +str(ep))
                return True
    return y + ": I do not know"
	
def check_notKh(y):
    if y in W:
        y = convert_to_123(y)
    y = str(y)
    for ep in equal_products:
        if ep.count(" "+y+" *")>1:
#            print(y + " : " + ep)
            return True
    return y + ": I do not know"

	
def is_simple_extreme(X):

    extreme = max(components(X))
    extreme_part = X.component[extreme]
    
    if len(extreme_part.keys())>1:
        return False
    
    if extreme_part[list(extreme_part.keys())[0]] != 1:
        return False
    
    return True	
	

# Lemma 8.28:

y = 434234123431
y = convert_from_123(y)
for x in W:
    if R_smaller(x,y.inverse()):
        d = [q for q in R_cell(x) if q in D][0]
        X = M(x.inverse(),x)
        m = KMM(d,y)
        
        if not(m==1 and a(x)==max(components(X)) and is_simple_extreme(X)):
            print("x="+convert_to_123(x))
            print("d="+convert_to_123(d))
            print("KMM(d,y)=%d"%m)
            print("a(x)=%d"%a(x))
            print_clean(X)
print("finished")


# Lemma 8.29:

y = 4341234231
y = convert_from_123(y)
print("y="+convert_to_123(y))
print("a(y)=%d"%a(y))
print_clean(M(y.inverse(),y))
d = [q for q in R_cell(w0*y.inverse()) if q in D][0]
print("d="+convert_to_123(d))
print("a(d)=%d"%a(d))
print("\n****\n")
for x in W:
    if R_smaller(x,y.inverse()):
        z = x.inverse()*w0
        if z != w0:
            
            m = KMM(d,z)
            if m != 1:
                print("z="+convert_to_123(z))
                print("KMM(d,z)=%d"%m)
print("finished")








################################################################################

def decompositions(X,summands):
    '''Tries to write "X" as a sum of elements from "summands".'''
    
    if X == char_0():
        return [[]]
    
    results = []
    
    for Y in summands:
        X_ = X-Y
        if X_.is_true_character():
            for smaller_sum in decompositions(X_,summands):
                X_sum = [Y] + smaller_sum
                results.append ( [Y] + smaller_sum )
    
    return results

def decompositions_theta_L(x,y):
    '''Tries to decompose theta_x(L_y) as a sum of other thetaL's.'''
    
    if type(x) == Integer or type(x) == str:
        x = convert_from_123(x)
    if type(y) == Integer or type(y) == str:
        y = convert_from_123(y)
    
    if x == e or y == w0 or not R_smaller(x,y.inverse()):  # Exclude trivial simple and tilting modules.
        return []
    
    X = M(x,y) 
    summands = []
    
    for x_ in L_cell(x):
        for y_ in R_cell(y):
            if R_smaller(x_,y_.inverse()):
                X_ = M(x_,y_)
                if X_ not in summands:
                    summands.append(X_)
                    
#    print("mark")    
            
    return decompositions(X,summands)



i=0

for x in W:
    if x !=e:
        for y in W:
            if y != w0:
                if R_smaller(x,y.inverse()):
                    
                    decs = decompositions_theta_L(x,y)
                    
                    for dec in decs:
                        if len(dec)>1:
                            if len(dec)>2:
                                print("ALERT!")
                            print("theta_%s(L(%s)) = "%(convert_to_123(x),convert_to_123(y)) + " + ".join([convert_to_123_long(s.name) for s in dec  ]))
                            
            i+=1
            print(i)