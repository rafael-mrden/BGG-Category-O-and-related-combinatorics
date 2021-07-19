'''This code searches for elementary isomorphisms and elementary decompositions of theta_a(L_b).
See my notebook from 2020-05-12.'''


def summands(s,x):
    '''Needed: s simple and x<s*x.
    Returns the summands of theta_x theta_s as a list.'''
    
    result = [s*x]  
    for z in W.bruhat_interval(e,x):
        if (s*z).bruhat_le(z):
            result += [z]*mu(z,x)
    
    return result


def middle(s,y):
    '''Needed: s simple and y*s<y.
    Returns the middle of theta_s L_y as a list.'''
    
    result = [y*s] 
    for z in W.bruhat_interval(y,w0):
        if z.bruhat_le(z*s):
            result += [z]*mu(y,z)
    
    return result


def check(x,y,s):
    
    if (s*x).bruhat_le(x) or y.bruhat_le(y*s):
        return []

    if R_smaller(x,y.inverse()):
        return []
    
    surviving_summands = []    
    for z in summands(s,x):

        if R_smaller(z,y.inverse()):
            surviving_summands.append(z)
    
    surviving_middle = []
    for z in middle(s,y):

        if R_smaller(x,z.inverse()):
            surviving_middle.append(z)
    
    if len(surviving_summands) >= 1 and len(surviving_middle) >= 1:        
        equation = [set(),set()]

        for z in surviving_summands: 
            equation[0].add( (z,y) )
        
        for z in surviving_middle:  
            equation[1].add( (x,z) )
        
        return equation
    
    return []

    
def print_equation(equation):
    string1 = ""
    for p in equation[0]:
        string1 += " + M(%s,%s)"%( convert_to_123(p[0]) , convert_to_123(p[1]) )        

    string2 = ""
    for p in equation[1]:
        string2 += " + M(%s,%s)"%( convert_to_123(p[0]) , convert_to_123(p[1]) )
        
    string3 = string1[3:] + " == " + string2[3:]
    print(string3)
    
    

    
    
############################################
        
elem = []  

for x in W:
    for y in W:
        for s in W.simple_reflections():
            eq = check(x,y,s)
            
            if eq != [] and eq not in elem and list(reversed(eq)) not in elem:
                elem.append( eq )
                print_equation(eq)
            

for eq in elem:
    KReq = [
        set([ (pair[1].inverse()*w0, w0*pair[0].inverse()) for pair in eq[0]]),
        set([ (pair[1].inverse()*w0, w0*pair[0].inverse()) for pair in eq[1]]),
    ]
    if KReq not in elem and list(reversed(KReq)) not in elem:
        elem.append(KReq)
        print_equation(KReq)
 
 

 
elem12 = []    # Equations A+B=C or A=B+C

for eq in elem:
    if (len(eq[0]) == 1 and len(eq[1]) == 2) or (len(eq[0]) == 2 and len(eq[1]) == 1):
        elem12.append(eq)

elem22 = []    # Equations A+B=C+D

for eq in elem:
    if (len(eq[0]) == 2 and len(eq[1]) == 2):
        elem22.append(eq)

for eq in elem:
    if (len(eq[0]) > 2 or len(eq[1]) > 2):
        print("We have a big relation:")
        print_equation(eq) 
 
     
############################################ 
############################################
 
 
 
G_dict = {}

for x in W:
    for y in W:
        if x != e and y != w0:
            if (x,y.inverse()) in R:
                G_dict[(x,y)] = []
				
				

for eq in elem:
    if (len(eq[0]) == 1 and len(eq[1]) == 1):
        
        p = list(eq[0])[0]
        q = list(eq[1])[0]
        
        if p in G_dict:
            G_dict[p].append(q)
        else:
            G_dict[p]=[q]

        if q in G_dict:
            G_dict[q].append(p)
        else:
            G_dict[q]=[p]

            
            
G_components = Graph(G_dict).connected_components()
print(len(elem), len(G_components))



for comp in G_components:
    
    print ("{ " + " , ".join(["%s * %s" %(convert_to_123(p[1]),convert_to_123(p[0]))   for p in  comp ]) + " }:")
    
    print("")
    
    X = M(comp[0][0],comp[0][1])
    X.name = ""
    
    print_clean(X)




########################################################

decomposables = []
for eq in elem12:
    lef,rig = eq
    if len(lef)==1:
        dec = list(lef)[0]
    if len(rig)==1:
        dec = list(rig)[0]
    if dec not in decomposables:
        decomposables.append(dec)
		
		
f=open("???.txt", "r")    # Change file
Equal_Products =f.read()
f.close()




Equal_Products_lines = Equal_Products.split("\n")




for pair in decomposables:
    x,y = pair
    prod = " %s * %s "%(convert_to_123(y),convert_to_123(x))
    
   
    for line in Equal_Products_lines:
        if prod in line:
            i = Equal_Products_lines.index(line)
            
            if "decomposable" not in line:
                line = line[:-1] + " (decomposable):"
                Equal_Products_lines[i] = line




string = ""
for line in Equal_Products_lines:
    string += line + "\n"
	
	
print(string)