def Equal_Products():
    '''Returns a pair of two lists L1 and L2 of the same length.
    For each i, L1[i] is the set of pairs (x,y) such that theta_x(L(y)) == L2[i].
    L2 consists of all mutually distinct non-trivial charaters that can be obtained as theta_x(L(y)),
    excluding the simple modules ( theta_e(L(y)) ) and the tilting modules ( theta_x(L(w0)) ).
    
    Works for A4 and B3, but not above.'''
    
    table = ([],[])

    for x in W:
        for y in W:

            if x != e and y != w0:

                if R_smaller(x,y.inverse()):   # Maybe want to edit here?, e.g., "if (x,y.inverse()) in R"?
                    
                    X = M(x,y)
                    X.name =""

                    if X in table[1]:

                        ind = (table[1]).index(X)
                        table[0][ind].append( (x,y) )

                    else:

                        (table[1]).append(X)
                        (table[0]).append([(x,y)])
    return table  
	
	
	
'''This code prints nicely equal products.'''


EP = Equal_Products()

for i in range(len(EP[1])):
    
    equal_factors = EP[0][i]
    product = EP[1][i]
    
    print ("{ " + " , ".join(["%s * %s" %(convert_to_123(p[1]),convert_to_123(p[0]))   for p in  equal_factors ]) + " }:")
    
    print("")
    
    print_clean(product)