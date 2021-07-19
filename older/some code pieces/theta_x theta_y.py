def organize_multiply(d,y):
    product = multiply(d,y)
    
    part1 = ProjFunctor()    # theta_z<a> where z sim_H y
    part2 = ProjFunctor()    # theta_z<a> where z sim_J y
    part3 = ProjFunctor()    # theta_z<a> where a = a(d)
    part4 = ProjFunctor()    # theta_z<a> where a = a(x)
    remainder = ProjFunctor()
    
    if CartanType(W)[0] == 'A':
        H = {y}
    else:
        H = set(L_cell(y)).intersection(set(R_cell(y)))
    
    
    J = two_cell(y)
        
    ad = a(d)
    ay = a(y)
    
    
    
    for t in product.content:
        if t[0] in H:
            part1.content[t] = product.content[t]
        else:
            if t[0] in J:
                part2.content[t] = product.content[t]
            else:
                if t[1] == ad:
                    part3.content[t] = product.content[t]
                else:
                    if t[1] == ay:
                        part4.content[t] = product.content[t]
                    else:
                        remainder.content[t] = product.content[t]
    
    print("H_y={%s}:"%(",".join([convert_to_123(w) for w in list(H)])))                
    print(part1)
    
    print("J_y:")
    print(part2)
    
    print("a(d)=%d:"%ad)
    print(part3)
    
    print("a(y)=%d:"%ay)
    print(part4)
    print(remainder)

            
    
    
    