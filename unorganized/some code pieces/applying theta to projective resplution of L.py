x = s2*s3*s1*s2*s3*s2
y = s3*s2*s3 

i = 0
PR = proj_resolution_param_L(x)


for P in PR:
    print("P_%d:"%i)
    for p in PR[P]:
        print( "%s%s<%d>"%(mult_one(PR[P][p]),convert_to_123(p),-i)  )
    print()
    i+=1

    
###################################### 
    

theta_PR = {}

for i in PR:
    theta_PR[i] = {}
    
    for p in PR[i]:
        
        m = PR[i][p]
        
        prod = multiply(y,p).content
        
        for q in prod:
            qq = (q[0], q[1]-i)
            
            if qq not in theta_PR[i].keys():
                theta_PR[i][qq] = m*prod[q]
            else:
                theta_PR[i][qq] += m*prod[q]
    
    
######################################    
    
C = char_0()
for i in theta_PR:
    for p in theta_PR[i]:
        C += (-1)^i * shift(char_P(p[0]),p[1]) * theta_PR[i][p]

C == M(y,x)


###################################### 


for i in sorted(list(theta_PR.keys())):
    print("P_%s:"%i)
    for p in theta_PR[i]:
        print("%s%s<%d>"%(mult_one(theta_PR[i][p]),convert_to_123(p[0]),p[1]))
    print()
    
    
theta_PR_cleaned = deepcopy(theta_PR)


###################################### 


for i in reversed(range(1,13)):
    for p in set(theta_PR_cleaned[i].keys()).intersection(set(theta_PR_cleaned[i-1].keys())):
        a = theta_PR_cleaned[i][p]
        b = theta_PR_cleaned[i-1][p]
        m = min(a,b)
        theta_PR_cleaned[i][p] = a-m
        theta_PR_cleaned[i-1][p] = b-m


for i in range(13):
    remove_keys_with_value(theta_PR_cleaned[i]  , 0)

remove_keys_with_value(theta_PR_cleaned  , {})