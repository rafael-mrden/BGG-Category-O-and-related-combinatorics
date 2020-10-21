all_paths = cells_graph("left").all_simple_paths(trivial=True)
L = L_relation()
R = R_relation() # a quicker way is just to invert L.



glpd = 2*w0.length()+1

m = len(W)*len(W[1:])
i = 0
problem = 0


for x in W[1:]:
    for y in W:
        
        if (x,y.inverse()) in R:
            X = M(x,y)
            PR = proj_resolution_param(X,glpd)
            
            if len(PR)-1 == a(w0*x) + b(y.inverse()*w0, w0*x.inverse()):
                
                P0 = PR[0]
                for w in P0:
                    if not ( (w[0],y) in R and (y,w[0]) in R):
                        print ("P0 problem with:")
                        print(x,y)
                        print(P0)
                        problem += 1
                
                for P in PR[:-1]:
                    for w in P:
                        if not ( (x,w[0]) in L ):
                            print ("Problem with:")
                            print(x,y)
                            print(PR)
                            problem += 1
                            
                Plast = PR[-1]
                for w in Plast:
                    if not ( (w[0],x) in L and (x,w[0]) in L):
                        print ("P_end problem with:")
                        print(x,y)
                        print(Plast)
                        problem += 1
        i+=1
        print("%d/%d"%(i,m))
        
print("Finished with %d problems!"%problem)
            