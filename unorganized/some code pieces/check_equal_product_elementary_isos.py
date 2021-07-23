'''This code checks whether all equalities of graded characters of thetas come from composition of elementary isomorphisms.
See my notebook from 2020-05-12.'''


all_paths = cells_graph("left").all_simple_paths(trivial=True)
R = R_relation()

EP = Equal_Products()
EP1 = [set(t) for t in EP[0] if len(t)>1]


i = 0


relation3 = set()

for x in W:
    for y in W:
        
        print("%d/%d"%(i,len(W)^2))
        i += 1
        
        if M(x,y) == char_0():
            
            for s in W.simple_reflections():
                if (y*s).bruhat_le(y) and (x).bruhat_le(s*x):

                    condition1 = 1
                    for z in W.bruhat_interval(e,x):
                        if (s*z).bruhat_le(z):
                            if mu(z,x) != 0:
                                if M(z,y) != char_0():
                                    condition1 = 0

                    if condition1 == 1:
                                             
                        middle = list((M(s,y).component[0].keys()))
                        x_middle = [z for z in middle if M(x,z) != char_0() ]
                        
                        if len(x_middle) == 1:
                            z = x_middle[0]
                            relation3.add( ( (s*x,y ),(x,z) ) )
                        
                        if len(x_middle) > 1:
                            print("Alert")
                            print(x,y,s,x_middle)
                        

                        
                        
vertices = [ p for p in cartesian_product([W,W]) if M(p[0], p[1]) != char_0() ]                       


G_dict = {}

for v in vertices:
    G_dict[v] = [p[1] for p in relation3 if p[0] != p[1] and p[0]==v]+[p[0] for p in relation3 if p[0] != p[1] and p[1]==v]


EP2 = []

for g in Graph(G_dict).connected_components_subgraphs():
    if len(g.vertices())>1:
        EP2.append(set(g.vertices()))
        
        
for x in EP1:
    if not x in EP2:
        print(x)

print("*********")
for x in EP2:
    if not x in EP1:
        print(x)