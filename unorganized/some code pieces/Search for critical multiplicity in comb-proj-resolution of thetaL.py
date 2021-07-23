D = Duflo_Involutions()

for x in D:
    ax = earliest_occurence(x)
    for d in D:
#        y = w0*d
        ad = earliest_occurence(d)
        X = M(d,w0*x)
        if X != char_0():
            res = proj_resolution_param(X,30)
            critical = res[ax].get((w0,ad-ax),0)
            aw0d = a(w0*d)
            bxw0d = b(x,w0*d)
            if critical != 1 0 or (critical == 1 and aw0d + bxw0d != len(res)-1):
                print("x=%s, d=%s, mult=%d, pd=%d, c.pd=%d" %(convert_to_123(x),convert_to_123(d),critical,aw0d + bxw0d,len(res)-1))
        
print("finished")