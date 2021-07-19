# Checks wether theta_x L_y is a sum of two other such characters in the Grothendieck group.
# Note: does not check for more than two summands, and also theta_x L_y coud be equal to a sum of objects which are not of the form theta_x' L_y'



thetas = []
for aa in W:
    for b in W:
        if (aa,b.inverse()) in R:
            X = M(aa,b)
            if X not in thetas:
                thetas.append(X)

def is_decomposable(X):
    for aa in thetas:
        for b in thetas:
            if aa+b==X:
                return True
    return False