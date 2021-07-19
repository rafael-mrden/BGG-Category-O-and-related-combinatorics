#######################################################################

# Projective functors combinatorics
# Here we have "Th", independent of the previous "theta", "M", "th".


def mult_one(m):
    '''Used in printing of the ProjFunctor class.'''
    
    if m == 1:
        return ""
    else:
        return str(m)+"*"


class ProjFunctor:
    '''A class representing a graded projective functor on the graded version of category O.
    self.content a dictionary {(w,i):m for w in W, i and m integers} representing m*theta_w<i>.'''
    
    def __init__(self):
        self.content = {}
           
    def __str__(self): # Temporary
        
        grades = {}
        for p in self.content:
            if p[1] not in grades.keys():
                grades[p[1]] = [p]
            else:
                grades[p[1]].append(p)
        
        string = ""
        for i in reversed(sorted(grades.keys())):
            string += " + ".join( ["%s%s<%d>"%(mult_one(self.content[p]), convert_to_123(p[0]), p[1]  )   for p in grades[i]  ]   ) + "\n"
        
        return string
            

    def cleaning(self):
        '''Removes summands with multiplicity zero from the dictionary.'''

        remove = []
        for p in self.content.keys():
            if self.content[p] == 0:
                remove.append(p)
        for p in remove:
            self.content.pop(p)
        
        return self
        

#    def __eq__(self, other): TODO


    def __add__(self,other):    
        summa = ProjFunctor()
        summa.content = copy(self.content)
        
        for p in other.content.keys():
            
            if p in summa.content.keys():
                summa.content[p] += other.content[p]
            
            else:
                summa.content[p] = other.content[p]
        
        return summa.cleaning()
   
    
    def __mul__(self,other):
        prod = ProjFunctor()
        prod.content = copy(self.content)
        
        for p in prod.content.keys():
            prod.content[p] = other*(prod.content[p])
        
        return prod.cleaning()
    
    
    def __rmul__(self,other):
        return self*other
    
    
    def __sub__(self,other):
        return self+((-1)*other)
    
    
    def __neg__(self):
        return (-1)*self
    
    
def sh(fun,i):
    '''Shifts the projective functor fun by i.'''

    res = ProjFunctor()
    
    for p in fun.content:
        res.content[p[0], p[1]+i] = fun.content[p]

    return res
    
    
def Th(w,i):
    '''This returns theta_w<i>'''
    
    result = ProjFunctor()
    result.content = {(w,i):1}
    return result


def multiply(xs,w):
    '''This returns the decomposition of theta_xs * theta_w.'''
    
    factors = xs.reduced_word()
    
    if len(factors) == 0:  # xs == e
        return Th(w,0)
    
    s = W.simple_reflections()[factors[len(factors)-1]]
    
    if len(factors) == 1:  # xs == s
        
        if (w*s).bruhat_le(w):
            return Th(w,1) + Th(w,-1)
        
        else:
            result = Th(w*s,0)
            for y in W.bruhat_interval(e,w):
                if (y*s).bruhat_le(y):
                    result = result + mu(y,w)*Th(y,0)
            return result
    
    x = xs * s    
    step = multiply(x,w)    
    subtract = [ y for y in W.bruhat_interval(e,x) if (y*s).bruhat_le(y)]
    
    res = ProjFunctor()
    
    for k in step.content.keys():
        res += step.content[k] * sh(  multiply(s,k[0]) , k[1])
    
    for y in subtract:
        res -= mu(y,x) * multiply(y,w)
    
    return res
    


def action(P,X):
    '''Returns the graded character of P(X), where P is a projective functor and X a graded character.'''
    
    result = char_0()
    for p in P.content:
        result = result + P.content[p] * shift( theta(p[0],X), p[1])
    
    return result