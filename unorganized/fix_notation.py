def fix_notation(l):
    '''This function converts a reduced word from GAP3 notation, to the corresponding word in the SageMath notation.
    To be run in SageMath, with initial setting e.g.:
    W = WeylGroup("E7", prefix="s")
    [s1,s2,s3,s4,s5,s6,s7] = W.simple_reflections()'''
    
    if CartanType(W)[0] in ['A', 'E', 'F', 'G']:
        return W.from_reduced_word(l)
    
    if CartanType(W)[0] in ['B', 'C', 'D', 'H']:
        new_l = []
        for d in l:
            new_l.append(n+1-d)
        return W.from_reduced_word(new_l)
    
       
    raise ValueError("Conversion in this type not yet implemented.")