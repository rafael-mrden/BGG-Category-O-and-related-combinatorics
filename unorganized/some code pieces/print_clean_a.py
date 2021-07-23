def print_clean_a(char):
    '''Prints cleaner output: simple reflections are denoted by 1,2,3, ...
    Prints also the value of a function in brackets next to each composition factor.'''
    
    char_dummy = GradedChar()
    
    char_dummy.name = (char.name).replace("(1)","(e)").replace("s","").replace("*","")
    
    for k in char.component:
        char_dummy.component[k] = {}
        
        for w in char.component[k]:
            char_dummy.component[k][convert_to_123(w)+"("+str(a(w))+")"] = char.component[k][w]
    
    print(char_dummy)   

    
def print_clean_a_equals_grade(char):
    '''Prints cleaner output: simple reflections are denoted by 1,2,3, ...
    Prints also the value of a function in brackets next to each composition factor.'''
    
    char_dummy = GradedChar()
    
    char_dummy.name = (char.name).replace("(1)","(e)").replace("s","").replace("*","")
    
    for k in char.component:
        char_dummy.component[k] = {}
        
        for w in char.component[k]:
            aw = a(w)
            if aw == k:
                char_dummy.component[k][convert_to_123(w)] = char.component[k][w]
    
    print(char_dummy) 