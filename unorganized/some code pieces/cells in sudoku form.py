D = Duflo_Involutions()

ListW_lex = [convert_to_123(w) for w in W if w != e]
ListW_lex.sort()
ListW_lex = ['e']+ListW_lex
def sort_lex(w):
    return ListW_lex.index(convert_to_123(w))

ListW_len_lex = [convert_to_123(w) for w in W if w != e]
ListW_len_lex.sort()
ListW_len_lex.sort(key=len)
ListW_len_lex = ['e']+ListW_len_lex
def sort_len_lex(w):
    return ListW_len_lex.index(convert_to_123(w))




def l(w):
    return w.length()



def list_descends(S):
    result = [convert_to_123(s) for s in S]
    result.sort()
    return ", ".join(result)


def two_cells_html():
    
    print("<b>KL-cells:</b>\n<br /><br />\n")
    
    W_dummy = list(W)
    
    W_dummy.sort(key=sort_len_lex)

    while W_dummy != []:
              
        Cell = two_cell(W_dummy[0])
        DufCell = [d for d in Cell if d in D]
        DufCell.sort(key=sort_lex)

        DufCell_dict = {}
        for d in DufCell:
            Ld = L_cell(d)
            Rd = [x.inverse() for x in Ld]
            DufCell_dict[d]=(Rd,Ld)

        Cell_table = [[[] for x in range(len(DufCell))] for x in range(len(DufCell))]
        for i in range(len(DufCell)):
            for j in range(len(DufCell)):

                intersect = [x for x in DufCell_dict[DufCell[i]][0] if x in DufCell_dict[DufCell[j]][1] ]
                intersect.sort(key=l)
                if i == j:
                    intersect.remove(DufCell[i])
                    intersect = [DufCell[i]] + intersect
    
                Cell_table[j][i] = "<br />".join([convert_to_123(x) for x in intersect])
        
        
        print('<table border="1">')
        print("<tr>")
        print("<td>a=%d</td>"%a(DufCell[0]))
        
        for i in range(len(DufCell)):
            print("<td>%s</td>"%list_descends(DL(DufCell[i])))
        
        print("</tr>")
        
        
        for i in range(len(Cell_table)):
            
            row = Cell_table[i]
            print("<tr>")
            print("<td nowrap>%s</td>"%list_descends(DR(DufCell[i])))
            for cel in row:
                print("<td>%s</td>"%cel)
            
            print("</tr>")
        print("</table>")
        
        
        print("<br /><br />\n")


        W_dummy = [x for x in W_dummy if x not in Cell]


    print("\n")




def two_cells_latex():
    
    W_dummy = list(W)
    
    W_dummy.sort(key=sort_len_lex)

    while W_dummy != []:
              
        Cell = two_cell(W_dummy[0])
        DufCell = [d for d in Cell if d in D]
        DufCell.sort(key=sort_lex)

        DufCell_dict = {}
        for d in DufCell:
            Ld = L_cell(d)
            Rd = [x.inverse() for x in Ld]
            DufCell_dict[d]=(Rd,Ld)

        Cell_table = [[[] for x in DufCell] for x in DufCell]
        for i in range(len(DufCell)):
            for j in range(len(DufCell)):

                intersect = [x for x in DufCell_dict[DufCell[i]][0] if x in DufCell_dict[DufCell[j]][1] ]
                intersect.sort(key=l)
                if i == j:
                    intersect.remove(DufCell[i])
                    intersect = [DufCell[i]] + intersect
    
                if len(intersect)==1:
                    Cell_table[j][i] = "$"+convert_to_123(intersect[0])+"$"
                else:
                    Cell_table[j][i] = "\\begin{tabular}{@{}l@{}}" + (" \\"+"\\ ").join(["$"+convert_to_123(x)+"$" for x in intersect]) + "\\end{tabular}"
                    
                    
        
        print("\\begin{tabular}{|%s|}"%("|").join(["l"]*(len(DufCell)+1) ))
        
        print("\\hline")
        print("a=%d &"%a(DufCell[0]) + " & ".join([ " %s "%list_descends(DL(DufCell[i])) for i in range(len(DufCell))] ) )
        print(" \\"+"\\")
        
        print("\\hline")
        for row in Cell_table:

            print(" & ".join([list_descends(DR(DufCell[i]))] + row) + " \\"+"\\"+" \\hline")
     
        print("\\end{tabular}\n")
        
        print("\\hskip 1cm \n")


        W_dummy = [x for x in W_dummy if x not in Cell]




def penultimate_two_cell_latex():
	
	Cell = two_cell(s1*w0)
	DufCell = [d for d in Cell if d in D]
	DufCell.sort(key=sort_lex)

# manual sorting:	DufCell = [ DufCell[2],  DufCell[3], DufCell[4], DufCell[1], DufCell[0], ]

	DufCell_dict = {}
	for d in DufCell:
		Ld = L_cell(d)
		Rd = [x.inverse() for x in Ld]
		DufCell_dict[d]=(Rd,Ld)

	Cell_table = [[[] for x in DufCell] for x in DufCell]
	for i in range(len(DufCell)):
		for j in range(len(DufCell)):

			intersect = [x for x in DufCell_dict[DufCell[i]][0] if x in DufCell_dict[DufCell[j]][1] ]
			intersect.sort(key=l)
			if i == j:
				intersect.remove(DufCell[i])
				intersect = [DufCell[i]] + intersect

			if len(intersect)==1:
				Cell_table[j][i] = "$"+str(KLP(e,intersect[0]))+"$"
			else:
				Cell_table[j][i] = "\\begin{tabular}{@{}l@{}}" + (" \\"+"\\ ").join(["$"+str(KLP(e,x))+"$" for x in intersect]) + "\\end{tabular}"



	print("\\begin{tabular}{|%s|}"%("|").join(["l"]*(len(DufCell)+1) ))

	print("\\hline")
	print("a=%d &"%a(DufCell[0]) + " & ".join([ " %s "%list_descends(DL(DufCell[i])) for i in range(len(DufCell))] ) )
	print(" \\"+"\\")

	print("\\hline")
	for i in range(len(Cell_table)):
		row = Cell_table[i]


		print(" & ".join([list_descends(DR(DufCell[i]))] + row) + " \\"+"\\"+" \\hline")

	print("\\end{tabular}\n")
	
	



def penultimate_two_cell_latex_v2():
    
    result = ''

    Cell = two_cell(w0*s1)
    DufCell = [d for d in Cell if d in D]
    DufCell.sort(key=sort_mult)
    
    # Some manual sorting:
#    temp = DufCell[1]
#    DufCell[1] = DufCell[3]
#    DufCell[3] = temp

    DufCell_dict = {}
    for d in DufCell:
        Ld = L_cell(d)
        Rd = [x.inverse() for x in Ld]
        DufCell_dict[d]=(Rd,Ld)

    Cell_table = [[[] for x in DufCell] for x in DufCell]
    for i in range(len(DufCell)):
        for j in range(len(DufCell)):

            intersect = [x for x in DufCell_dict[DufCell[i]][0] if x in DufCell_dict[DufCell[j]][1] ]
            intersect.sort(key=l)
            if i == j:
                intersect.remove(DufCell[i])
                intersect = [DufCell[i]] + intersect

            if len(intersect)==1:
                Cell_table[j][i] = "$"+convert_to_123(intersect[0])+" \\ (%s)"%KLP(e,intersect[0])+"$"
            else:
                Cell_table[j][i] = "\\begin{tabular}{@{}l@{}}" + (" \\"+"\\ ").join(["$"+convert_to_123(x)+" \\ (%s)"%KLP(e,x)+"$" for x in intersect]) + "\\end{tabular}"



    result += ("\\begin{tabular}{|%s|}"%("|").join(["l"]*len(DufCell) ))
    result += ("\\hline ")
    for row in Cell_table:

        result += (" & ".join(row) + " \\"+"\\"+" \\hline")

    result += ("\\end{tabular}")
    
    print(result)
	






### KLP's for the penultimate cell in the table format

def sort_asc(w):
    return eval(convert_to_123(list(AL(w))[0]))

def l(w):
    return w.length()

def list_descends(S):
    result = [convert_to_123(s) for s in S]
    result.sort()
    return ", ".join(result)

def coeff_one(i):
    if i != 1:
        return str(i)
    else:
        return ""
    
def KLPv_e(w):
    p = KLP(e,w)
    l = w.length()
    pv = sum(p.coefficients()[i] * q^(l - 2*p.exponents()[i]) for i in range(len(p.exponents())))
    
    # Format pv nicely for latex:
    
    return " + ".join([ coeff_one(pv.coefficients()[i])+"v^{%s}"%pv.exponents()[i] for i in range(len(pv.exponents())) ])

def penultimate_two_cell_latex():
	
	Cell = two_cell(s1*w0)
	DufCell = [d for d in Cell if d in Duflo_Involutions()]
	DufCell.sort(key=sort_asc)

	DufCell_dict = {}
	for d in DufCell:
		Ld = L_cell(d)
		Rd = [x.inverse() for x in Ld]
		DufCell_dict[d]=(Rd,Ld)

	Cell_table = [[[] for x in DufCell] for x in DufCell]
	for i in range(len(DufCell)):
		for j in range(len(DufCell)):

			intersect = [x for x in DufCell_dict[DufCell[i]][0] if x in DufCell_dict[DufCell[j]][1] ]
			intersect.sort(key=l)
			if i == j:
				intersect.remove(DufCell[i])
				intersect = [DufCell[i]] + intersect

			if len(intersect)==1:
				Cell_table[j][i] = "$"+KLPv_e(intersect[0])+"$"
			else:
				Cell_table[j][i] = "\\begin{tabular}{@{}l@{}}" + (" \\"+"\\ ").join(["$"+KLPv_e(x)+"$" for x in intersect]) + "\\end{tabular}"



	print("\\begin{tabular}{|l||%s|}"%("|").join(["l"]*(len(DufCell)) ))

	print("\\hline")
	print( " $%s_%d$ & "%(CartanType(W)[0],CartanType(W)[1]) + " & ".join([ " $%s$ "%convert_to_123(AL(DufCell[i])) for i in range(len(DufCell))] ) )
	print(" \\"+"\\")

	print("\\hline"*2)
	for i in range(len(Cell_table)):
		row = Cell_table[i]


		print(" & ".join([" $%s$ "%convert_to_123(AR(DufCell[i]))] + row) + " \\"+"\\"+" \\hline")

	print("\\end{tabular}\n")
	
	