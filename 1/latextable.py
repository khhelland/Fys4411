

E = open("6particles3shells.mat","r")
T = open("6particles3shells.tex","w")

endl = "\\\ \n"
sep = "&"
T.write("$\\left(\\begin{matrix}")
next(E)
next(E)
for line in E:
    cols = line.split()
    T.write("%.2f" %(float(cols[0])))
    for col in cols[1:]:
        #print type(val)
        T.write(sep)
        val = float(col)
        T.write("%.2f" %(val))
    T.write(endl)

T.write("\\end{matrix}\\right)$")
