import subprocess

delta = [0.01, 0.001]
degs=[1, 2, 3, 0]
mus = [0.1, 0.01]
Cs = [1, 10, 1.4]

for C in Cs:
    for mu in mus:
        print("\\begin{tabular}{ |l|l|l|}")
        print("\\hline")
        if C == 1.4:
            print("\\multicolumn{3}{|c|}{$\\mu =", str(mu), ",p(\\rho) = \\rho^{1,4}$}\\\\")
        else:
            print("\\multicolumn{3}{|c|}{$\\mu =", str(mu), ",p(\\rho) =", str(C), "\\rho$}\\\\")
        print("\\hline")
        #print("$ $ & tau=h=$0.1$& tau=h=$0.01$ & tau=h=$0.001$\\\\")
        print("$ $ & tau=h=$0.01$ & tau=h=$0.001$\\\\")
        print("\\hline")
        for k in degs:
            if k != 0:
                print("$h - h^{", str(k), "}$ ", end="", sep="")
            else:
                print("$h - \\rho$ ", end="", sep="")
            norms = {}
            for h in delta:
                result = subprocess.run([
                    "./a.out", 
                    str(h), 
                    str(h), 
                    str(mu), 
                    str(C), 
                    "1" if C != 1.4 else "0", 
                    str(k)
                ], capture_output=True, text=True)
                for line in result.stdout.split('\n'):
                    words = line.split()
                    for word in words:
                        if word == "C_norm_H":
                            norms[("C_norm_H", str(h))] = words[2]
                        if word == "L_norm_H":
                            norms[("L_norm_H", str(h))] = words[2]
                        if word == "W_norm_H":
                            norms[("W_norm_H", str(h))] = words[2]
            for h in delta:
                print("& $", norms[("C_norm_H", str(h))], "$", " ", sep="", end="")
            print("\\\\")
            for h in delta:
                print("& $", norms[("L_norm_H", str(h))], "$", " ", sep="", end="")
            print("\\\\")
            for h in delta:
                print("& $", norms[("W_norm_H", str(h))], "$", " ", sep="", end="")
            print("\\\\")
            print("\\hline")
        print("\\end{tabular}")
        print()

