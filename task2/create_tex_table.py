import subprocess

delta = [0.1, 0.01, 0.001]
degs=[1, 2, 3]
mus = [0.001]
Cs = [1.4]

for C in Cs:
    for mu in mus:
        for h in [0.01]:
            for tau in [0.001]:
                print("\\begin{tabular}{ |l|l|l|l|l|}")
                print("\\hline")
                if C == 1.4:
                    print("\\multicolumn{5}{|c|}{$\\mu =", str(mu), ", p(\\rho) = \\rho^{1.4},", "h = ", str(h), ", \\tau = ", str(tau), "$}\\\\")
                else:
                    print("\\multicolumn{5}{|c|}{$\\mu =", str(mu), ", p(\\rho) =", str(C), " \\rho,", "h = ", str(h), ", \\tau = ", str(tau), "$}\\\\")
                print("\\hline")
                result = subprocess.run([
                    "./a.out", 
                    str(h), 
                    str(tau), 
                    str(mu), 
                    str(C), 
                    "1" if C != 1.4 else "0", 
                    "-1",
                    "25",
                    "0"
                ], capture_output=True, text=True)
                stab_time = 1.
                stab_norm = None
                delta_massa = None
                stab_index = 1
                for line in result.stdout.split('\n'):
                    words = line.split()
                    for word in words:
                        if word == "stab_norm":
                            stab_norm = words[2]
                        if word == "stab_time":
                            stab_time = float(words[2])
                        if word == "delta_massa":
                            delta_massa = words[2]
                        if word == "stab_index":
                            stab_index = int(words[2])
                stab_indexes = [int(stab_index / 4), int(stab_index / 2), int(3 * stab_index / 4), stab_index]
                results = {}
                for i, index in enumerate(stab_indexes):
                    if (i == len(stab_indexes) - 1):
                        results[(i, "stab_norm")] = str(stab_norm)
                        results[(i, "delta_massa")] = str(delta_massa)
                        break
                    result = subprocess.run([
                        "./a.out", 
                        str(h), 
                        str(tau), 
                        str(mu), 
                        str(C), 
                        "1" if C != 1.4 else "0", 
                        str(index),
                        "25",
                        "0"
                    ], capture_output=True, text=True)
                    for line in result.stdout.split('\n'):
                        words = line.split()
                        for word in words:
                            if word in words:
                                if word == "stab_norm":
                                    results[(i, "stab_norm")] = words[2]
                                if word == "delta_massa":
                                    results[(i, "delta_massa")] = words[2]
                print("$ $ & $n_{st}/4$ & $n_{st}/2$ & $3n_{st}/4$ & $n_{st}, (", str(stab_time), ")$\\\\")
                print("\\hline")
                print("$norm$", end="")
                for i in range(len(stab_indexes)):
                    print(" & $", results[(i, "stab_norm")], "$", sep="", end="" )
                print("\\\\")
                print("\\hline")
                print("$\\Delta_{massa}$", end="")
                for i in range(len(stab_indexes)):
                    print(" & $", results[(i, "delta_massa")], "$", sep="",end="")
                print("\\\\")
                print("\\hline")
                print("\\end{tabular}")
                print("\\centering")
                print("\\begin{tabular}{ |l|l|}")
                print("\\hline")
                if C == 1.4:
                    print("\\multicolumn{2}{|c|}{$\\mu =", str(mu), " ,p(\\rho) = \\rho^{1,4}$}\\\\")
                else:
                    print("\\multicolumn{2}{|c|}{$\\mu =", str(mu), " ,p(\\rho) = ", str(C), "\\rho$}\\\\")
                print("\\hline")
                print("$ $ & tau=", str(tau), ",h=", str(h), "\\\\")
                print("\\hline")
                for k in degs:
                    print("$h - h^{", str(k), "}$ ", end="", sep="")
                    norms = {}
                    result = subprocess.run([
                        "./a.out", 
                        str(h), 
                        str(tau), 
                        str(mu), 
                        str(C), 
                        "1" if C != 1.4 else "0", 
                        str(int(stab_index / 10)),
                        "25",
                        str(k)
                    ], capture_output=True, text=True)
                    for line in result.stdout.split('\n'):
                        words = line.split()
                        for word in words:
                            if word == "C_norm_H":
                                norms[("C_norm_H")] = words[2]
                            if word == "L_norm_H":
                                norms[("L_norm_H")] = words[2]
                            if word == "W_norm_H":
                                norms[("W_norm_H")] = words[2]
                    print("& $", norms[("C_norm_H")], "$", " ", sep="", end="")
                    print("\\\\")
                    print("& $", norms[("L_norm_H")], "$", " ", sep="", end="")
                    print("\\\\")
                    print("& $", norms[("W_norm_H")], "$", " ", sep="", end="")
                    print("\\\\")
                    print("\\hline")
                print("\\end{tabular}")
                print()



