import subprocess
import functools

print = functools.partial(print, flush=True)

mus = [0.1]
Cs = [10]
#Cs = [1.4]
init_data = [1, 2, 3, 4, 5, 6, 7]

for mu in mus:
    for C in Cs:
        for h in [0.001]:
            for tau in [0.0001]:
                print("\\begin{tabular}{*{", len(init_data) + 2, "}{|c}}")
                print("\\hline")
                if C == 1.4:
                    print("\\multicolumn{", len(init_data) + 1, "}{|c|}{$","\\mu = ", mu, ", p(\\rho) = \\rho^{1.4}", ", h = ", h, ", \\tau = ", tau , "$} \\\\")
                else:
                    print("\\multicolumn{", len(init_data) + 1, "}{|c|}{$","\\mu = ", mu, ", p(\\rho) = ", C, "\\cdot", "\\rho", ", h = ", h, ", \\tau = ", tau , "$} \\\\")
                print("\\hline")
                print("\\diagbox{$\\rho$}{u}", end="")
                for u in init_data:
                    print("& ", u, end=" ")
                print("\\\\")
                print("\\hline")
                
                for rho in init_data:
                    print("$", rho, "$", end="")
                    for u in init_data:
                        result = subprocess.run([
                                "./a.out",
                                str(h),
                                str(tau),
                                str(mu),
                                str(C),
                                "1" if C != 1.4 else "0",
                                str(u),
                                str(rho) 
                            ], capture_output=True, text=True)
                        stab_time = None
                        for line in result.stdout.split('\n'):
                            words = line.split()
                            for word in words:
                                if word == "stab_time":
                                    stab_time = words[2]
                        print(" & $", stab_time, "$", end="")
                    print("\\\\")
                    print("\\hline")
                print("\\end{tabular}")
                print()
