import subprocess
import functools

print = functools.partial(print, flush=True)

mus = [0.1, 0.01, 0.001]
#Cs = [1, 10, 100]
Cs = [1.4]

for h in [0.001]:
    for tau in [0.0001]:
        print("\\begin{tabular}{*{20}{|c}}")
        print("\\hline")
        print("$k$", end="")
        for C in Cs:
            for mu in mus:
                if C != 1.4:
                    print(" & $\\begin{array}{c}\\mu = ", mu, "\\\\", "p(\\rho) = ", C, "\\cdot \\rho", "\\\\", "\\tau = ", tau, "\\\\", "h = ", h, "\\end{array}$", sep="", end="")
                else:
                    print(" & $\\begin{array}{c}\\mu = ", mu, "\\\\", "p(\\rho) = \\rho^{1.4}", "\\\\", "\\tau = ", tau, "\\\\", "h = ", h, "\\end{array}$", sep="", end="")
        print("\\\\")
        print("\\hline")

        for k in range(1, 20):
            if k <= 10:
                print("$", k, "$", end="")
            else:
                print("$10 +\\frac{",k - 10, "\\cdot M}{10}$", end="") 
            for C in Cs:
                for mu in mus:
                    result = subprocess.run([
                        "./a.out",
                        str(h),
                        str(tau),
                        str(mu),
                        str(C),
                        "1" if C != 1.4 else "0",
                        str(k),
                        "28"
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
