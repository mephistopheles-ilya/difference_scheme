file_name = "res.txt"

table = {}

with open(file_name, 'r') as file:
    key = []
    nums = []
    mu = None
    pp = None
    for line in file:
        words = line.split()
        if len(words) > 9 and words[9] == 'pp' and words[6] == 'mu':
            mu = words[8]
            pp = words[11]
        if len(words) > 4 and words[0] == 'h' and words[3] == 'tau':
            key.append(words[2])
            key.append(words[5])
        if len(words) > 2 and words[0] == 'C_norm_H':
            nums.append(words[2])
        if len(words) > 2 and words[0] == 'L_norm_H':
            nums.append(words[2])
        if len(words) > 2 and words[0] == 'W_norm_H':
            nums.append(words[2])
        if len(words) > 2 and words[0] == "time":
            nums.append(words[2])
            key.append(mu)
            key.append(pp)
            table[tuple(key)] = nums.copy()
            nums.clear()
            key.clear()


for pp in ('1.000000', '10.000000', '100.000000', '1.400000'):
    for mu in ('0.100000', '0.010000', '0.001000'):
        print("\\begin{tabular}{ |l|l|l|l|l| }")
        print("\\hline")
        if (pp == '1.400000'):
            print("\\multicolumn{5}{|c|}{$\\mu =", mu.rstrip('0').rstrip('.'), ", p(\\rho) = \\rho^{1.4}$}\\\\")
        else:
            print("\\multicolumn{5}{|c|}{$\\mu =", mu.rstrip('0').rstrip('.'), ", p(\\rho) =", pp.rstrip('0').rstrip('.'), "\\rho$}\\\\")
        print("\\hline")
        print("$\\tau\\setminus h$ & $0.1$ & $0.01$ & $0.001$ & $0.0001$\\\\")
        print("\\hline")
        for tau in ("0.100000", "0.010000", "0.001000", "0.000100"):
            print("$", tau, "$", sep="", end=" ")
            for i in (0, 1, 2, 3):
                for h in ("0.100000", "0.010000", "0.001000", "0.000100"):
                    print("&", " ", "$", table[(h, tau, mu, pp)][i], "$", end=" ", sep="")
                print("\\\\")
            print("\\hline")
        print("\\end{tabular}")
        print()






