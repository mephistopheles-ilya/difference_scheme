#! /bin/bash

delta=(0.1 0.01 0.001 0.0001)
mus=(0.001 0.01 0.1)
Cs=(1 10 100)

for C in "${Cs[@]}"; do
    for mu in "${mus[@]}"; do
        for h in "${delta[@]}"; do
            for tau in "${delta[@]}"; do
                echo "===================="
                ./a.out "$h" "$tau" "$mu" "$C" 1 0
                echo "===================="
                echo
            done
        done
    done
done

for mu in "${mus[@]}"; do
    for h in "${delta[@]}"; do
        for tau in "${delta[@]}"; do
            echo "===================="
            ./a.out "$h" "$tau" "$mu" 1.4 0 0
            echo "==================="
            echo
        done
    done
done
