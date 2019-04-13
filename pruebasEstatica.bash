#!/bin/bash

# FIXING K
for B in 10000000 50000000 100000000
do
    for L in 256 512 1024
    do 
        cp ./L$L\_K26_B$B/GRCh38.align ./
        echo "ARCHIVO L'$L'_K26_B$B/GRCh38.align COPIADO" >> "estaticaTest.txt"
        ./EncoderEst -N 88 >> "estaticaTest.txt"
        rm GRCh38.align
        echo "TERMINADO L'$L'_K26_B$B/"
    done
done

for B in 10000000 50000000 100000000
do
    for L in 256 512 1024
    do 
        cp ./L$L\_K128_B$B/GRCh38.align ./
        echo "ARCHIVO L'$L'_K128_B$B/GRCh38.align COPIADO" >> "estaticaTest.txt"
        ./EncoderEst -N 88 >> "estaticaTest.txt"
        rm GRCh38.align
        echo "TERMINADO L'$L'_K128_B$B/"
    done
done