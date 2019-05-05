#!/bin/bash

# FIXING K
for B in 10000000 50000000 100000000
do
    cp ./FPL1024\_K26_B$B/GRCh38.align ./
    echo "ARCHIVO L1024_K26_B$B/GRCh38.align COPIADO" >> "secuencialTest.txt"
    echo "PRUEBAS SECUENCIAL" >> "secuencialTest.txt"
    ./EncoderSec >> "secuencialTest.txt"
    rm GRCh38.align
    echo "TERMINADO L1024_K26_B$B/"
done

for B in 10000000 50000000 100000000
do
    cp ./FPL1024\_K128_B$B/GRCh38.align ./
    echo "ARCHIVO L1024_K128_B$B/GRCh38.align COPIADO" >> "secuencialTest.txt"
    echo "PRUEBAS SECUENCIAL" >> "secuencialTest.txt"
    ./EncoderSec >> "secuencialTest.txt"
    rm GRCh38.align
    echo "TERMINADO L1024_K128_B$B/"
done

for B in 10000000 50000000 100000000
do
    cp ./FPL1024\_K256_B$B/GRCh38.align ./
    echo "ARCHIVO L1024_K256_B$B/GRCh38.align COPIADO" >> "secuencialTest.txt"
    echo "PRUEBAS SECUENCIAL" >> "secuencialTest.txt"
    ./EncoderSec >> "secuencialTest.txt"
    rm GRCh38.align
    echo "TERMINADO L1024_K256_B$B/"
done
