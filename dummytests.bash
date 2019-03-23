#!/bin/bash

cp ./L1024_K26_B50000000/GRCh38.align ./
echo "ARCHIVO L1024_K26_B50000000/GRCh38.align COPIADO" >> "PruebaDummy.txt"
echo "PRUEBAS SECUENCIAL" >> "PruebaDummy.txt"
./EncoderSec >> "PruebaDummy.txt"
echo "PRUEBAS ESTATICA" >> "PruebaDummy.txt"
for N in 2 4 6 8 16 32 48; do
    echo "PRUEBAS ESTATICA CON N = $N" >> "PruebaDummy.txt"
    ./EncoderEst -N $N >> "PruebaDummy.txt"
done
echo "PRUEBAS DINAMICAS" >> "PruebaDummy.txt"
echo "PRUEBAS DINAMICAS CON C = 10000" >> "PruebaDummy.txt"
for N in 2 4 6 8 16 32 48; do
    echo "PRUEBAS DINAMICAS CON N = $N" >> "PruebaDummy.txt"
    ./EncoderDin -N $N -C 10000 >> "PruebaDummy.txt"
done
echo "PRUEBAS DINAMICAS CON C = 20000" >> "PruebaDummy.txt"
for N in 2 4 6 8 16 32 48; do
    echo "PRUEBAS DINAMICAS CON N = $N" >> "PruebaDummy.txt"
    ./EncoderDin -N $N -C 20000 >> "PruebaDummy.txt"
done