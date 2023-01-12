#!/bin/bash

for i in {21..1000}
do
    ./benchmark -N 1000 -mu 0.2 -k 10 -maxk 100
    mv ./network.dat ./sbm/$i.dat
done