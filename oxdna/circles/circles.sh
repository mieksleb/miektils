#!/bin/bash

echo "Enter a circle length in base pairs"
read N
echo "Enter a linking number"
read L

dir="circle_$N-Lk$L"
mkdir $dir
cd $dir

python3 $MIEKTILS/oxdna/circles/circular_double_generator.py $N
python2 $OXDNA/oxDNA/UTILS/generate.py $N sequence_$N $L

mv generated.dat circle_$N-Lk$L.dat
mv generated.top circle_$N-Lk$L.top

#python3 $MIEKTILS/oxdna/circles/string_forces.py circle_$N-Lk$L.dat
python $MIEKTILS/oxdna/circles/CPU_input.py circle_$N-Lk$L.top circle_$N-Lk$L.dat trajectory_$N-Lk$L.dat

mv input.dat cpu_input

python $MIEKTILS/oxdna/circles/GPU_input.py circle_$N-Lk$L.top circle_$N-Lk$L.dat trajectory_$N-Lk$L.dat
mv input.dat gpu_input
