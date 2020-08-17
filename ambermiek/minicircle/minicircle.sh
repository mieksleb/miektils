#!/bin/bash

export PATH=$PATH:.

echo "enter raw sequence"

read seq_file

# call the python script ambermiek.py which takes a raw sequence file and generates everything needed to create a circular piece of dna with a thymine dimer
N=$(python $AMBERMIEK/minicircle/ambermiek.py $seq_file)  # N stores the number of base pairs
M=$[$N+2]                                                 # M is the number of base pairs for circularisation process (extra g and c)

# make a folder to store all generated files
dir_name="TT_${N}"
mkdir $dir_name
mv dv${N}.nab $dir_name/dv${N}.nab
cd $dir_name

# find M which is twice the number of base pairs in the system
P=$[$M * 2]

$AMBERHOME/bin/nab -O "dv${N}.nab"
mv a.out dv${N}.exe

./dv${N}.exe

egrep -v H dv${N}linear.pdb > dv${N}linearnoH.pdb

echo "enter a linking number"

read L

$AMBERMIEK/minicircle/circularise.exe -n $P -t $L < dv${N}linearnoH.pdb > dvec${N}t${L}.pdb

cd ../


dir_name2="TT_${N}_${L}"                # now we name master directory to include the linking number
mv $dir_name $dir_name2 

cd $dir_name2

mkdir BUILD 
mkdir SIMULATION
cd SIMULATION

echo "Enter the base pair index i such that the thymine dimer is between bases i and i+1"

# first we must check that there are in fact thymines at the locations specified



# note that we index base pairs starting from 1, oxDNA topologies index from 0.

read TT

python $AMBERMIEK/minicircle/ambermiek2.py $N $L $TT


#now we must run the tleap scripts, although first we have to move some scripts from the minicircle home directory

cp $AMBERMIEK/minicircle/parmBSC1.lib .
cp $AMBERMIEK/minicircle/leaprc.ff14SBcirc .
cp $AMBERMIEK/minicircle/frcmod.parmbsc1 .
cp $AMBERMIEK/minicircle/frcmod.TT .
cp $AMBERMIEK/minicircle/TT.prepin .
cp $AMBERMIEK/minicircle/$dir_name2/* .


source $AMBERHOME/amber.sh
# leapscript_0 created the non-TD pdb
tleap -f leapscript_0

# we then run ambermiek3.py which replaces the residue types from DT to TT for $TT and $TT+1
python $AMBERMIEK/minicircle/ambermiek3.py dv${N}t${L}.pdb dv${N}t${L}tt.pdb $TT
# output is dvNtLtt.pdb

tleap -f leapscript_1

tleap -f leapscript_2

# call submission.py to generate the submission scripts
python $AMBERMIEK/minicircle/submission.py $N $L $AMBERMIEK/minicircle/sub_arcus.sh $AMBERMIEK/minicircle/sub_csd3.sh


# now we create 2 directories, one for the TD simulation, and one for the undamaged dna simulation. Then move appropriate files into each one.
mkdir MTTD
mkdir MDNA
mv *.in MTTD
mv dv${N}t${L}tt.* MTTD

mv dv${N}t${L}.* MDNA
cp MTTD/*.in MDNA

pwd
ls


mv *tt MTTD
mv *tt.sh MTTD
mv *.sh MDNA


# then clean up all directories, moving everything unmoved into the build directory
cd ../
mv *.* BUILD
