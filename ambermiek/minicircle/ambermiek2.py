#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 14:19:26 2020

This script creates all the tleapfiles and amber input files required, including a SLURM submission script



@author: michaelselby
"""
import sys

leap0 = "leapscript_0"
leap1 = "leapscript_1"
leap2 = "leapscript_2"

N = sys.argv[1]
L = sys.argv[2]
TT = sys.argv[3]


leapfile0 = open(leap0, "w")
leapfile1 = open(leap1, "w")
leapfile2 = open(leap2, "w")

leapfile0.write(
'source leaprc.ff14SBcirc'+'\n'
'addAtomTypes {' +'\n'
'  { "C1" "C" "sp2" }' +'\n'
'  { "C2" "C" "sp2" }' +'\n'
'  { "CI" "C" "sp3" }' +'\n'
'  { "CE" "C" "sp3" }' +'\n'
' }'+'\n'
' loadoff parmBSC1.lib' +'\n'
' loadamberparams frcmod.parmbsc1' +'\n'
"addPdbAtomMap {{OP1 O1P}{OP2 O2P}{H5' H5'1}{H5'' H5'2}{H2' H2'1}{H2'' H2'2}}" +'\n'
'set default PBradii mbondi3' +'\n'

'dna = loadpdb '+'dvec'+str(N)+'t'+str(L)+'.pdb' +'\n'

'savepdb dna '+ 'dv'+str(N)+'t'+str(L)+'.pdb'+'\n'

'quit'


)
leapfile0.close()


# now we generate the second leapscript
leapfile1.write(
        
'source leaprc.ff14SBcirc'  +'\n'
'addAtomTypes {'+'\n'
'  { "C1" "C" "sp2" }' +'\n'
'  { "C2" "C" "sp2" }' +'\n'
'  { "CI" "C" "sp3" }' +'\n'
'  { "CE" "C" "sp3" }' +'\n'
' }' +'\n'
' loadoff parmBSC1.lib' +'\n'
' loadamberparams frcmod.parmbsc1' +'\n'
"addPdbAtomMap {{OP1 O1P}{OP2 O2P}{H5' H5'1}{H5'' H5'2}{H2' H2'1}{H2'' H2'2}}" +'\n'
'set default PBradii mbondi3' +'\n'

'dna = loadpdb ' + 'dv'+str(N)+'t'+str(L)+'.pdb' +'\n'

'bond dna.1.P dna.'+str(N)+".O3' S" +'\n'
'bond dna.'+str(int(N)+1)+'.P dna.'+str(2*int(N))+".O3' S" +'\n'
'saveamberparm dna '+'dv'+str(N)+'t'+str(L)+'.top dv'+str(N)+'t'+str(L)+'.crd'+'\n'

'quit'

)

leapfile1.close()




leapfile2.write(

'source leaprc.ff14SBcirc'  +'\n'
'addAtomTypes {'+'\n'
'  { "C1" "C" "sp2" }' +'\n'
'  { "C2" "C" "sp2" }' +'\n'
'  { "CI" "C" "sp3" }' +'\n'
'  { "CE" "C" "sp3" }' +'\n'
' }' +'\n'
' loadoff parmBSC1.lib' +'\n'
' loadamberprep TT.prepin' +'\n'
' loadamberparams frcmod.parmbsc1' +'\n'
' loadamberparams frcmod.TT' +'\n'
"addPdbAtomMap {{OP1 O1P}{OP2 O2P}{H5' H5'1}{H5'' H5'2}{H2' H2'1}{H2'' H2'2}}" +'\n'
'set default PBradii mbondi3' +'\n'

'dna = loadpdb dv'+str(N)+'t'+str(L)+'tt.pdb' +'\n'

"bond dna.1.P dna."+str(N)+".O3' S" +'\n'
"bond dna."+str(int(N)+1)+".P dna."+str(2*int(N))+".O3' S" +'\n'
'bond dna.'+str(TT)+'.C5 dna.'+str(int(TT)+1)+'.C5 S' +'\n'
'bond dna.'+str(TT)+'.C6 dna.'+str(int(TT)+1)+'.C6 S'+'\n'

'saveamberparm dna dv'+str(N)+'t'+str(L)+'tt.top dv'+str(N)+'t'+str(L)+'tt.crd'+'\n'

'quit'
)

leapfile2.close()


min1 = "min1.in"
md1 = "md1.in"
md2 = "md2.in"
md3 = "md3.in"
md4 = "md4.in"
md5 = "md5.in"

min1file = open(min1, "w")
md1file = open(md1, "w")
md2file = open(md2, "w")
md3file = open(md3, "w")
md4file = open(md4, "w")
md5file = open(md5, "w")


min1file.write(

"initial minimization whole system" +'\n'
" &cntrl" +'\n'
"  imin=1," +'\n'
"  maxcyc=10000," +'\n'
"  ncyc=2500," +'\n'
"  ntb=0," +'\n'
"  ntp=0," +'\n'
"  ntr=0," +'\n'
"  igb=1," +'\n'
"  cut=1000.," +'\n'
" /"      
        
)
min1file.close()


md1file.write(
        
"MD on the water and ions about the DNA" +'\n'
" &cntrl" +'\n'
'        ntf=2, ntc=2, ntb=0,cut=1000.,' +'\n'
'        nstlim=5000, dt=0.002,' +'\n'
'        tempi=100.0, temp0=300.0,' +'\n'
'        ntt =3, gamma_ln=0.01,' +'\n'
'        imin=0, irest=0, ntx=1,' +'\n'
'        igb=8, gbsa=0, saltcon=0.2,' +'\n'
'        ntr=1,' +'\n'
'        restraint_wt=50.0,' +'\n'
"        restraintmask=':1-"+str(2*int(N))+"'" +'\n'
' &end'
            
)
md1file.close()


md2file.write(

'MD on the water and ions about the DNA'+'\n'
' &cntrl' +'\n'
'        ntf=2, ntc=2, ntb=0, ,cut=1000.,' +'\n'
'        nstlim=50000, dt=0.002, ntpr=500,' +'\n'
'        temp0=300.0, ntt=3, gamma_ln=0.01,' +'\n'
'        imin=0, irest=1, ntx=5,' +'\n'
'        igb=8, gbsa=0, saltcon=0.2,' +'\n'
'        ntr=1, ' +'\n'
'        restraint_wt=10.0,' +'\n'
"        restraintmask=':1-"+str(2*int(N))+"'"+'\n'
' &end '       
              
)
md2file.close()



md3file.write(
        
'MD on the water and ions about the DNA'+'\n'
' &cntrl' +'\n'
'        ntf=2, ntc=2, ntb=0, ,cut=1000.,' +'\n'
'        nstlim=100000, dt=0.002, ntpr=500,' +'\n'
'        temp0=300.0, ntt=3, gamma_ln=0.01,' +'\n'
'        imin=0, irest=1, ntx=5,' +'\n'
'        igb=8, gbsa=0, saltcon=0.2,' +'\n'
'        ntr=1, ' +'\n'
'        restraint_wt=1.0,' +'\n'
"        restraintmask=':1-"+str(2*int(N))+"'"+'\n'
' &end '       
              
)
md3file.close()



md4file.write(
        
'MD on the water and ions about the DNA'+'\n'
' &cntrl' +'\n'
'        ntf=2, ntc=2, ntb=0, ,cut=1000.,' +'\n'
'        nstlim=5000000, dt=0.002, ntpr=5000,' +'\n'
'        ntwx=5000,ntwe=5000,ntwr=5000,' +'\n'
'        temp0=300.0, ntt=3, gamma_ln=0.01,' +'\n'
'        imin=0, irest=1, ntx=5,' +'\n'
'        igb=8, gbsa=0, saltcon=0.2,' +'\n'
'        ntr=0, nmropt=1,'+'\n'
' &end '       
              
)
md4file.close()



md5file.write(
                
'MD on the water and ions about the DNA'+'\n'
' &cntrl'+'\n'
'        ntf=2, ntc=2, ntb=0, ,cut=1000.,'+'\n'
'        nstlim=5000000, dt=0.002, ntpr=5000,'+'\n'
'        ntwx=5000,ntwe=5000,ntwr=5000,'+'\n'
'        temp0=300.0, ntt=3,gamma_ln=0.01,'+'\n'
'        imin=0, irest=1, ntx=5,'+'\n'
'        igb=8, gbsa=0, saltcon=0.2,'+'\n'
'        ntr=0,'+'\n'
' &end'

)
md5file.close()


