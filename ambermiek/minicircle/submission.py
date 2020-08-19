#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 11:25:44 2020


This script generates the submission scripts required to run amber all atom simulations on either arcus or csd3
but could easily be adapted for other clusters/machines.

All follow the same structure but can can be a mixture of sander, pmemd.cuda, or pmemed.cuda.MPI. This variable
is known as the string 'application'.

The cluster specific parts of the submission scripts will be generated using generic scripts for each cluster which must
be stored in directory $AMBERMIEK

@author: michaelselby
"""
import sys

sub_file_csd3 = "sub_csd3.sh"
sub_file_arcus = "sub_arcus.sh"



N = sys.argv[1]
L = sys.argv[2]
#sub_dir = sys.argv[3]
sub_file_arcus = sys.argv[3]
sub_file_csd3 = sys.argv[4]

TT_string = str(N)+'t'+str(L)+'tt'
hom_string = str(N)+'t'+str(L)


def command_block(input_string, application, cluster_name, template, minimization=False):
    # this function writes the machine independent command for the amber to process.
    
    sub_file_name = "sub_"+cluster_name+"_"+input_string+".sh"
    input_string = "dv" + input_string
    template_file = open(template,"r").readlines()
    sub_file = open(sub_file_name,"w")
    for line in template_file:
        if "jobinput" in line:
            line = "#SBATCH -J "+input_string
        elif "amberinput" in line:
            if minimization==True:
                if application=="pmemd.cuda.MPI":
	                sub_file.write(
	                    'pmemd.cuda -O -i min1.in -o min1.out -inf min1.inf \\' + '\n' 
	                    '\t\t-c '+input_string+'.crd -ref '+input_string+'.crd -r '+input_string+'.min1 \\' + '\n'
	                    '\t\t-p '+input_string+'.top -x '+input_string+'.min1.x -e '+input_string+'.min1.ene'+'\n'
	                    '\n'
	                    )
                else:
                    sub_file.write(
                            application+' -O -i min1.in -o min1.out -inf min1.inf \\' + '\n'
                            '\t\t-c '+input_string+'.crd -ref '+input_string+'.crd -r '+input_string+'.min1 \\' + '\n'
                            '\t\t-p '+input_string+'.top -x '+input_string+'.min1.x -e '+input_string+'.min1.ene'+'\n'
                            '\n'
                            )
           
                
                
            sub_file.write(

                application+' -O -i md1.in -o md1.out -inf md1.inf \\' + '\n' 
                '\t\t-c '+input_string+'.min1 -ref '+input_string+'.min1 -r '+input_string+'.md1 \\' + '\n'
                '\t\t-p '+input_string+'.top -x '+input_string+'.md1.x -e '+input_string+'.md1.ene'+'\n'
                '\n'+
                application+' -O -i md2.in -o md2.out -inf md2.inf \\' + '\n'
                '\t\t-c '+input_string+'.md1 -ref '+input_string+'.md1 -r '+input_string+'.md2 \\' + '\n'
                '\t\t-p '+input_string+'.top -x '+input_string+'.md2.x -e '+input_string+'.md2.ene'+'\n'
                '\n'+
                application+' -O -i md3.in -o md3.out -inf md3.inf \\' + '\n'
                '\t\t-c '+input_string+'.md2 -ref '+input_string+'.md2 -r '+input_string+'.md3 \\' + '\n'
                '\t\t-p '+input_string+'.top -x '+input_string+'.md3.x -e '+input_string+'.md3.ene'+'\n'
                '\n'+
                application+' -O -i md5.in -o md5.out -inf md5.inf \\' + '\n'
                '\t\t-c '+input_string+'.md3 -ref '+input_string+'.md3 -r '+input_string+'.md5 \\' + '\n'
                '\t\t-p '+input_string+'.top -x '+input_string+'.md5.x -e '+input_string+'.md5.ene' +'\n'

                        )
               

            
            
            
        elif "cpptrajinput" in line:
            line = "\n"
            sub_file.write(
                "cpptraj.MPI -p "+input_string+".top -y "+input_string+".md5.x -x "+input_string+"_output-file.pdb"
                )
        else:
            sub_file.write(line)            
    sub_file.close()
    
            
    return


command_block(TT_string,'pmemd.cuda.MPI', 'arcus', sub_file_arcus, True)
command_block(hom_string,'pmemd.cuda.MPI', 'arcus', sub_file_arcus, True)

command_block(TT_string,'pmemd.cuda.MPI', 'csd3', sub_file_csd3, False)
command_block(hom_string,'pmemd.cuda.MPI', 'csd3', sub_file_csd3, False)





                 
