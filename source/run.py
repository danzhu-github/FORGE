#!/usr/bin/env python

import main

pathin ='/your_directory/INPUT/'  # the path where the input files are stored
pathout='/your_directory/OUTPUT/' # the path where the output file will be generated
res='year'  # temporal resolution in the output file: 'year','month','daily'
run_year=300  # run how many years

pars=[0.015,0.1,1.0,200,2.5] # default parameter values,corresponding to ΦV,ΦA,emax,c,q in Table 1. 
phi_v =pars[0]
phi_a =pars[1]
eG_cst=pars[2]
eH_cst=pars[2]
c_dep =pars[3]
q     =pars[4] 

filein_pre='N42W123_'  # for a single point (42N, 123W) as a test. As for a global 2-degree simulation with run_year=300, the cpu-time is about 6 hours with 8 cores in parallel. 
fileout='outname.nc'  # name of the output file

main.run(pathin,pathout,res,run_year,phi_v,phi_a,eG_cst,eH_cst,c_dep,q,filein_pre,fileout)
