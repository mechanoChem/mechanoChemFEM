#!/bin/bash

thisroot=$PWD

# 3
mesh_list=(300)
# 5
# M_list=(1e-5 1e-4 1e-3 1e-2)
M_list=(1e-2 1e-1 1e-0)

# 5
# kappa_list=(1e-5 1e-4 1e-3 1e-2)
kappa_list=(1e-3 1e-2 1e-1 1e-0 1e1)

# 5
#flux_list=(-0.0003 -0.001 -0.003 -0.01 -0.03)
flux_list=(-0.03 -0.05 -0.1 -0.2)

# 5
omega_list=(1e-0 1e-1 1e-2 1e-3)

# total 3*6*6*10=1080
# 10GB -> 320 jobs

for mesh0 in "${mesh_list[@]}"
do
  for M0 in "${M_list[@]}" 
  do 
    for kappa0 in "${kappa_list[@]}" 
    do
        for flux0 in "${flux_list[@]}" 
        do
            for omega0 in "${omega_list[@]}" 
            do 
                cd $thisroot
                newfolder=m$mesh0+M$M0+k$kappa0+o$omega0+f$flux0
                echo $newfolder
                mkdir -p $newfolder
                cd $thisroot/$newfolder
                cp ../../template.prm parameters.prm
                cp ../job.script tmp.script
                sed -i "s/mesh_list/$mesh0/g" parameters.prm
                sed -i "s/M_list/$M0/g" parameters.prm
                sed -i "s/kappa_list/$kappa0/g" parameters.prm
                sed -i "s/flux_list/$flux0/g" parameters.prm
                sed -i "s/omega_list/$omega0/g" parameters.prm
                sbatch tmp.script
                sleep 1
#                exit
            done 
        done
    done 
  done
done
