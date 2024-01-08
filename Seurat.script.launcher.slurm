#!/bin/bash
#----------------------------------------------------
# Sample Slurm job script
#   for TACC Frontera CLX nodes
#
#   *** Serial Job in Small Queue***
# 
# Last revised: 22 June 2021
#
# Notes:
#
#  -- Copy/edit this script as desired.  Launch by executing
#     "sbatch clx.serial.slurm" on a Frontera login node.
#
#  -- Serial codes run on a single node (upper case N = 1).
#       A serial code ignores the value of lower case n,
#       but slurm needs a plausible value to schedule the job.
#
#  -- Use TACC's launcher utility to run multiple serial 
#       executables at the same time, execute "module load launcher" 
#       followed by "module help launcher".
#----------------------------------------------------

#SBATCH -J jobname           # Job name
#SBATCH -o jobname.o%j       # Name of stdout output file
#SBATCH -e jobname.e%j       # Name of stderr error file
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -p normal
#SBATCH -t 08:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A MCB21061       # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=yan3.liu@utsouthwestern.edu

# Any other commands must follow all #SBATCH directives...

# Launch serial code...
PI="Kim"
Code_DIR="/home1/08270/yanliu74/Seurat_scRNASeq/"
org="human"
res=1.4

prefix="six_samples"
samples="Y270L Y270R Y276L Y276R Y277L Y277R Y278L Y278R Y280L Y280R Y284L Y284R"
In_Dir="/scratch/08270/yanliu74/Harry_Kim/six_samples/"

#Stage_0: QC before filter
Rscript ${Code_DIR}/Stage_0_merging_QC_Before_filter.R ${In_Dir} ${prefix}
wait

#Stage_1: QC after filter
Rscript ${Code_DIR}/Stage_1_merging_QC_After_filter.R ${In_Dir} ${prefix}
wait

#Stage_2
Rscript ${Code_DIR}/Stage_2_phasing_splitsamples.R ${In_Dir} ${prefix} ${org}
wait

#Stage_3
for i in ${samples}; do Rscript ${Code_DIR}/Stage_3_splitsamples_downstream.R ${In_Dir} ${i} ${res} ${org}; done
wait

#Stage_4
Rscript ${Code_DIR}/Stage_4_phasing_integration.R ${In_Dir} ${prefix} ${org}
wait

#Stage_5
Rscript ${Code_DIR}/Stage_5_integration_downstream.R ${In_Dir} ${prefix} ${res} ${In_Dir}/Seurat_input_data/SamplesSheet.csv ${In_Dir}/Seurat_output/data/results/
wait

tar -zcvf six_samples_Res1.4_Seurat_output.tar.gz Seurat_output
