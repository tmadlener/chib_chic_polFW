#!/bin/bash
#SBATCH -J chib_toymc
#SBATCH -D /scratch/jnecker/toymc
#SBATCH -o /scratch/jnecker/toymc/log/log_%A.out
#/afs/hephy.at/work/j/jnecker/data/chib_results/masterthesis/toymc

# Parameter:
# $1 outfile
# $2 nevents
# $3 chibstate
# $4 R1
# $5 R2
code_base_dir="/afs/hephy.at/user/j/jnecker/code/chib_chic_polFW/toy_mc_generation"

${code_base_dir}/chibgen \
   --outfile "$1" --nevents "$2" --chibstate "$3" \
   --helicity1 "$4" --helicity2 "$5" \
   --ptmin 7.5 --ptmax 60 --absrapmin 0 --absrapmax 1.25
