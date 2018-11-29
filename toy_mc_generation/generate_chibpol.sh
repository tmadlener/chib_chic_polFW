#!/bin/bash
#SBATCH -J chib_toymc
#SBATCH -D /afs/hephy.at/work/j/jnecker/data/chib_results/masterthesis/toymc
#SBATCH -o /afs/hephy.at/work/j/jnecker/data/chib_results/masterthesis/toymc/log/log_%A.out

# Parameter:
# $1 outfile
# $2 nevents
# $3 chibstate
# $4 R1
# $5 R2
# $6 1=do NOT apply selection, default: selection is applied

apply_selection="true"
if [[ $# -gt 5 && $6 -eq 1 ]]; then apply_selection="false"; fi
echo $apply_selection

code_base_dir="/afs/hephy.at/user/j/jnecker/code/chib_chic_polFW/toy_mc_generation"

${code_base_dir}/chibgen \
   --outfile "$1" --nevents "$2" --chibstate "$3" \
   --helicity1 "$4" --helicity2 "$5" \
   --ptmin 0 --ptmax 60 --absrapmin 0 --absrapmax 2 --applyselection $apply_selection
 