#!/bin/bash
test="c"
# Wrong chibstate
if [[ $test = *"a"* ]]; then
./chibgen --outfile "foo.root" --nevents 1000 --seed 1 --chibstate 3
./chibgen --outfile "foo.root" --nevents 1000 --seed 1 --chibstate -3
fi 
# Wrong helicity, first should be fine
if [[ $test = *"b"* ]];then
./chibgen --outfile "foo.root" --nevents 1000 --seed 1 --chibstate 1 --helicity1 0.8 --helicity2 0.7
./chibgen --outfile "foo.root" --nevents 1000 --seed 1 --chibstate 1 --helicity1 1.8
./chibgen --outfile "foo.root" --nevents 1000 --seed 1 --chibstate 2 --helicity1 0.8 --helicity2 0.7
./chibgen --outfile "foo.root" --nevents 1000 --seed 1 --chibstate 2 --helicity1 0.8 --helicity2 0.2
./chibgen --outfile "foo.root" --nevents 1000 --seed 1 --chibstate 2 --helicity1 0.7 --helicity2 -0.2
fi 
# Wrong kinematics
if [[ $test = *"c"* ]]; then
./chibgen --outfile "foo.root" --nevents 1000 --seed 1 --chibstate 2 --helicity1 0.333 --helicity2 0.111 --ptmin -5 --ptmax 50 --absrapmin -1 --absrapmax 1.45 
./chibgen --outfile "foo.root" --nevents 1000 --seed 1 --chibstate 2 --helicity1 0.333 --helicity2 0.111 --ptmin 10 --ptmax 5 --absrapmin 10 --absrapmax 5
 fi 