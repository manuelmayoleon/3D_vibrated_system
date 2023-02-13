#!/bin/bash

R=0.5 # R is 0.5 to have sigma = 1
Nx=16
Ny=16
Nz=2
t0=0.0
# tf=1000000
tf=500000
alpha=(0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95)

# alpha=(0.30)

omega=0.001
ns=(1 2 3 4 5 6 7 8 9) #number of repetitions
#ns=(12 13 14 15 16 17 18 19 20) #number of repetitions


prog=$(basename $0)                    # name of this file
progn=$(basename $prog .${prog##*.})   # name of this file without the .sh extension
log=$progn.log                         # log  file
#dat=$progn.dat                         # data file

date                                                        >>  $log
echo "#===================================================" >>  $log
echo "# Compiling 3Dgas.cpp                               " >>  $log
g++ -O3  3Dgas.cpp -o 3Dgas                               >>   $log
echo "# Done compilation 3Dgas.exe was executed           " >>  $log

# var=1.23
# echo "${var:2}" 

mkdir DataSim

# -O3 in g++ turns on all optimizations specified by -O2. see
# https://gcc.gnu.org/onlinedocs/gcc-9.3.0/gcc/Optimize-Options.html#Optimize-Options

# Run for default values:

 echo "# Executing 3Dgas.exe ${n} times                   " >>  $log
# for k in ${alpha[@]};do

# mkdir /DataSim/alpha_p${k}

for k in ${alpha[@]};do

mkdir DataSim/alpha_p${k:2}

for n in ${ns[@]};do
  ./3Dgas <<EOF >> log
  $R $Nx $Ny $Nz
  $t0 $tf
  $k $omega
EOF

echo " # Calculating R${n}                              " >>  $log

mkdir R${n}
mv 3Dgas_* R${n}
mv R${n} DataSim/alpha_p${k:2}

done

done

#  ---------------------------------------------------------------------
#  Copyright by Konstantinos N. Anagnostopoulos (2004-2014)
#  Physics Dept., National Technical University,
#  konstant@mail.ntua.gr, www.physics.ntua.gr/~konstant
#  
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, version 3 of the License.
#  
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#  
#  You should have received a copy of the GNU General Public Liense along
#  with this program.  If not, see <http://www.gnu.org/licenses/>.
#  -----------------------------------------------------------------------

