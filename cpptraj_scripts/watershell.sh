#! /bin/bash                      

. /software/amber20/amber.sh
rm trajin.in rmsf.dat rmsf.lig.dat
echo "reference ./GPX6Hcys.inpcrd" >> trajin.in
for i in `seq 0 20 10000`;do
bzip2 -df GPX6Hcys.md$i.rst.bz2
echo "trajin ./GPX6Hcys.md'$i'.rst" >> trajin.in
done
echo "autoimage
center origin :1-197
solvent byres :1-197
watershell :49 watershell500.dat
go" >> trajin.in

cpptraj -p ./GPX6Hcys.prmtop < trajin.in
