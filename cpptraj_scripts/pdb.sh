#! /bin/bash

. /software/amber20/amber.sh
rm trajin.in rmsd.dat rmsd.lig.dat
echo "reference ./GPX6Hcys.inpcrd" >> trajin.in
for i in `seq 0 20 5000`;do
bzip2 -df GPX6Hcys.md$i.rst.bz2
echo "trajin ./GPX6Hcys.md'$i'.rst" >> trajin.in
done
echo "autoimage
center origin :1-197
rms first :1-197@CA out rmsd250.dat mass
solvent byres :WAT
closestwaters 200 :46
#strip :WAT,K+,Cl-,Na+
trajout MD250.pdb PDB
go" >> trajin.in

cpptraj -p ./GPX6Hcys.prmtop < trajin.in
