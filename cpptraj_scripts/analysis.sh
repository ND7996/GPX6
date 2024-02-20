#! /bin/bash                      

. /software/amber20/amber.sh
rm trajin.in rmsf.dat rmsf.lig.dat
echo "reference ./mousecysgssg.inpcrd" >> trajin.in
for i in `seq 0 20 5000`;do
bzip2 -df mousecysgssg.md$i.rst.bz2
echo "trajin ./mousecysgssg.md'$i'.rst" >> trajin.in
done
echo "autoimage
center origin :1-189
atomicfluct first :1-189@CA out rmsf.dat byres
go" >> trajin.in

cpptraj -p ./mousecysgssg.prmtop < trajin.in


#! /bin/bash                      

. /software/amber20/amber.sh
rm trajin.in rmsf.dat rmsf.lig.dat
echo "reference ./mousecysgssg.inpcrd" >> trajin.in
for i in `seq 0 20 5000`;do
bzip2 -df mousecysgssg.md$i.rst.bz2
echo "trajin ./mousecysgssg.md'$i'.rst" >> trajin.in
done
echo "autoimage
center origin :1-189
hbond donormask :1-189 acceptormask :49 avgout hbond2.dat dist 3
go" >> trajin.in

cpptraj -p ./mousecysgssg.prmtop < trajin.in



#! /bin/bash                      

. /software/amber20/amber.sh
rm trajin.in rmsf.dat rmsf.lig.dat
echo "reference ./mousecysgssg.inpcrd" >> trajin.in
for i in `seq 0 20 5000`;do
bzip2 -df mousecysgssg.md$i.rst.bz2
echo "trajin ./mousecysgssg.md'$i'.rst" >> trajin.in
done
echo "autoimage
center origin :1-189
distance dist1 :49@SG :83@CB out distance.dat
go" >> trajin.in

cpptraj -p ./mousecysgssg.prmtop < trajin.in




#! /bin/bash                      

. /software/amber20/amber.sh
rm trajin.in rmsf.dat rmsf.lig.dat
echo "reference ./mousecysgssg.inpcrd" >> trajin.in
for i in `seq 0 20 5000`;do
bzip2 -df mousecysgssg.md$i.rst.bz2
echo "trajin ./mousecysgssg.md'$i'.rst" >> trajin.in
done
echo "autoimage
center origin :1-189
solvent byres :1-189
watershell :49 watershell.dat
go" >> trajin.in



#! /bin/bash                      

. /software/amber20/amber.sh
rm trajin.in rmsf.dat rmsf.lig.dat
echo "reference ./mousecysgssg.inpcrd" >> trajin.in
for i in `seq 0 20 5000`;do
bzip2 -df mousecysgssg.md$i.rst.bz2
echo "trajin ./mousecysgssg.md'$i'.rst" >> trajin.in
done
echo "autoimage
center origin :1-189
surf :1-189@CA out sasa.dat	
go" >> trajin.in

cpptraj -p ./mousecysgssg.prmtop < trajin.in


#! /bin/bash                      

. /software/amber20/amber.sh
rm trajin.in rmsf.dat rmsf.lig.dat
echo "reference ./mousecysgssg.inpcrd" >> trajin.in
for i in `seq 0 20 5000`;do
bzip2 -df mousecysgssg.md$i.rst.bz2
echo "trajin ./mousecysgssg.md'$i'.rst" >> trajin.in
done
echo "autoimage
center origin :1-189
secstruct :1-189 out sumout secstr.dat totalout secstr.dat	
go" >> trajin.in

cpptraj -p ./mousecysgssg.prmtop < trajin.in

#! /bin/bash                      

. /software/amber20/amber.sh
rm trajin.in rmsf.dat rmsf.lig.dat
echo "reference ./mousecysgssg.inpcrd" >> trajin.in
for i in `seq 0 20 5000`;do
bzip2 -df mousecysgssg.md$i.rst.bz2
echo "trajin ./mousecysgssg.md'$i'.rst" >> trajin.in
done
echo "autoimage
center origin :1-189
radgyr R1 :1-189@CA out radiusgyr.dat
go" >> trajin.in

cpptraj -p ./mousecysgssg.prmtop < trajin.in


