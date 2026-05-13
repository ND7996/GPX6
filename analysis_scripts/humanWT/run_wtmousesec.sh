rm -rf replicamousesec*
rm -rf relax
relax=./analysis_scripts/mouseWT/mousesec
# Qtools is installed on PATH in the Docker image.
dir=`pwd` 
echo '##################################################'
echo 'running qprep5'
echo '##################################################'
qprep5 <prepgpx6mousesec.inp
echo '##################################################'
echo 'running makefep'
echo '##################################################'
python makeFEP.py
#does not create the file with atoms corresponsing to topolgy and needs to be manually created 
# creates generatedGPX6_wt.fep - the atom number corresponds to q atom numbers in qmap.fep, manually need to be edited to topology number

#INPUT RELAX

echo '##################################################'
echo 'running qgenrelax'
echo '##################################################'

#mkdir relax
#cp ./analysis_scripts/mouseWT/relax
q_genrelax.py genrelax.proc \
          --top GPX6sec_mouse.top \
          --rs  run_qdyn.sh \
          --pdb GPX6sec_mouse.pdb  \
          --fep GPX6_wtmousesecnew.fep \
          --outdir relax
cd relax
cd $relax















