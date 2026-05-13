rm -rf ./analysis_scripts/mouseWT/relax
rm -rf replicamousecys*
rm -rf relax
relax=./analysis_scripts/mouseWT/relax
# Qtools is installed on PATH in the Docker image.
dir=`pwd` 
echo '##################################################'
echo 'running qprep5'
echo '##################################################'
qprep5 <prep5.inp
echo '##################################################'
echo 'running makefep'
echo '##################################################'
python makeFEPmousecys.py
#does not create the file with atoms corresponsing to topolgy and needs to be manually created 
# creates generatedGPX6_wt.fep - the atom number corresponds to q atom numbers in qmap.fep, manually need to be edited to topology number

#INPUT RELAX

echo '##################################################'
echo 'running qgenrelax'
echo '##################################################'

#mkdir relax
q_genrelax.py genrelax.proc \
          --top GPX6cys_mouse.top \
          --rs  run_qdyn_5.sh \
          --pdb GPX6cys_mouse.pdb  \
          --fep corrected.fep\
          --outdir relax
cd relax
cd $relax














