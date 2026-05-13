rm -rf ./analysis_scripts/humanWT/relax
rm -rf replicahumancys*
rm -rf relax
relax=./analysis_scripts/humanWT/relax
# Qtools is installed on PATH in the Docker image.
dir=`pwd` 
echo '##################################################'
echo 'running qprep5'
echo '##################################################'
qprep5 <prep5.inp
echo '##################################################'
echo 'running makefep'
echo '##################################################'
python makeFEPhumancys.py

#INPUT RELAX

echo '##################################################'
echo 'running qgenrelax'
echo '##################################################'

#mkdir relax
q_genrelax.py genrelax.proc \
          --top GPX6cys_human.top \
          --rs  run_qdyn_5.sh \
          --pdb GPX6cys_human.pdb  \
          --fep corrected.fep\
          --outdir relax
cd relax
cd $relax














