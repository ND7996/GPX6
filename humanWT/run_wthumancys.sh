rm -rf /home/nayanika/github/PhD_Thesis/EVB/EVB-Q6/openMM-OPLS/protein/GPX6WT/humancys/relax/correctedGPX6_humancys.fep/relax
rm -rf replicahumancys*
rm -rf relax
relax=/home/nayanika/github/PhD_Thesis/EVB/EVB-Q6/openMM-OPLS/protein/GPX6WT/humancys/relax/correctedGPX6_humancys.fep/relax
source $HOME/github/qtools/qtools_init.sh 
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
$HOME/github/qtools/qscripts-cli/q_genrelax.py genrelax.proc \
          --top GPX6cys_human.top \
          --rs  run_qdyn_5.sh \
          --pdb GPX6cys_human.pdb  \
          --fep corrected.fep\
          --outdir relax
cd relax
cd $relax













