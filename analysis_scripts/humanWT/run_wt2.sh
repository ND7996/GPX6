echo '##################################################'
echo 'running qgenrelax'
echo '##################################################'
rm -rf $relax
cp correctedGPX6_wtnewcys /$HOME/Onedrive/Thesis/EVB/GPX_EVB/CALCS/WT $relax
$HOME/github/qtools/qscripts-cli/q_genrelax.py ../genrelax.proc \
          --top GPX_wtnewcys.top \
          --rs ../../COMMON/run_qdyn.sh \
          --pdb GPX_wtnew.pdb \
          --fep correctedGPX6_wtnewcys.fep \
          --outdir $relax
cd $relax
echo '##################################################'
echo 'running run_qdyn'
echo '##################################################'
ls
sh run_qdyn.sh
echo $dir
cd $dir
ls
echo '##################################################'
echo 'running qgenfep'
echo '##################################################'
$HOME/github/qtools/qscripts-cli/q_genfeps.py ../genfeps.proc $relax/relax_008.inp relax \
          --rs ../../COMMON/run_qdyn.sh \
           --repeats 1 \
           --frames 51 \
           --fromlambda 0.90 \
           --pdb GPX_wtnewcys.pdb \
           --prefix replicamousesec




#$HOME/github/qtools/qscripts-cli/q_genrelax.py ../genrelax.proc --top GPX6_wt.top --rs ../../COMMON/run_qdyn.sh --pdb GPX6_wt.pdb --fep correctedGPX6_wt.fep --outdir relax



#~/github/qtools/qscripts-cli/q_genfeps.py ../genfeps.proc ~/scratchq6/relax_wt/relax_008.inp relax --rs ../../COMMON/run_qdyn.sh --repeats 10 --frames 51 --fromlambda 0.90 --pdb GPX_wtnew.pdb --fep correctedGPX6_wtnew.fep --prefix replica


