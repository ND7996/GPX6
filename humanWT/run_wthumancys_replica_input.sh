#INPUT FOR FEP
rm -rf replicamousecys*
source $HOME/github/qtools/qtools_init.sh 

echo '##################################################'
echo 'running qgenfep'
echo '##################################################'
$HOME/github/qtools/qscripts-cli/q_genfeps.py genfeps.proc /home/nayanika/github/PhD_Thesis/EVB/protein_stepwise/GPX6WT/mousecys/1-prep/1-prep/relax/relax_008.inp relax \
          --rs run_qdyn_5.sh \
           --repeats 5 \
           --frames 51 \
           --fromlambda 0.90 \
           --pdb GPX6cys_mouse.pdb \
           --prefix replicamousecys


