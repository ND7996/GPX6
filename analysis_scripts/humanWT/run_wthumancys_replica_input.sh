#INPUT FOR FEP
rm -rf replicamousecys*
# Qtools is installed on PATH in the Docker image.

echo '##################################################'
echo 'running qgenfep'
echo '##################################################'
q_genfeps.py genfeps.proc ./analysis_scripts/mouseWT/relax/relax_008.inp relax \
          --rs run_qdyn_5.sh \
           --repeats 5 \
           --frames 51 \
           --fromlambda 0.90 \
           --pdb GPX6cys_mouse.pdb \
           --prefix replicamousecys



