Open in maestro
Click Protein Preparation - Preparation workflow
- **Preprocess** - Cap termini and fill the missing side chains
- **check structures** - add hydrogens if required from EDIT - ADD HYDROGENS
- **3D builder** - Adds carbon, oxygen, charges etc
- **optimize** - this step checks the propka and adds flips whereever required
- export the structure as .mae (maestroGPX6_wt.mae) file from maestro
4. Convert the .mae file into parameter from maestro using ffld_server.

 Command - /opt/schrodinger2022-1/utilities/ffld_server -imae maestroGPX6_wt.mae -version 14 -print_parameters -out_file GPX_PARAM.log

 /opt/schrodinger2022-1/utilities/ffld_server -ipdb h2o2.pdb -print_parameters -version 14 > h2o2_param.ffld11

5. Creating the .prm , .lib , .prm.chk using the q_ffld2q.py

 Required:
  ffld_output  ffld_server output  ,   pdb  PDB structure file WHICH WAS USED TO CREATE THE FFLD_OUTPUT

**Command : q_ffld2q.py  h2o2_param.ffld11 h2o2.pdb**

