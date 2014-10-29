#!/bin/bash
#Test LArSoft code with "prodsingle_uboone.fcl".
# echurch@fnal.gov


cp  ${UBOONECODE_DIR}/job/standard_g4_uboone.fcl .
echo "services.user.FileCatalogMetadataExtras.RenameTemplate: '' " >> ./standard_g4_uboone.fcl

lar -c ./standard_g4_uboone.fcl -s ../lar_ci_hitana_prod/hitana_uboone_prod.root -n -1 -T hitana_uboone_g4_hist.root -o hitana_uboone_g4.root 
