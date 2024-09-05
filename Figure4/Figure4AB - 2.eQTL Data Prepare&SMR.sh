reference_file=/COVIDMDD/Data/g1000_eur/g1000_eur
gwas_dir=/COVIDMDD/Data/GWAS
eqtl_dir=/COVIDMDD/Data/Blood_cis_eQTL/

# 1. Mecs, Blood eQTL meta analysis
cd $eqtl_dir
cat file.list
# ./eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense
# ./CAGE-Blood/cage_eqtl_data/CAGE.sparse

smr --besd-flist file.list --mecs --thread-num 5 --out Blood_cis_eQTL

# SMR 
smr --bfile $reference_file --gwas-summary $gwas_dir/smr_input_finn_dep_gwas.txt --beqtl-summary $eQTL_dir/Blood_cis_eQTL --thread-num 20 --out SMR_output_eqtlblood_finndep > SMR_output_eqtlblood_finndep.log 2>&1 &
smr --bfile $reference_file --gwas-summary $gwas_dir/smr_input_finn_scz_gwas.txt --beqtl-summary $eQTL_dir/Blood_cis_eQTL --thread-num 20 --out SMR_output_eqtlblood_finnscz > SMR_output_eqtlblood_finnscz.log 2>&1 &
smr --bfile $reference_file --gwas-summary $gwas_dir/smr_input_pgc_asd_gwas.txt --beqtl-summary $eQTL_dir/Blood_cis_eQTL --thread-num 20 --out SMR_output_eqtlblood_pgcasd > SMR_output_eqtlblood_pgcasd.log 2>&1 &
smr --bfile $reference_file --gwas-summary $gwas_dir/smr_input_pgc_bip_gwas.txt --beqtl-summary $eQTL_dir/Blood_cis_eQTL --thread-num 20 --out SMR_output_eqtlblood_pgcbip > SMR_output_eqtlblood_pgcbip.log 2>&1 &
smr --bfile $reference_file --gwas-summary $gwas_dir/smr_input_pgc_mdd_gwas.txt --beqtl-summary $eQTL_dir/Blood_cis_eQTL --thread-num 20 --out SMR_output_eqtlblood_pgcmdd > SMR_output_eqtlblood_pgcmdd.log 2>&1 &
smr --bfile $reference_file --gwas-summary $gwas_dir/smr_input_pgc_scz_gwas.txt --beqtl-summary $eQTL_dir/Blood_cis_eQTL --thread-num 20 --out SMR_output_eqtlblood_pgcscz > SMR_output_eqtlblood_pgcscz.log 2>&1 &
