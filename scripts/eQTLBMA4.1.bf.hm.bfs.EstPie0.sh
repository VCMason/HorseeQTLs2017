#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -pe smp 4

module add R/3.2.2;

/home/vmason/bin/bin/eqtlbma_bf --geno list_GenotypeFiles.txt --scoord 82_RAO_HD.2M.sort.trim2eQTLIndivs.ne1e3.eQTLBMA.AllSNPLoc.Unix.tab.gz --exp list_GenExpFiles.txt --gcoord features.NCBI.tx.residuals.inter.Unix.tab.gz --anchor TSS --analys join --error hybrid --cis 1000000 --out out_eqtlbma --thread=4 --bfs all --outss sum_stats.join.tx --gridL grid_phi2_oma2_general.txt --gridS grid_phi2_oma2_with-configs.txt >& stdout_eqtlbma_bf.txt;
/home/vmason/bin/bin/eqtlbma_hm --data "out_eqtlbma_l10abfs_raw.txt.gz" --thread=4 --nsubgrp 4 --dim 15 --ngrid 10 --out out_eqtlbma_hm.txt.gz >& stdout_eqtlbma_hm.txt;
zcat out_eqtlbma_hm.txt.gz | grep "#grid" | cut -f2 > grid_weights.txt;
zcat out_eqtlbma_hm.txt.gz | grep "#config" | awk '{split($1,a,"."); print a[2]"\t"$2}' > config_weights.txt;
/home/vmason/bin/bin/eqtlbma_avg_bfs --in "out_eqtlbma_l10abfs_raw.txt.gz" --thread=4 --gwts grid_weights.txt --nsubgrp 4 --dim 15 --cwts config_weights.txt --save bf --out out_eqtlbma_avg_bfs_genes.txt.gz >& stdout_eqtlbma_avg_bfs_genes.txt;
Rscript /data2/users/pacholewska/vmason/eQTLBMA/dataIN/hybridTSS_MAF0.05_VSD_binary_TagSNPs_residuals_no43_2Mb/EstimatePie0WithEBF02.R;