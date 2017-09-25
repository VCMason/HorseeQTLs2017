#!/bin/bash
#$ -q all.q
#$ -cwd

module add R/3.2.2;

/home/vmason/bin/bin/eqtlbma_avg_bfs --thread 1 --in "out_eqtlbma_l10abfs_raw.txt.gz" --gwts grid_weights.txt --nsubgrp 4 --dim 15 --cwts config_weights.txt --bestsnp 1 --save bf+post --pi0 0.554662 --post a+b+c+d --bestdim --alldim --out out_eqtlbma_avg_bfs.1.txt.gz >& stdout_eqtlbma_avg_bfs.1.txt;