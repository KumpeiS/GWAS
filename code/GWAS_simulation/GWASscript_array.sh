#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=70G
#$ -t 1-1000

module use /usr/local/package/modulefiles/
module load R

Rscript GWASscript_array.R ${SGE_TASK_ID}
