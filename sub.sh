#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
source /public1/soft/modules/module.sh
source /public1/home/sc30439/NUCAS_source/new/env.sh

cd /public1/home/sc30439/xiuli20210521/BEPS_PF/realsimulation/US_Ne2_corn
bld/prepare_case.sh --org --expdir pkg_test2003/bepsorg
