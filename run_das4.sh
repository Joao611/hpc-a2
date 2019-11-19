#!/bin/bash
sbatch --nodes=2 main.job /var/scratch/hpcl1905/mouse_gene.mtx
squeue