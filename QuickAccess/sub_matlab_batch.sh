#!/bin/bash
#SBATCH -p v6_384
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 48
export PATH=/public1/home/sch2396/software-sch2396/matlab2021a-install/matlab2021a/bin:$PATH
matlab < QuickAccess_Mtd2_PorousInfillOpti.m
