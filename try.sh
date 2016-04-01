#!/bin/bash
#$ -N try
#$ -cwd
#$ -j Y
#$ -V
#$ -m be
#$ -M andrew.haoyu@gmail.com

R CMD BATCH try.R