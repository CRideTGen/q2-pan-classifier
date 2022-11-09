#!/bin/env bash

#export TMPDIR=/scratch/sps266/tmp
#module load anaconda3
#conda activate qiime2-2022.8

query_string=$1
forward_primer=$2
reverse_primer=$3
min_length=$4
max_length=$5

echo "Query String: $query_string"
echo "Forward Primer: $forward_primer"
echo "Reverse Primer: $reverse_primer"
echo "Min Length: $min_length"
echo "Max Length: $max_length"

qiime

#qiime rescript get-ncbi-data \
#       --p-query "$query_string" \
#       --output-dir references

