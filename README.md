# q2-pan-classifier
 # Enterovirus Classification Tutorial

## Purpose

Using the Qiime2 framework to develop a pipeline to classify pan-viral samples against a database of known viruses. The database will be generated by the user. 



## Pipeline Steps Overview

1. Create Custom Classifier
   1. Following steps found on Qiime2 website: https://docs.qiime2.org/2021.4/tutorials/feature-classifier/
      1. import reference sequences
      2. Extract reference reads
      3. Train Classifier
      4. Test Classifier
2. Classifying Samples (Modified Procedure from Qiime2 Website)
   1. Preparing Sequences
      1. Trim sequences
      2. Import sequences
      3. Check quality of reads
   2. Clustering Sequences (Dada2)
   3. Classify Clusters
   4. Generate Visualizations

Below is a tutorial classifying enteroviruses using pan-entero primers. 



# Entero Classification Tutorial



Note: This tutorial expects that you have Qiime2 already installed. If you do not, please follow the directions on their website, please: 

https://docs.qiime2.org/2022.8/install/



## Clone git repository

```
git clone https://github.com/CRideTGen/q2-pan-classifier

cd test_analysis
```


## Primers Used
| Primer Name | Sequence             |                 
|-------------| ---------------------|
| forward     | CAAGCACTTCTGTTTCCCCGG|
| reverse     | ATTGTCACCATAAGCAGCCA |



## Part 1: Create Custom Classifier

1. ### Retrieve Reference Sequences

```bash
   qiime rescript get-ncbi-data \
       --p-query "Picornaviridae[ORGANISM] AND 5000:10000000[SLEN]" \
       --output-dir references
```

     

2. ### Extract Reference Reads

```bash
  qiime feature-classifier extract-reads \
      --i-sequences references/sequences.qza \
      --p-f-primer CAAGCACTTCTGTTTCCCCGG \
      --p-r-primer ATTGTCACCATAAGCAGCCA \
      --p-min-length 0 --p-max-length 500 \
      --o-reads references/sequences_trimmed.qza
  
  qiime rescript dereplicate \
      --i-sequences references/sequences_trimmed.qza \
      --i-taxa references/taxonomy.qza \
      --p-derep-prefix \
      --o-dereplicated-sequences references/sequences_dereplicated.qza \
      --o-dereplicated-taxa references/taxonomy_dereplicated.qza
```

      

3. ### Train Classifier

```bash
  qiime rescript evaluate-fit-classifier \
      --i-sequences references/sequences_dereplicated.qza \
      --i-taxonomy references/taxonomy_dereplicated.qza \
      --output-dir classifier
```


## Part 2: Classifying Samples

### objective:

We will be using the classifier generated in part1 to classify the 5 example samples found in the folder **example_sequences**



## Preparing Sequences


1. ### Importing Sequences into Qiime2

   1. #### Create manifest file

   This file is used by qiime2 to import the sequences into a qiime2 specific .qza file. Documentation on this import can be found at: [link]( https://docs.qiime2.org/2020.6/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq)
   scroll down to the section titled *“Fastq manifest” formats*

   The make_manifest_file.sh bash script will be used to generate your manifest:

   ```bash
   bash ./scripts/make_manifest_file.sh -d example_sequences -o entero_manifest
   
   #check the generated manifest file
   
   cat entero_manifest | column -t
   ```

   

2. #### 	qiime2 import command

   ```bash
   module load qiime2
    
   qiime tools import \
     --type 'SampleData[PairedEndSequencesWithQuality]' \
     --input-path entero_manifest \
     --output-path paired_end_demux_entero.qza \
     --input-format PairedEndFastqManifestPhred33V2
   ```

   1. #### Notes on inputs: 

      - --type: tells qiime2 what type of sequence data is being imported, e.g., single or paired end

      - --input-path: path to manifest file

      - --output-path: name of output file

      - --input-format tells qiime2 what format the sequencing data is. See link above for more information
   

      
3. ####  Trim Sequences using q2_cutadapt

         cutadapt will trim off the primers by matching the supplied primer sequences to the first 30bp of the sequencing reads.
   ```bash
    qiime cutadapt trim-paired \
      --i-demultiplexed-sequences paired_end_demux_entero.qza \
      --p-front-f CAAGCACTTCTGTTTCCCCGG \
      --p-front-r ATTGTCACCATAAGCAGCCA \
      --o-trimmed-sequences paired_end_demux_entero_trimmed.qza
   ```

4. ### Visualization of sequence quality

   This step produces plots of the average quality per base for the forward and reverse reads. This information is used to trim off the bases of low quality within the reads.

   To produce the qiime2 visualization file, .qzv, use the following command:

   ```bash
   #make directory to store visualizations
   
   mkdir visual_output
   
   #creating visualization   
   qiime demux summarize \
        --i-data paired_end_demux_entero_trimmed.qza \
        --o-visualization visual_output/demux_entero.qzv
   ```

            To visualize the .qzv file you have to visit https://view.qiime2.org/

   

5. ## Clustering sequences using DADA2

   1. ​	This section is where we cluster the reads into 100% OTUs using the algorithm implemented by DADA2

      Using the demux_AA.qzv we will decide where to trim and truncate our forward and reverse reads. 

      ```
      qiime dada2 denoise-paired \
           --i-demultiplexed-seqs paired_end_demux_entero_trimmed.qza \
           --p-trim-left-f 0 \
           --p-trim-left-r 0 \
           --p-trunc-len-f 200 \
           --p-trunc-len-r 200 \
           --o-representative-sequences rep_seqs_dada2_entero.qza \
           --o-table table_dada2_entero.qza \
           --o-denoising-stats stats_dada2_entero.qza
      
      ```

      #### Notes on inputs:

   - --p-trim-left-f Trims all of the bases upto the base given for the forward reads. e.g., --p-trim-left-f 25 will trim base position 1-25 off all the reads.

   *  --p-trim-left-r same as above but for reverse reads.
   *  --p-trunc-len-f Truncates all of the bases after the base given for the forward reads. e.g., --p-trunc-len-f 200 will remove bases 201-end of read.
   *  --p-trunc-len-r same as above but for reverse reads.

6. ### Generate Metadata file using Keemei  

   1. See instructions on the qiime2 website: https://docs.qiime2.org/2021.4/tutorials/metadata/
   2. save file as ***ENT_meta.tsv***

   

7. ### Visualizing feature table 

   1. ```bash
      qiime feature-table summarize \
        --i-table table_dada2_entero.qza \
        --o-visualization visual_output/table_entero.qzv \
        --m-sample-metadata-file ENT_meta.tsv
      
      qiime feature-table tabulate-seqs \
        --i-data rep_seqs_dada2_entero.qza \
        --o-visualization visual_output/rep_seqs_entero.qzv
        
      ```



## Classification

1. Using a naïve Bayes classifier, the features are classified against all publicly available virus genomes.

   1. ```bash
      qiime feature-classifier classify-sklearn   \
      --i-classifier ../classifier/classifier.qza  \
       --i-reads rep_seqs_dada2_entero.qza   \
       --o-classification taxonomy_entero.qza
      ```

      

2. ## Visualizing Results

   1. #### Create Taxa Bar plot

   2. ```bash
       qiime taxa barplot \
       --i-table table_dada2_entero.qza \
       --i-taxonomy taxonomy_entero.qza \
       --m-metadata-file ENT_meta.tsv \
       --o-visualization visual_ouput/taxa_bar_plots_entero.qzv
      
      ```

3. ## Create Combined Table 

   This section we will produce a table that contains the feature id, sequences, and feature counts per sample.

   1. we have to reorient the feature table so the columns and rows match the taxonomy_AA.qza and paired-rep-seqs-dada2.qza. This is accomplished using the 
     feature-table transpose command:

   ```bash
   qiime feature-table transpose \
     --i-table table_dada2_entero.qza \
     --o-transposed-feature-table transposed_table_dada2_entero.qza \
   
   ```

   Now we generate the table:

   ```bash
   qiime metadata tabulate 
     --m-input-file taxonomy_entero.qza \
     --m-input-file rep_seqs_dada2_entero.qza \
     --m-input-file transposed_table_dada2_entero.qza \
     --o-visualization visual_output/feat-tax-rep.qzv
   ```

   
