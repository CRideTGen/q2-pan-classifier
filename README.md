# q2-pan-classifier
 # Flavivivirus Analysis of potential co-infected WNV and SLEV
 
 ## Making classifier
 pan-flavivirus analysis, here are the primers: PFlav-fAAR (TACAACATGATGGGAAAGAGAGAGAARAA from 9040 to 9068 of AF196835) and PFlav-rKR (GTGTCCCAKCCRGCTGTGTCATC from positions 9305 to 9283 of AF196835
 
 ```{bash}
 
 qiime rescript get-ncbi-data --p-query "[Flaviviridae[ORGANISM] AND 5000:10000000[SLEN]" --output-dir refrences
 qiime feature-classifier extract-reads --i-sequences refrences/sequences.qza --p-f-primer TACAACATGATGGGAAAGAGAGAGAARAA --p-r-primer GTGTCCCAKCCRGCTGTGTCATC --p-min-length 200 --p-max-length 300 --o-reads refrences/sequences_trimmed.qza
 qiime rescript dereplicate --i-sequences sequences_trimmed.qza --i-taxa taxonomy.qza --p-derep-prefix --o-dereplicated-sequences sequences_dereplicated.qza --o-dereplicated-taxa taxonomy_dereplicated.qza
 
 qiime rescript evaluate-fit-classifier --i-sequences refrences/sequences_dereplicated.qza --i-taxonomy refrences/taxonomy_dereplicated.qza --output-dir classifier
 ```
 
 Untrimmed references count: 25086
 trimmed references count: 17241
 dereplicated references: 3790
 
 
 ## Prepping Sequences
 
 directory location: /TGenNextGen/TGN-MiSeq1248/VECTR-PanFlavi/
 ```{bash}
 bash /scratch/cridenour/Crystal/Adenovirus/Analysis/HAdV41_classification_16May2022/scripts/make_manifest_file.sh -d /TGenNextGen/TGN-MiSeq1248/VECTR-PanFlavi/ -o manifest-file
 qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest-file --output-path paired-end-demux.qza --input-format PairedEndFastqManifestPhred33V2
 qiime cutadapt trim-paired --i-demultiplexed-sequences paired-end-demux.qza --p-front-f TACAACATGATGGGAAAGAGAGAGAARAA --p-front-r GTGTCCCAKCCRGCTGTGTCATC --o-trimmed-sequences paired-end-demux-trimmed.qza
 qiime dada2 denoise-paired --i-demultiplexed-seqs paired-end-demux-trimmed.qza --p-trunc-len-f 0 --p-trunc-len-r 0 --output-dir dada2_results
 qiime feature-classifier classify-sklearn --i-reads dada2_results/representative_sequences.qza --i-classifier classifier/classifier.qza --o-classification classification.qza

 
 ```
 