#! /bin/sh

#pip install -e .
#conda install bcftools
#conda install bgzip
#conda install r-base -c conda-forge

mkselect -i test_mkprimer_snp/ref/test_ref.fasta.fai \
         -V test_mkprimer_snp/vcf/Merged_filtered_variants_selected_primer_added.vcf \
         -n 10000
mkselect -i test_mkprimer_snp/ref/test_ref.fasta.fai \
         -V test_mkprimer_snp/vcf/Merged_filtered_variants_selected_primer_added.vcf \
         -n 20
mkselect -i test_mkprimer_snp/ref/test_ref.fasta.fai \
         -V test_mkprimer_snp/vcf/Merged_filtered_variants_selected_primer_added.vcf \
         -n 20 -d density.tsv

