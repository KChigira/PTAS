#! /bin/sh

###For developmant. Please ignore.
#pip install -e .
#conda activate mkdesigner
#cd test

#SNP
mkprimer -r test_ref.fasta \
         -V test_mkvcf/vcf_2nd/Merged_filtered_variants.vcf \
         -n1 lineA -n2 lineB \
         -p test_mkprimer_snp -t SNP \
         --mindep 5 --maxdep 120 --mismatch_allowed 5 --cpu 6 \
         --min_prodlen 150 --opt_prodlen 180 --max_prodlen 280 \
         --search_span 140 --primer_min_size 24 \
         --primer_opt_size 24 --primer_max_size 24

#INDEL
mkprimer -r test_ref.fasta \
         -V test_mkvcf/vcf_2nd/Merged_filtered_variants.vcf \
         -n1 lineA -n2 lineB \
         -p test_mkprimer_indel -t INDEL \
         --mindep 5 --maxdep 120 --mismatch_allowed 5 --cpu 4 \
         --search_span 250 --primer_min_size 18 \
         --primer_opt_size 20 --primer_max_size 26