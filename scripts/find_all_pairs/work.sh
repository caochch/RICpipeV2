#process mapped results
samtools view -@ 10 -b -S -o patientA_rep1_read1_toGenomeAligned.out.bam patientA_rep1_read1_toGenomeAligned.out.sam 
samtools view -@ 10 -b -S -o patientA_rep1_read2_toGenomeAligned.out.bam patientA_rep1_read2_toGenomeAligned.out.sam 
samtools sort -@ 10 patientA_rep1_read1_toGenomeAligned.out.bam patientA_rep1_read1_toGenomeAligned.out.sort 
samtools sort -@ 10 patientA_rep1_read2_toGenomeAligned.out.bam patientA_rep1_read2_toGenomeAligned.out.sort 

samtools view -@ 10 -h -q 30 -F 256 -b -o patientA_rep1_read1_toGenomeAligned.out.sort.uniq.bam patientA_rep1_read1_toGenomeAligned.out.sort.bam 
samtools view -@ 10 -h -q 30 -F 256 -b -o patientA_rep1_read2_toGenomeAligned.out.sort.uniq.bam patientA_rep1_read2_toGenomeAligned.out.sort.bam
samtools view -h -o patientA_rep1_read1_toGenomeAligned.out.sort.uniq.sam patientA_rep1_read1_toGenomeAligned.out.sort.uniq.bam 
samtools view -h -o patientA_rep1_read2_toGenomeAligned.out.sort.uniq.sam patientA_rep1_read2_toGenomeAligned.out.sort.uniq.bam
perl precess_Chimeric_sam.pl patientA_rep1_read1_toGenomeChimeric.out.sam > patientA_rep1_read1_toGenomeChimeric.out.processed.sam 
perl precess_Chimeric_sam.pl patientA_rep1_read2_toGenomeChimeric.out.sam > patientA_rep1_read2_toGenomeChimeric.out.processed.sam

#find pairs
perl obtain_pairs_from_pair.pl patientA_rep1_read1_toGenomeAligned.out.sort.uniq.sam patientA_rep1_read2_toGenomeAligned.out.sort.uniq.sam
samtools view -@ 10 -b -S -o interaction_from_pair_mapped_reads_1.bam interaction_from_pair_mapped_reads_1.sam
samtools view -@ 10 -b -S -o interaction_from_pair_mapped_reads_2.bam interaction_from_pair_mapped_reads_2.sam
samtools sort -n -@ 10 interaction_from_pair_mapped_reads_1.bam interaction_from_pair_mapped_reads_1.sort
samtools sort -n -@ 10 interaction_from_pair_mapped_reads_2.bam interaction_from_pair_mapped_reads_2.sort
samtools view -h -o interaction_from_pair_mapped_reads_1.sort.sam interaction_from_pair_mapped_reads_1.sort.bam
samtools view -h -o interaction_from_pair_mapped_reads_2.sort.sam interaction_from_pair_mapped_reads_2.sort.bam

perl obtain_pairs_from_gapped_reads.pl patientA_rep1_read1_toGenomeAligned.out.sort.uniq.sam patientA_rep1_read2_toGenomeAligned.out.sort.uniq.sam patientA_rep1_read1_toGenomeChimeric.out.processed.sam patientA_rep1_read2_toGenomeChimeric.out.processed.sam interaction_from_gapped_reads.sam

perl merge_interaction.pl num_of_interactions_from_part.list interaction_from_pair_mapped_reads_1.sort.sam interaction_from_pair_mapped_reads_2.sort.sam interaction_from_gapped_reads.sam patientA_rep1_read1_toGenomeChimeric.out.processed.sam patientA_rep1_read2_toGenomeChimeric.out.processed.sam patientA_rep1_interaction.sam
perl count_link_for_each_kind.pl interaction_from_pair_mapped_reads_1.sort.sam interaction_from_pair_mapped_reads_2.sort.sam num_of_interactions_from_part.list patientA_rep1_read1_toGenomeChimeric.out.processed.sam patientA_rep1_read2_toGenomeChimeric.out.processed.sam > num_of_interactions.list
