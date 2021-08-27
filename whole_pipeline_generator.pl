#!/usr/bin/perl

my @samples=(
"1.METTL3KO_rep1",
"2.METTL3KO_rep2",
);

my $STAR_index_dir="/mnt/pub/work/caochch/NCBIM37.mm9";
my $genome_ref_fasta="/mnt/pub/work/caochch/genome.NCBIM37.clean.fa";

open(TSH,">run.sh") || die;
foreach my $d (@samples){
	my $prefix=$d;
	$prefix=~s/^\d+\.//;

	#`mkdir ./$d`;
        my $in_read_1="../../6.split_chrM_and_otherReads/$d/z1.second_round_byBWA/read1_torRNA_Unmapped_really.FutherByBwa.fq";
        my $in_read_2="../../6.split_chrM_and_otherReads/$d/z1.second_round_byBWA/read2_torRNA_Unmapped_really.FutherByBwa.fq";
        #`ln -s $in_read_1 ./$d`;
        #`ln -s $in_read_2 ./$d`;

	#print TSH "STAR --runMode alignReads --genomeDir $STAR_index_dir --readFilesIn ./$d/read1_torRNA_Unmapped_really.FutherByBwa.fq --outFileNamePrefix ./$d/$prefix","_read1_toGenome --outReadsUnmapped Fastx --outFilterMultimapNmax 100 --outSAMattributes All --alignIntronMin 1 --scoreGapNoncan -4 --scoreGapATAC -4 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --limitOutSJcollapsed 10000000 --limitIObufferSize 1500000000 --runThreadN 16 --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --alignSJstitchMismatchNmax 5 -1 5 5 --outFilterMatchNminOverLread 0.5 --outFilterScoreMinOverLread 0.5 \n";
	#print TSH "STAR --runMode alignReads --genomeDir $STAR_index_dir --readFilesIn ./$d/read2_torRNA_Unmapped_really.FutherByBwa.fq --outFileNamePrefix ./$d/$prefix","_read2_toGenome --outReadsUnmapped Fastx --outFilterMultimapNmax 100 --outSAMattributes All --alignIntronMin 1 --scoreGapNoncan -4 --scoreGapATAC -4 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --limitOutSJcollapsed 10000000 --limitIObufferSize 1500000000 --runThreadN 16 --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --alignSJstitchMismatchNmax 5 -1 5 5 --outFilterMatchNminOverLread 0.5 --outFilterScoreMinOverLread 0.5\n";

	#print TSH "perl /media/ibm_disk/work/caochch/project_My/RICpipeV2_20210713/7.mapping_and_pairs/select_not_fully_aligned.pl ./$d/read1_torRNA_Unmapped_really.FutherByBwa.fq ./$d/$prefix","_read1_toGenomeAligned.out.sam ./$d/$prefix","_read1_toGenomeChimeric.out.sam > ./$d/$prefix","_read1_toGenomeUnmapped_really.fq &\n";
	#print TSH "perl /media/ibm_disk/work/caochch/project_My/RICpipeV2_20210713/7.mapping_and_pairs/select_not_fully_aligned.pl ./$d/read2_torRNA_Unmapped_really.FutherByBwa.fq ./$d/$prefix","_read2_toGenomeAligned.out.sam ./$d/$prefix","_read2_toGenomeChimeric.out.sam > ./$d/$prefix","_read2_toGenomeUnmapped_really.fq\n";

	#`mkdir ./$d/z1.second_round_byBWA`;
	#print TSH "/media/ibm_disk/work/caochch/project_My/RICpipeV2_20210713/bwa-master/bwa mem -t 16 -k 12 -T 15 -o ./$d/z1.second_round_byBWA/read1_futher_by_bwa.sam $genome_ref_fasta ./$d/$prefix","_read1_toGenomeUnmapped_really.fq\n";
	#print TSH "/media/ibm_disk/work/caochch/project_My/RICpipeV2_20210713/bwa-master/bwa mem -t 16 -k 12 -T 15 -o ./$d/z1.second_round_byBWA/read2_futher_by_bwa.sam $genome_ref_fasta ./$d/$prefix","_read2_toGenomeUnmapped_really.fq\n";

	#print TSH "perl /media/ibm_disk/work/caochch/project_My/RICpipeV2_20210713/7.mapping_and_pairs/collect_chimeric_ligation_from_sam.pl ./$d/z1.second_round_byBWA/read1_futher_by_bwa.sam 1 ./$d/z1.second_round_byBWA > ./$d/z1.second_round_byBWA/out1.read1.chimeric.sam\n";
	#print TSH "perl /media/ibm_disk/work/caochch/project_My/RICpipeV2_20210713/7.mapping_and_pairs/collect_chimeric_ligation_from_sam.pl ./$d/z1.second_round_byBWA/read2_futher_by_bwa.sam 2 ./$d/z1.second_round_byBWA > ./$d/z1.second_round_byBWA/out1.read2.chimeric.sam \n";

	#print TSH "perl /media/ibm_disk/work/caochch/project_My/RICpipeV2_20210713/7.mapping_and_pairs/get_final_unmapped_reads.pl ./$d/$prefix","_read1_toGenomeUnmapped_really.fq ./$d/z1.second_round_byBWA/unmapped.1.readID.list > ./$d/z1.second_round_byBWA/read1_toGenome_Unmapped_really.FutherByBwa.fq &\n";
	#print TSH "perl /media/ibm_disk/work/caochch/project_My/RICpipeV2_20210713/7.mapping_and_pairs/get_final_unmapped_reads.pl ./$d/$prefix","_read2_toGenomeUnmapped_really.fq ./$d/z1.second_round_byBWA/unmapped.2.readID.list > ./$d/z1.second_round_byBWA/read2_toGenome_Unmapped_really.FutherByBwa.fq\n";

	#you can stop if you do not care about the interaction in this given reference
	
	#run find_all_pairs first

	#`mkdir ./$d/z2.update_interaction_sam`;
	#print TSH "perl /media/ibm_disk/work/caochch/project_My/RICpipeV2_20210713/7.mapping_and_pairs/supp_bwaSam_and_update.pl ./$d/find_all_pairs/$prefix","_interaction.sam ./$d/find_all_pairs/num_of_interactions_from_part.list ./$d/find_all_pairs/num_of_interactions.list ./$d/z1.second_round_byBWA/out1.read1.chimeric.sam ./$d/z1.second_round_byBWA/out1.read2.chimeric.sam ./$d/z2.update_interaction_sam > ./$d/z2.update_interaction_sam/mergeBwa.$prefix","_interaction.sam\n";	
	#mkdir ./$d/z3.update_alignPair`;
	#rint TSH "perl /media/ibm_disk/work/caochch/project_My/RICpipeV2_20210713/7.mapping_and_pairs/replace_AlignPair_by_ChimericFragment.pl ./$d/z2.update_interaction_sam/mergeBwa.$prefix","_interaction.sam ./$d/find_all_pairs/$prefix","_read1_toGenomeChimeric.out.processed.sam ./$d/find_all_pairs/$prefix","_read2_toGenomeChimeric.out.processed.sam ./$d/z1.second_round_byBWA/out1.read1.chimeric.sam ./$d/z1.second_round_byBWA/out1.read2.chimeric.sam > ./$d/z3.update_alignPair/update.mergeBwa.$prefix","_interaction.sam &\n";

 

}
