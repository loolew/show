#!/usr/bin/perl
# by liulif

use strict;
use warnings;
use FindBin;
use Time::Local;
use Parallel::ForkManager;

#====================================================================================================================
# 套餐
my ($Pid);
# 0.Samples_run_in_parallel
my ($Samples_run_in_parallel);
# 1.Reads type
my ($ReadtypeOptions);
# 2.获得参数样本信息
my (@sample,$Samples,$KeyName,$Fastq1,$Fastq1_file,$Fastq2,$Fastq2_file,$Fastq,$Fastq_file,%KEYNAMES);
# 3.QC
my ($Adapter_3,$Adapter_5,$FilterType,$Trim_para_1,$Trim_para_2);
# 4.mapping
my ($Bowtie2Options,$Ref_gene);
# 5.depth
my ($Samtools_depth);
# 6.depth data
my ($Mutalyzer,$D_Bed,$D_Bed_format,$GFF,$IEM_gene_phe);
# 7.Gather_info_each_sample
my ($Sampleinfo);
# 8.samtools mpileup
my ($Mpileup,$Genome);
# 9.annovar
my ($S_annovar,$S_humandbdir,$Xref);
# 10.Local_freq
my ($Local_freq);
# 11.format1
my ($R_HGMD,$Mutalyzer_All,$EXAC,$gnomad,$PP5,$BP6,$Model,$GeneDisease,$Phe,$OMIM,$Clinvar,$GHR,$Hgmd);
#12.format1
my ($G_Trans_gb,$Gene_info_cds,$Verbose,$local_pd,$GeneList);
#====================================================================================================================
my $parameter=$ARGV[0];
get_params($parameter);

#====================================================================================================================
my $test_dir="$FindBin::Bin";
my $test_bin="$test_dir/bin";
my $version=$1 if($test_dir=~/.*\_(.+?)$/);

require("$test_bin/config.pl");
our ($samtools,$annovar);

#==================================================
#定义文件目录
my $output="Output";
my $Input="Input";

#========================================
# 流程开始时间，并创建开始时间的LOG文件
my $run_start_time;
$run_start_time = TIME();
my $run_cmd_log="run_cmd\_$run_start_time.log";
open(LOG,">$run_cmd_log") or die $!;


#======================================得到样本列表=============
foreach my $key (sort keys %KEYNAMES) {
	push @sample,$key;
	$Samples=join",",@sample;
}
#======================================初始化参数=============
CHECKING_PARA_SETS();
#======================================创建文件夹=============
Create_DIR($Samples);
#========================================质控=================
Trimmomatic_QC($Samples);
#========================================比对=================
Bowtie_2($Samples,$Bowtie2Options,$Ref_gene);
#=====================================文库大小信息============
Distribution_of_Library($Pid);
#=====================================样本信息统计============
Sample_info_stat($Samples);
#Gather the information of each sample
Gather_info_each_sample($Samples,$Sampleinfo);
#=====================================样本信息统计============
CALL_SNP($Samples);
#=====================================合并信息================
Merge_sample_information($Samples);
#=====================================合并信息================
Summary($Samples);

#==============================复制筛查点和删除fastq===========
#复制筛选出的点到本地库内
#Local_freq = /work1/yx/TargetSeq/hg19_multianno_IEMS
#`cp ./Output/SNP/*.hg19_multianno.txt /work1/liulifeng/PerlScript_local_freq/hg19_multianno_$Pid`;	
`cp ./Output/SNP/*.hg19_multianno.txt $Local_freq`;	

#删除fq过程文件
	#`rm ./Input/Unpaired/*.*`;
	#`rm ./Input/*.*`;
	#`rm ./Output/Count/*/*.depth.txt`;
	
close LOG;

`chmod 600 $run_cmd_log`;


#============================================================================
#Summary
sub Summary{
	my ($Samples) = shift @_;
	my @samples = split (/,/,$Samples);
	my $pici = $1 if ($Sampleinfo =~ /(.+?)_SampleInfo.txt/);
	my $Summary_pici = $pici."_filter_multianno";
	mkdir("$output/Summary") unless (-e "$output/Summary");
	mkdir("$output/Summary/$Summary_pici") unless (-e "$output/Summary/$Summary_pici");
	foreach my $sample(@samples){
		my $summary = "$test_bin/8.1summary.pl";
		my $in_file = "$output/SNP/All_sample_sorted/${sample}_sorted_all.txt";
		my $job = "Get a summary file of $sample";
		my $cmd = "perl $summary $sample $Summary_pici";
		run_cmd($cmd,$job);
	
	}
	`cp $output/Count/${pici}_All_sample_state.txt $output/Summary`;
}


#============================================================================
#Sort filter
sub Merge_sample_information{
	my ($Samples) = shift @_;
	my @samples = split (/,/,$Samples);
	#Judging the related genes according to the clinical complaint
	#根据临床主诉，排序
	if (-e $Sampleinfo){
		my $sort_filter = "$test_bin/7.1sort_filter_llf.pl";
		my $in_file = "$output/SNP";
		my $job = "Judging the related genes according to the clinical complaint files";
		my $cmd = "perl $sort_filter $in_file";
		run_cmd($cmd,$job);
	}
	#mkdir("$output/Combined_result") unless (-e "$output/Combined_result");
	foreach my $sample(@samples){
		my $combine_txtfile2excle = "$test_bin/7.2combine_txtfile2excle-sheet.pl";
		my $combine_txt2xlsx_file = "$output/Combined_result/$sample".".combined.xls";
		#所有样本统计信息
		my $Gather_info_each_sample_dir = "$output/Count/";
		my $pici = $1 if ($Sampleinfo =~ /(.+?)_SampleInfo.txt/);
		my $Gather_info_each_sample_file = $Gather_info_each_sample_dir.$pici."_All_sample_state.txt";
		#疾病相关是位点文件
		my $sort_filter = "$output/SNP/All_sample_sorted/$sample"."_sorted_filter.txt";
		#新格式所有样本
		my $sort_all = "$output/SNP/All_sample_sorted/$sample"."_sorted_all.txt";
		#过滤后的文件
		my $filter_file = "$output/SNP/$sample".".filter.txt";
		#区间覆盖数文件
		my $covered_block_file = "$output/Low_depth/$sample.low-coverage-block.symptom.txt"; 
		
		my $job1 = "Merge all information to excel of $sample";
		my $cmd1 = "perl $combine_txtfile2excle $Gather_info_each_sample_file $sort_filter $sort_all $filter_file $covered_block_file $combine_txt2xlsx_file";
		run_cmd($cmd1,$job1);
	}
	
}

#============================================================================
#CALL SNP

sub CALL_SNP{
	my ($Samples) = shift @_;
	my @samples = split (/,/,$Samples);
	my $sample_num = $#samples+1;
	my $parallel = ($Samples_run_in_parallel >= $sample_num) ? $sample_num : $Samples_run_in_parallel;
	
	my $pm = Parallel::ForkManager->new($parallel);
	foreach my $sample(@samples){
		my $pid = $pm->start and next;
		#call BCF
		my $call_bcf = "$test_bin/6.1call_bcf.pl";
		my $Mpileup_tmp = $Mpileup;
		$Mpileup_tmp =~ s/\s/,/g;
		my $bam_sorted_dup = "$output/Bowtie2/$sample/$sample.sort.dup.bam";
		my $bcf_file = "$output/SNP/$sample.var.bcf";
		my $job = "Get the BCF file of $sample";
		my $cmd = "perl $call_bcf $Mpileup_tmp $Genome $D_Bed $bam_sorted_dup $bcf_file";
		run_cmd($cmd,$job,$bcf_file);
		#call VCF
		my $call_vcf = "$test_bin/6.2call_vcf.pl";
		my $vcf_file = "$output/SNP/$sample.vcf";
		my $job1 = "Get the BCF file of $sample";
		my $cmd1 = "perl $call_vcf $bcf_file $vcf_file";
		run_cmd($cmd1,$job1,$vcf_file);
		#annovar
		my $Annovar = "$test_bin/6.3annovar.pl";
		my $annovar_qz = "$output/SNP/$sample";#前缀
		my $annovar_file = "$output/SNP/$sample.hg19_multianno.txt";
		my $job2 = "Variation annotation of $sample";
		my $cmd2 = "perl $Annovar $vcf_file $S_humandbdir $annovar_qz $Xref";
		run_cmd($cmd2,$job2,$annovar_file);
		#run_cmd($cmd2,$job2);
		
		#Local_freq and Filter
		my $Filter = "$test_bin/6.4filter_snp.pl";
		my $Filter_file = "$output/SNP/$sample.filter.txt";
		my $job3 = "Get filter and local frequency files of $sample";
		my $cmd3 = "perl $Filter $Local_freq $annovar_file $Filter_file";
		run_cmd($cmd3,$job3,$Filter_file);
		
		#Get the common format 1 大致格式  后面修改格式，从这里开始
		my $format1 = "$test_bin/6.5_format_llf.pl";
		my $format1_file = "$output/SNP/$sample.filter.format1.temp.txt"; #2018030057-IEMS.filter.format1.temp.txt
		my $job4 = "Get the common format.temp files of $sample";
		#my ($R_HGMD,$Mutalyzer_All,$EXAC,$gnomad,$PP5,$BP6,$Model,$GeneDisease,$Phe,$OMIM,$Clinvar,$GHR,$Hgmd);
		my $cmd4 = "perl $format1 $R_HGMD $Filter_file $Mutalyzer_All $EXAC $gnomad $PP5 $BP6 $Model $GeneDisease $Phe $OMIM $Clinvar $GHR $Hgmd";
		run_cmd($cmd4,$job4,$format1_file);
		
		#20180626增加插入缺失重复突变——frameshift
		#################################################
		my $frameshift = "$test_bin/6.6_format_frameshift_0712.pl";
		my $frameshift_file = "$output/SNP/$sample.format.frameshift.txt"; #2018030057-IEMS.format.frameshift.txt
		my $job5 = "Deletion, deletion, repeat mutation corresponding to amino acid change of $sample";
		#my ($R_HGMD,$Mutalyzer_All,$EXAC,$gnomad,$PP5,$BP6,$Model,$GeneDisease,$Phe,$OMIM,$Clinvar,$GHR,$Hgmd);
		my $cmd5 = "perl $frameshift $G_Trans_gb $Gene_info_cds $format1_file $frameshift_file";
		run_cmd($cmd5,$job5);
		
		
		
		#Get the common format 2 增加内含子的位置信息，和是否与疾病相关的基因
		my $format2 = "$test_bin/6.7_format_verbose.pl";
		my $format_file = "$output/SNP/$sample.format.txt";
		my $format_filter_file = "$output/SNP/$sample.format.filter.txt";
		my $job6 = "Get the common format files of $sample";
		my $cmd6 = "perl $format2 $Verbose $GeneList $local_pd $IEM_gene_phe $Sampleinfo $frameshift";
		run_cmd($cmd6,$job6);
		
		$pm->finish;
	}
	$pm->wait_all_children;#等待所有子进程
	
}


#============================================================================
#Gather the information of each sample
sub Gather_info_each_sample{
	my ($Samples,$Sampleinfo) = @_;
	my $Gather_info_each_sample = "$test_bin/5.1Gather_info_each_sample.pl";
	my $Gather_info_each_sample_dir = "$output/Count/";
	
	my $pici = $1 if ($Sampleinfo =~ /(.+?)_SampleInfo.txt/s);
	my $Gather_info_each_sample_file = $Gather_info_each_sample_dir.$pici."_All_sample_state.txt";
	
	my $job = "Gather the information of each sample";
	my $cmd = "perl $Gather_info_each_sample $Pid $Samples $Sampleinfo $Gather_info_each_sample_dir";
	run_cmd($cmd,$job);
	
	#Gene coverage of all samples
	my $Gene_cov_all_samples = "$test_bin/5.2Gene_coverage_of_all_samples.pl";
	my $Gene_cov_all_samples_file = $Gather_info_each_sample_dir.$pici."_All_sample_gene_cov_infor.txt";
	my $job1 = "Gene coverage of all samples";
	my $cmd1 = "perl $Gene_cov_all_samples $Pid $Gather_info_each_sample_file $Gene_cov_all_samples_file";
	run_cmd($cmd1,$job1);
}

#Sample information statistics
sub Sample_info_stat{
	my ($Samples) = shift @_;
	my @samples = split (/,/,$Samples);
	my $sample_num = $#samples+1;
	my $parallel = ($Samples_run_in_parallel >= $sample_num) ? $sample_num : $Samples_run_in_parallel;
	my $pm = Parallel::ForkManager->new($parallel);
	foreach my $sample(@samples){
		
		my $pid = $pm->start and next;
		#pos_depth_info
		my $sort_bam = "$output/Bowtie2/$sample/$sample.sort.bam";
		my $pos_depth_file = "$output/Count/$sample/$sample.depth.txt";
		my $pos_depth_info="$test_bin/4.1pos_depth_info.pl";
		my $Samtools_depth_tmp = $Samtools_depth;
		$Samtools_depth_tmp =~ s/\s/,/g;
		my $job = "Site depth information of $sample";
		my $cmd = "perl $pos_depth_info $Samtools_depth_tmp $D_Bed $sort_bam $pos_depth_file";
		run_cmd($cmd,$job,$pos_depth_file);
		
		#得到分段后是详细深度信息文件：CX2017090091-IEMF-Sciclone.low-coverage-block.txt
		
		my $covered_block="$test_bin/4.2low-covered-block_of_each-exon_of_panel-gene-2.pl";
		my $covered_block_dir = "$output/Low_depth";
		my $covered_block_file = "$output/Low_depth/$sample.low-coverage-block.txt"; 
		my $job1 = "Detailed depth information of $sample";
		my $cmd1 = "perl $covered_block -depth $pos_depth_file -cds $Mutalyzer -bed $D_Bed -len 2 -o $covered_block_dir";
		run_cmd($cmd1,$job1,$covered_block_file);
		
		##症状相关基因的覆盖度信息 20180507 增加
		my $symptom_covered_block = "$test_bin/4.2symptom_covered_block.pl";
		my $symptom_covered_block_file = "$output/Low_depth/$sample.low-coverage-block.symptom.txt"; 
		my $job1a = "Increase the coverage of clinical related genes of $sample";
		my $cmd1a = "perl $symptom_covered_block $Sampleinfo $IEM_gene_phe $covered_block_file $symptom_covered_block_file";
		run_cmd($cmd1a,$job1a);
		
	
		#得到比对q20,q30比对率文件2018030057-IEMS_map_state.txt统计文件 和gene_exon_cov.txt文件
		
		my $bt2_log = "$output/Bowtie2/$sample/$sample.bt2.log";
		my $state_trimomatic_infor="$test_bin/4.3state_trimomatic_infor-bed.pl";
		
		my $gene_exon_cov = "$output/Count/$sample/$sample"."_gene_exon_cov.txt";
		my $map_state = "$output/Count/$sample/$sample"."_map_state.txt";
		
		my $job2 = "Statistical data after mapped of $sample";
		my $cmd2 = "perl $state_trimomatic_infor $sample $bt2_log $pos_depth_file $D_Bed_format $GFF $gene_exon_cov $map_state";
		run_cmd($cmd2,$job2,$map_state);
		
		#性别判断
		my $Predict_sex="$test_bin/4.4Predict-sex.pl";
		my $bam_sorted_dup = "$output/Bowtie2/$sample/$sample.sort.dup.bam";
		my $Predict_sex_file = "$output/Count/$sample/$sample".".sex.txt";
		my $job3 = "Predicting the sex of $sample";
		my $cmd3 = "perl $Predict_sex $sample $bam_sorted_dup $Predict_sex_file";
		run_cmd($cmd3,$job3,$Predict_sex_file);
	
		$pm->finish;
	}
	$pm->wait_all_children;#等待所有子进程
	
	
}

#============================================================================
#dis_lib
sub Distribution_of_Library{
	my ($Pid) = shift @_;
	my $Distribution_of_Library="$test_bin/3.1count_library_distribution_R.pl";
	my $dis_lib_file = "$output/Count/${Pid}_All_sample_library_distribution.pdf";
	my $in_dir = "$output/Bowtie2";
	my $out_dir = "$output/Count";
	
	my $job = "Distribution of Library";
	my $cmd = "perl $Distribution_of_Library $in_dir $Pid $out_dir";
	run_cmd($cmd,$job,$dis_lib_file);
	`rm Rplots.pdf`;
}

#============================================================================
#mapping
sub Bowtie_2{
	my ($Samples,$Bowtie2Options,$Ref_gene) = @_;
	my @samples = split (/,/,$Samples);
	my $sample_num = $#samples+1;
	my $parallel = ($Samples_run_in_parallel >= $sample_num) ? $sample_num : $Samples_run_in_parallel;
	my $pm = Parallel::ForkManager->new($parallel);
	foreach my $sample(@samples){
		
		my $pid = $pm->start and next;
		
		my $raw_fastq1 = $KEYNAMES{$sample}->[11];
		my $raw_bianhao = $1 if ($raw_fastq1 =~ /.\/$sample(.+?)\.R1\.fastq/);
		my $Paired_fastq1 = "$Input/$sample"."$raw_bianhao".".R1.fastq";
		my $Paired_fastq2 = "$Input/$sample"."$raw_bianhao".".R2.fastq";
		
		# bowtie2_build
		my $bowtie_build_bin="$test_bin/2.1bowtie_build.pl";
		my $gonome_bt2_version=$1 if( $Ref_gene=~/.*\/(.+?)\.f(ast)?a$/i);
		my $gonome_bt2_file = "/work1/liulifeng/Genome/Mapped_index/bowtie2/$gonome_bt2_version";
		
		my $gonome_rev_bt2 = $gonome_bt2_file.".rev.2.bt2"; #最后生成的文件
		my $job=" bowtie2 Building genome index";
		my $cmd="perl $bowtie_build_bin $Ref_gene $gonome_bt2_file";
		run_cmd($cmd,$job,$gonome_rev_bt2);
		
		# bowtie2_mapping
		my $Bowtie_Mapping="$test_bin/2.2Bowtie_Mapping.pl";
		my $Bowtie2Options_tmp = $Bowtie2Options;
		$Bowtie2Options_tmp =~ s/\s/,/g;
		
		my $sample_bt2_file = "$output/Bowtie2/$sample/$sample.sam";#./Output/Bowtie2/2018030057-IEMS/2018030057-IEMS.sam
		my $bt2_log = "$output/Bowtie2/$sample/$sample.bt2.log";
		my $job2 = "Mapping sample to genome of $sample";
		my $cmd2 = "perl $Bowtie_Mapping $sample $gonome_bt2_file $Paired_fastq1 $Paired_fastq2 $Bowtie2Options_tmp $sample_bt2_file $bt2_log";
		run_cmd($cmd2,$job2,$bt2_log);
		
		# bowtie2_sam2bam
		my $bowtie2_sam2bam = "$test_bin/2.3bowtie2_sam2bam.pl";
		my $sam2bam_file = "$output/Bowtie2/$sample/$sample.bam";
		my $job3 = "Sam files are converted to BAM of $sample";
		my $cmd3 = "perl $bowtie2_sam2bam $sample_bt2_file $sam2bam_file";
		run_cmd($cmd3,$job3,$sam2bam_file);
		
		#sam2bam_sort
		my $sam2bam_sort = "$test_bin/2.4sam2bam_sort.pl";
		my $s2b_sort = "$output/Bowtie2/$sample/$sample.sort";
		my $s2b_s_bam = "$output/Bowtie2/$sample/$sample.sort.bam";
		my $job4 = "BAM files are sorted of $sample";
		my $cmd4 = "perl $sam2bam_sort $sam2bam_file $s2b_sort";
		run_cmd($cmd4,$job4,$s2b_s_bam);
		
		#sam2bam_sort_dup and index
		my $bam_sort_dup = "$test_bin/2.5bam_sort_dup.pl";
		my $s2b_sort_bam = $s2b_sort.".bam";
		my $bam_sorted_dup = "$output/Bowtie2/$sample/$sample.sort.dup.bam";
		my $job5 = "BAM files are repeated and indexed of $sample";
		my $cmd5 = "perl $bam_sort_dup $s2b_sort_bam $bam_sorted_dup";
		run_cmd($cmd5,$job5,$bam_sorted_dup);
		
		# Get the library size file
		my $get_lib_size = "$test_bin/2.6count_library_distribution_perl.pl";
		my $lib_size_file = "$output/Bowtie2/$sample.length.txt";
		my $job6 = "Get the library size file of $sample";
		my $cmd6 = "perl $get_lib_size $sample_bt2_file $lib_size_file";
		run_cmd($cmd6,$job6,$lib_size_file);
		
		#删除sam和bam文件
		#`rm -f $sample_bt2_file $sam2bam_file`;
		
		$pm->finish;
	}
	$pm->wait_all_children;#等待所有子进程
}

#============================================================================
#QC
sub Trimmomatic_QC{
	my $Samples = shift @_;
	my @samples = split (/,/,$Samples);
	my ($raw_fastq1,$raw_fastq2,$raw_bianhao) = ("","","");
	
	my $sample_num = $#samples+1;
	my $parallel = ($Samples_run_in_parallel >= $sample_num) ? $sample_num : $Samples_run_in_parallel;

	my $pm = Parallel::ForkManager->new($parallel);
	foreach my $sample(@samples){
		
		my $pid = $pm->start and next;
		
		$raw_fastq1 = $KEYNAMES{$sample}->[11];#q1=/work1/liulifeng/WSH_learn_20180408/IEMS/raw_read/2018030057-IEMS_S15.R1.fastq
		$raw_fastq2 = $KEYNAMES{$sample}->[22];
		$raw_bianhao = $1 if ($raw_fastq1 =~ /.\/$sample(.+?)\.R1\.fastq/);
		
		my $Paired_fastq1 = "$Input/$sample"."$raw_bianhao".".R1.fastq";
		my $Paired_fastq2 = "$Input/$sample"."$raw_bianhao".".R2.fastq";
		my $Unpaired_fastq1 = "$Input/Unpaired/$sample"."$raw_bianhao".".R1.fastq";
		my $Unpaired_fastq2 = "$Input/Unpaired/$sample"."$raw_bianhao".".R2.fastq";
		
		##运行Trimmomatic软件
		my $Trimmomatic_QC="$test_bin/./1.1Trimmomatic_QC.pl";
		my $job="filter flow quality reads of $sample";
		
		my $Trim_para_1_tmp = $Trim_para_1;
		$Trim_para_1_tmp =~ s/\s/,/g;
		my $Trim_para_2_tmp = $Trim_para_2;
		$Trim_para_2_tmp =~ s/\s/,/g;
		
		my $cmd="perl $Trimmomatic_QC $ReadtypeOptions $Trim_para_1_tmp $raw_fastq1 $raw_fastq2 $Paired_fastq1 $Unpaired_fastq1 $Paired_fastq2 $Unpaired_fastq2 $Trim_para_2_tmp";
		run_cmd($cmd,$job,$Unpaired_fastq2);
		
		#计算reads质量情况
		my $QC_count_trimmomatic = "$test_bin/./1.2QC_count_trimmomatic.pl";
		my $job2="Statistical reads quality of $sample";
		my $QC_count_outfile = "$output/Count/$sample/$sample"."_QC_count.txt";
		
		my $phred33 = $1 if ($Trim_para_1_tmp =~ /phred(.+?)$/);
		my $cmd2="perl $QC_count_trimmomatic $sample $raw_fastq1 $raw_fastq2 $Paired_fastq1 $Paired_fastq2 $phred33 $QC_count_outfile";
		run_cmd($cmd2,$job2,$QC_count_outfile);
		
		$pm->finish;
	}
	$pm->wait_all_children;#等待所有子进程结束
}
#============================================================================
#创建文件夹
sub Create_DIR{
	my $Samples = shift @_;
	mkdir("$Input") unless (-e "$Input");
	mkdir("$Input/Unpaired") unless (-e "$Input/Unpaired");
	mkdir("$output") unless (-e "$output");
	mkdir("$output/Bowtie2") unless (-e "$output/Bowtie2");
	mkdir("$output/Count") unless (-e "$output/Count");
	mkdir("$output/Low_depth") unless (-e "$output/Low_depth");
	mkdir("$output/SNP") unless (-e "$output/SNP");
	mkdir("$output/Combined_result") unless (-e "$output/Combined_result");

	
	my @samples = split(/,/,$Samples);
	
	foreach my $sample(@samples){
		mkdir("$output/Bowtie2/$sample") unless (-e "$output/Bowtie2/$sample");
		mkdir("$output/Count/$sample") unless (-e "$output/Count/$sample");
	}
}

#===================================run_cmd运行命令==========================
sub run_cmd {
	my ($cmd, $job, $out_result) = @_;
	
	my ($start,$start_h) = get_start_time();
	print STDERR "Start $job at $start_h ...\n";
	print LOG "#$job:\n";
	print LOG "$cmd\n";
	
	my $ret = 0;
	if ($out_result && -e $out_result) {
		print STDERR "$out_result exists, $job skipped!\n";
	}
	else {
		$ret = system($cmd);
	}
	
	my ($end,$end_h) = get_end_time();
	my $dur = get_duration($start,$end);
	
	my $status = "Done";
	if ($ret != 0) {
		$status = "Stop";
	}
	
	print STDERR "$status at $end_h.\n";
	print STDERR "Duration: $dur\n\n";
	print LOG "#Run status: $status\n\n";

	if ($ret != 0) {
		exit(-1);
	}
}

# 该步骤持续时间
sub get_duration {
	my ($start,$end)=@_;
	my $sec = $end - $start;
	my $days = int($sec/(24*60*60));
	my $hours = ($sec/(60*60))%24;
	my $mins = ($sec/60)%60;
	my $secs = $sec%60;
	
	return "$days days $hours hours $mins minutes $secs seconds";
}

sub get_start_time {
	my $start = time();
	my $start_h = localtime();
	return ($start,$start_h);
}

sub get_end_time {
	my $end = time();
	my $end_h = localtime();
	return ($end,$end_h);
}


#============================================================================
#检查参数
sub CHECKING_PARA_SETS{
	
	my $CHECKING_PARA_SETS="$test_bin/CHECKING_PARA_SETS.pl";
	my $job="CHECKING_PARA_SETS";
	my $cmd="perl $CHECKING_PARA_SETS $parameter";
	run_cmd($cmd,$job);
	
	my $check_para_stat_log="./check_para_sets.log";
	
	my $pass_flag=0;
	open(IN,$check_para_stat_log) or die "No file of check_para_stat_log:$check_para_stat_log be found! $!(ACGT101_lncRNA)";
	while(<IN>){
		chomp;
		if(/unpassed/){
			print STDERR "ERRORS detected, please checking the log file of ./check_para_sets.log!\n";		
			print STDERR "Quiting the Pipeline!\n";
			
			$pass_flag=1;
			exit 1;
		}
	}
	close(IN);
	
	if($pass_flag==0){
		`rm $check_para_stat_log`;  # 参数正常，删除参数
		print STDERR "Starting pipeline!\n";
		print STDERR "#====================================\n";
	}
}


#============================================================================
sub TIME{
	my ($min,$hour,$day,$month,$year)=(localtime)[1,2,3,4,5];
	($min,$hour,$day,$month,$year) = (sprintf("%02d", $min),sprintf("%02d", $hour),sprintf("%02d", $day),sprintf("%02d", $month + 1),$year + 1900);
	my $run_start_time="$year$month$day$hour$min";
	return $run_start_time;
}

#============================================================================
sub get_params{
	my $parameter=shift;
	my %paramfile_hash=();
	my @sample= ();
	my $IN=open_IN($parameter);
	while(<$IN>){
		next if ($_ =~ /^\#/);
		next unless ($_ =~ /=/);
		chomp $_;
		my($key, $value) = split('=', $_,2);
		$key=~s/\s//g;
		$value=~s/\s+$//;
		$value=~s/^\s+//;
		$paramfile_hash{$key} = $value;
		
		#获得样本信息
		if (/^KeyName\s*=\s*(\S+)/i) { #\S匹配非空白字符
			$KeyName= $1;
			push @sample, $KeyName;
			print STDERR "Extracting $KeyName and fastq file information\n";
		}
		if (/^q1\s*=\s*(.*\/(.+))$/xms) {
			$Fastq1_file = $1; #/work1/liulifeng/WSH_learn_20180408/IEMS/raw_read/2018030057-IEMS_S15.R2.fastq
			$Fastq1 = $2;      #2018030057-IEMS_S15.R1.fastq
			die "can't find file: $Fastq1_file" if(!-e $Fastq1_file);
			$KEYNAMES{$KeyName}->[1] = $Fastq1;				#没啥用
			$KEYNAMES{$KeyName}->[11] = $Fastq1_file;
		}elsif (/^q2\s*=\s*(.*\/(.+))$/xms) {
			$Fastq2_file = $1;
			$Fastq2 = $2;
			die "can't find file: $Fastq2_file" if(!-e $Fastq2_file);
			$KEYNAMES{$KeyName}->[2] = $Fastq2;
			$KEYNAMES{$KeyName}->[22] = $Fastq2_file;
		}elsif (/^q\s*=\s*(.*\/(.+))$/xms) {
			$Fastq_file = $1;
			$Fastq = $2;
			die "can't find file: $Fastq_file" if(!-e $Fastq_file);
			$KEYNAMES{$KeyName}->[3] = $Fastq;
			$KEYNAMES{$KeyName}->[33] = $Fastq_file;
		}
	}
	close $IN;
	my $Samples;
	$Samples=join",",@sample;
	print STDERR "done!\n";
# panel
	$Pid = $paramfile_hash{Pid};
# Samples_run_in_parallel
	$Samples_run_in_parallel = $paramfile_hash{Samples_run_in_parallel};
# Reads type
	$ReadtypeOptions = $paramfile_hash{ReadtypeOptions};
# QC
	$Adapter_3 = $paramfile_hash{Adapter_3};
	$Adapter_5 = $paramfile_hash{Adapter_5};
	$FilterType = $paramfile_hash{FilterType};
	$Trim_para_1 = $paramfile_hash{Trim_para_1};
	$Trim_para_2 = $paramfile_hash{Trim_para_2};
# mapping
	$Bowtie2Options = $paramfile_hash{Bowtie2Options};
	$Ref_gene = $paramfile_hash{Ref_gene};
# depth	
	$Samtools_depth = $paramfile_hash{Samtools_depth};
# depth data
	$Mutalyzer = $paramfile_hash{Mutalyzer};
	$D_Bed = $paramfile_hash{D_Bed};
	$D_Bed_format = $paramfile_hash{D_Bed_format};
	$GFF = $paramfile_hash{GFF};
	$IEM_gene_phe = $paramfile_hash{IEM_gene_phe};
# 7.Gather_info_each_sample
	$Sampleinfo = $paramfile_hash{Sampleinfo};
# 8.samtools mpileupMpileup
	$Mpileup = $paramfile_hash{Mpileup};
	$Genome = $paramfile_hash{Genome};
# 9.annovar
	$S_annovar = $paramfile_hash{S_annovar};
	$S_humandbdir = $paramfile_hash{S_humandbdir};
	$Xref = $paramfile_hash{Xref};
# 10.Local_freq
	$Local_freq = $paramfile_hash{Local_freq};
# 11.format1
	$R_HGMD = $paramfile_hash{R_HGMD};
	$Mutalyzer_All = $paramfile_hash{Mutalyzer_All};
	$EXAC = $paramfile_hash{EXAC};
	$gnomad = $paramfile_hash{gnomad};
	$PP5 = $paramfile_hash{PP5};
	$BP6 = $paramfile_hash{BP6};
	$Model = $paramfile_hash{Model};
	$GeneDisease = $paramfile_hash{GeneDisease};
	$Phe = $paramfile_hash{Phe};
	$OMIM = $paramfile_hash{OMIM};
	$Clinvar = $paramfile_hash{Clinvar};
	$GHR = $paramfile_hash{GHR};
	$Hgmd = $paramfile_hash{Hgmd};
# 12.format1
	$G_Trans_gb = $paramfile_hash{G_Trans_gb};
	$Gene_info_cds = $paramfile_hash{Gene_info_cds};
	$Verbose = $paramfile_hash{Verbose};
	$local_pd = $paramfile_hash{local_pd};
	$GeneList = $paramfile_hash{GeneList};

####
}

sub open_IN {
	my ($file) = @_;
	
	my $IN = ();
	if ($file =~ /.gz$/i) {
		open($IN, "gzip -cd $file |") or die "could not open $file $!\n";		
	}else {
		open($IN, "$file") or die "could not open $file $!\n";
	}
	return $IN;
}