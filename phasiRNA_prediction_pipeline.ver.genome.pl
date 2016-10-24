#!/usr/local/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use Data::Dumper;
use File::Basename;

use lib "/data2/homes/pingchuan/lib";
use Parallel::ForkManager;

######################################################################
#  PhasiRNA Prediction Pipeline (PPP)
#
#  Copyright 2013 University of Delaware
#
#  Pingchuan Li: li@dbi.udel.edu
#  @ Blake Meyers lab
#
#  Department of Plant Soil Science
#  Delaware Biotechnology Insititute
#  University of Delaware
#  Newark DE, 19711
#
#  Maintenanced by Pingchuan Li 
#  
#  This program is a FREE software; you can redistribute it and/or
#  modify as your wishes
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
######################################################################

my $version  = 0.296;
Ptime("games start...");

# make sure Perl is over 5.12
no 5.12;


# NOTE:
# In order to do the analysis, the format of samall RNA in fasta should be looks like >seq_xxxx|abun_3
# now the alignment engine only support the bowtie, because the output of bowtie is very common for further analysis
# The script is based under Ho-Ming Chen's original paper, PNAS, 2007, V104, No.9 3318-3323
# you'll be appreciated to report any of issue to Pingchuan Li at: lipingchuan@gmail.com 

# This perl script is developed to identify small RNAs clustered in 21-nt ~ 24-nt increments or any other size
# As for the original method for scoring the phasiRNA, please refer to Ho-Ming Chen et al., PNAS, 2007, V104, No.9 3318-3323

# 2013-05-30 normalized the phased siRNA abundance
# 2013-05-31 adding the 231bp sliding window to calc the p value, which give all the phased siRNA located in the same cluster pretty similar score, which benefit the boundary finding of cluster
# 2013-05-31 adding the hits number for each of reads in the output of cluster
# 2013-06-01 v.22 adding the sliding windows calc support for scoring a certain sRNA
# 2013-06-01 v.23 adding the filter for the interested phased siRNA abundance -n (noise)
# 2013-06-02 v.24 adding more levels for the cutoff
# 2013-06-05 v.25 adding the filter for the average of hits for phased siRNA
# 2013-06-10 v.26 adding 2 filter in the command line, t and ht
# 2013-06-14 v.27 enforce the k value to 21 when it's bigger than 21 otherwise the scritp will get error infor
# 2013-06-24 v.28 add a qualifiation filter for the cluster in terms of the sRNA amount, >=3 by default
# 2013-06-28 v.29 finish the whole thing
# 2013-06-29 v.291 bug fixes.
# 2013-07-01 v.294 change the data structure for the scoring which dramatically increase the scoring speed .
# 2013-07-11 v.295 bug fix, the default of value for n and t should be added the 'defined' 
# 2013-07-12 v.296 feature added, add bowtie alignment reports limitation by -k options

my %sum = ();

GetOptions(my $options = {},
              "-i=s","-f=s","-d=s","-p=s","-cpu=i","-rl=s","-px=s","-n=i","-g=i","-t=i","-ht=i","-q=s","-k=i"#"-o=s","-s=s","-a=s",
);
########### USAGE ##############
my $USAGE = <<USAGE;
Scripts: PhasiRNA Prediction Pipeline (PPP). The sRNA readswere mapped
by bowtie by default.

Version: $version, Written by Pingchuan Li \@ Blake Meyers Lab 04/2013
Arguments:

    Mandatory:
        -i        sRNA input file name
        -f        two files formats are supported,[t] for tagcount, [f] for fasta. 
        -d        indexed genome by bowtie, indicating the prefix of the genome
        -rl       register len(rl), such as 21 or 24, separated by comma like: 21,24 or 21,22,23,24 etc
        -px       prefix of all the output file, in order to distinguish for each other

    optional:
        -k        bowtie alignment hits report, default = all
        -q        degradome or PARE data in ONLY tagcount format
        -m        mismatches for both the sRNA and PARE/degradome alignment by bowtie, default = 0
        -p        p-value in decimal number, defult = 0.005;
        -cpu      cpu numbers for the bowtie alignment
        -n        noise, default = 1, for those which have abundance less equal than 1, properly increase
                  noise value for union dataset	
        -g        gap between two separate cluster, 300bp by default
        -t        minimal proportation of the interested register small RNA abundance, default = 85
        -ht       the maximal average of hits for the small RNA of a certain cluster, defaut = 10

USAGE
die $USAGE unless defined ($options->{rl});#
die $USAGE unless defined ($options->{px});

my $prefix = ($options->{px})?$options->{px}:'Results';
my @tchr = ();
my @temp_files = ();
my $bowtie_output = generate_rand_filename(11);
my $total_r_abundance = 0;  #only hit small RNA will be counted.
my @rl = split/\,/,$options->{rl};
my $p_cutoff = (defined $options->{p})?($options->{p} + 0):0.005;
my $noise = (defined $options->{n})?$options->{n}:1;
my $ta_cf =(defined $options->{t})?$options->{t}:85; #in a certain cluster, the minimal abundance of interested small RNA proportation.
my $avg_hit = (defined $options->{ht})?$options->{ht}:10; # the average hits of a phased siRNA within a cluster , in order to filter the TE cluster sites if they are not real
my $gap_len = (defined $options->{g})?$options->{g}:300;
$options->{a} = 'm';
my $pare_data = ($options->{q})?$options->{q}:0;
my $cpu = ($options->{cpu})?$options->{cpu}:2;
my $mm = (defined $options->{m})?$options->{m}:0;
my $input = $options->{i};
my $k = ($options->{k})?"-k $options->{k}" :"-a";  #bowtie options for hits reports limitation.


$options->{o} = ($options->{o})?$options->{o}:"output";
$options->{s} = ($options->{s})?$options->{s}:'score';

my %sh = ();#define the Small RNA Hash data. use the unique combination of information as the key.
my %hits = (); # this will record the hits number for each of the small RNAs
my %qualified_cluster = ();
my %boundary = (); # this will record the boundary of each predicted phasiRNA cluster
my %pre_scored_sRNA = (); #this will pre-saving the score of each siRNA in a sliding window


# the below two hash will record validate cluster by pare data and unvalided cluster if possible;
my %candidate = ();
my %miscandidate = ();

#when the motif is very close to both cluster,it's hard to judge which cluster has this cleavage site.
my %candidate_vague = ();

push(@temp_files,$prefix.".".$bowtie_output);
# -----------tag count processing if the input is tagcount---------------
my $t_i = 0;
if ($options->{f} =~ /t/i) {
	my $fasta_file = generate_rand_filename(11);
	push(@temp_files,$fasta_file);
	
	open(IN, $options->{i});
	open(OUT,">$fasta_file");
	while (<IN>){
		chomp;
		$t_i++;
		my($seq,$abun) = split/\s+/,$_;
		print OUT ">seq_$t_i|$abun\n$seq\n";
	}
	close IN;
	close OUT;
	push(@temp_files, $fasta_file);
	$input = $fasta_file;
	Ptime("seq conversion is done...");
}

Ptime("aligning...");
Ptime("script version is $version..");

my @parameters = join(" ", "-f",  "$k -n $mm", "-p $cpu", "$options->{d}", "$input", " >$prefix.$bowtie_output");
my @bowtie_out = system("bowtie @parameters");
Ptime("aligning is done");

#put all the signature and their information into the HASH %sh
#--------bowtie out put pre processing -------------------
open(DATA, "$prefix.$bowtie_output") || die "Cannot open the bowtie outputed DATA";
while(my $line=<DATA>) {# input small RNA mapping data
	chomp $line;
	my ($seqid, $strand, $target, $pos,  $seq, @other)=split(/\t/, $line);
	my $abun = 0;
	

	# calc the individual abundance
	if ($seqid =~ /\|(\d+)/) {
		$abun = $1;
	}
	elsif($seqid =~ /abun\_(\d+)/) {
		$abun  = $1;
	}
	else {
		die "we can't locate the abundance from the seq_id";
	}
	
	# calc the total mapped reads abundance, the noise will also be summed regardless of the sRNA abundance.
	$total_r_abundance += $abun unless ($hits{$seqid});
	$hits{$seqid}++;

	# in order to improve the sensitivity, run the below filter.
	next if ($abun <= $noise);
	
	($pos,$seq) = bowtie_process($strand,$pos,$seq);
		
	my $chr;
	if ($target =~ /(\d+)/) {
		$chr = $1; # incase $chr equal to '01' or '02', what we need is '1','2' etc
		$chr += 0;
	}else {
		print "make sure all the chr including the chloroplast and mitchondria are also initiated by digital number\n";
		die "the chromosome is somoehow not named as number, check it up\n";
	}
	#if the same position of the same chromosome contains multiple small RNA in different length or strand, they will be put together
	push(@{$sh{$chr}->{$pos}->{$strand}},[$seqid,$seq,$abun]);
	
}
Ptime("finished the loading of sRNA alignment data.");
# ------------end of preprocessing ---------------




my @pre_splitted_cluster = ();   #@pre_splitted_cluster = ( [[p1],[p2.1],[p2.2],[p3]], [[p1]],  [ [p1],[p2] ] );
#-----------------------sorting, divided  the cluster prior to the calculation ------------------------------
foreach my $chr (sort {$a <=> $b} keys %sh) {
	my $i = 0;
	my $temp_pos = 0;
	foreach my $pos (sort {$a <=> $b} keys %{$sh{$chr}}) {
		foreach my $strand ( keys %{$sh{$chr}->{$pos}}) {
			my @sRNA_ref = @{$sh{$chr}->{$pos}->{$strand}};
			foreach my $sRNA_ref (@sRNA_ref) {
				my ($seqid,$seq,$abun) = @{$sRNA_ref};
				$i++;
				
				#initialization
				if ($i == 1) { #the first sRNA in a chr
				    $temp_pos = $pos;
				    push(@pre_splitted_cluster, [[$chr,$pos,$strand,$seqid,$seq,$abun]]);
				}
				
				else {
					# judge the distance for two sRNA
					# if shorter than the pre-defined distance
					if ($pos - 400 <= $temp_pos) {
						
						push(@{$pre_splitted_cluster[-1]},[$chr,$pos,$strand,$seqid,$seq,$abun]);					
						$temp_pos = $pos;
						
					} else {
						push(@pre_splitted_cluster, [[$chr,$pos,$strand,$seqid,$seq,$abun]]);
						$temp_pos = $pos;
						
					}										
				}
			}
		}
	}
}

Ptime("finished the sorting of presplitted sRNA clusters data");

my $total_cluster = 0;
my @qualified_cluster;
my %sh2; # this is a selected sh from sh2, structure: push(@{$sh{$chr}->{$pos}->{$strand}},[$seqid,$seq,$abun]);

foreach my $cluster_ref (@pre_splitted_cluster) {
	my @record_ref = @{$cluster_ref};
	if ($#record_ref<= 7) {
		# this cluster is not qualified
	}else {
		$total_cluster++;
		#push(@qualified_cluster, $cluster_ref);
		foreach my $sRNA_ref (@{$cluster_ref}) {
				my ($chr,$pos,$strand,$seqid,$seq,$abun) = @{$sRNA_ref};
				push(@{$sh2{$chr}->{$pos}->{$strand}},[$seqid,$seq,$abun]);
		}
	}
}

Ptime("totoal $total_cluster presplitted sRNA clusters");

# calc the total chromosome number --------------
foreach my $chr (sort {$a <=> $b} keys %sh2) {
	push(@tchr,$chr);
}

# _____   _____   _____   _____    _   __   _   _____  
#/  ___/ /  ___| /  _  \ |  _  \  | | |  \ | | /  ___| 
#| |___  | |     | | | | | |_| |  | | |   \| | | |     
#\___  \ | |     | | | | |  _  /  | | | |\   | | |  _  
# ___| | | |___  | |_| | | | \ \  | | | | \  | | |_| | 
#/_____/ \_____| \_____/ |_|  \_\ |_| |_|  \_| \_____/ 
# -------start to score each of reads in different size---
foreach my $rl (sort {$a <=> $b} @rl) {#rl, register length
	
	
	
	# fork the processing, the nubmer is eqaul to the chromosomes
	my $pm = new Parallel::ForkManager($#tchr + 1); 
	foreach my $chr(sort {$a<=>$b} @tchr) {
		my $pid = $pm->start and next;
		
		phasing_analysis_mul($chr, $rl);
		
		$pm->finish;
	}
	$pm->wait_all_children;
	
	my @all_out_chr  = glob"$prefix.$options->{o}\_all_sRNA_$rl\_out.chr*.txt";
	my @qualified_out = glob"$prefix.$options->{s}\_p$p_cutoff\_sRNA_$rl\_out.chr*.txt";
	
	merge_files(\@all_out_chr,"$prefix.$options->{o}\_all_sRNA_$rl\_out.txt");
	merge_files(\@qualified_out,"$prefix.$options->{s}\_p$p_cutoff\_sRNA_$rl\_out.txt");

	#----------------cluster the outputed data-----------------------
	my @cutArrange = (5e-3,1e-3,5e-4,1e-4,5e-5,1e-5,5e-6,1e-6,5e-7,1e-7);
	foreach my $cut (sort {$b <=> $a} @cutArrange) {
		
		my %cluster = ();#record all the qualified small RNA under the p-value cutoff
		my $cluster = 0;
		my $phaseNumb = 0;
		my $i = 0;
		my $temp_chr = 0;
		my $temp_pos = 0;
		my $temp_row = '';
		
		# open the qualified scored output
		open(IN,"$prefix.$options->{s}\_p$p_cutoff\_sRNA_$rl\_out.txt");
		open(TEMP1,     ">$prefix.temp.sorted.$cut.cluster.minimal.$rl.requirement.txt");
		push(@temp_files,"$prefix.temp.sorted.$cut.cluster.minimal.$rl.requirement.txt");
		while (<IN>) {
			chomp;
			my @array = split/\t/,$_;
			if ($array[-1] <= $cut) {
				print TEMP1 "$_\n";
			}
		}
		close TEMP1;
		close IN;

		open(TEMP2,"$prefix.temp.sorted.$cut.cluster.minimal.$rl.requirement.txt");
		while (<TEMP2>) {
		
				chomp;$i++;
				my (@array) = split/\t/,$_;
				if ($i == 1) {#this sub will only execute once
						$cluster++;
						
						push(@{$cluster{$cluster}},[$array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$array[6],"n=$array[7]","k=$array[8]","hts=$hits{$array[3]}",$array[9],$array[10]]);
						
						$temp_chr = $array[0];
						$temp_pos = $array[2];
						next;
				}
				if ($temp_chr == $array[0]) {
						#judge the pos, 2 possibilities
						if ($array[2] - $gap_len <= $temp_pos) {
								
								push(@{$cluster{$cluster}},[$array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$array[6],"n=$array[7]","k=$array[8]","hts=$hits{$array[3]}",$array[9],$array[10]]);
								
								$temp_chr = $array[0];
								$temp_pos = $array[2];
						}
						else {#gap is bigger than predified gap length, such as 300bp
								
								$cluster++;

								push(@{$cluster{$cluster}},[$array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$array[6],"n=$array[7]","k=$array[8]","hts=$hits{$array[3]}",$array[9],$array[10]]);
								
								$temp_chr = $array[0];
								$temp_pos = $array[2];
						}
				}
				else {
						#output the data to the last, update the cluste #
						$cluster++;

						push(@{$cluster{$cluster}},[$array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$array[6],"n=$array[7]","k=$array[8]","hts=$hits{$array[3]}",$array[9],$array[10]]);

						$temp_chr = $array[0];
						$temp_pos = $array[2];
				}
		}
		close TEMP2;
		
		# the check will save those cluster which  mainly have the interested small RNA
		cluster_abundance_check(\%cluster,$rl,$ta_cf,$cut);
	}
}

Ptime( "total sRNA abundance = $total_r_abundance");
unlink "$prefix.temp.sorted.cluster.minimal.requirement.txt";

# -----------output the summary file ----------------
open(SUM,">$prefix.summary.sRNA.v$version.txt");
print SUM "Note:the algorithm is baed on the -$options->{a}\n";
foreach my $cutoff (sort {$b <=> $a}keys %sum) {
	print SUM "$cutoff\t";
	foreach my $size (sort {$a <=> $b} keys %{$sum{$cutoff}}) {
		my $C_num = ($sum{$cutoff}->{$size}->{'phaseClus'})?$sum{$cutoff}->{$size}->{'phaseClus'}:0;
		my $P_num = ($sum{$cutoff}->{$size}->{'phaseNumb'})?$sum{$cutoff}->{$size}->{'phaseNumb'}:0;
		print SUM "$size\t$C_num\t$P_num\t\t";
	}
	print SUM "\n";
}
close SUM;
close DATA;


open(CLUSTER_LITE,">$prefix.cluster.boundary.without.PARE.validation.list") or die "cant create the lite cluster file\n"; #$prefix.valid.cluster.by.PARE.txt
# -------------------------------- output the cluster-------------------------------------------
foreach my $cut (sort {$a <=> $b} keys %qualified_cluster) {
	foreach my $rl (sort {$a <=> $b} keys %{$qualified_cluster{$cut}}) {
		open(OUT_CLU,">$prefix.$options->{s}\_p$cut\_sRNA_$rl\_out.cluster") or die "We cant create the file\n";
		foreach my $i  (sort {$a <=> $b} keys %{$qualified_cluster{$cut}->{$rl}}) {

			my $cluster_ref  = $qualified_cluster{$cut}->{$rl}->{$i};

			# the below are only for v3.2
			my ($start,$end,$chr) = boundary($cluster_ref,$rl);
			push(@{$boundary{$rl}->{$cut}->{$chr}}, [$start,$end]);
			
			# ----------------------------
			print OUT_CLU ">cluster = $i | chr = $chr and pos between $start and $end\n";
			print CLUSTER_LITE "$cut|$rl|chr = $chr:$start..$end\n";

				foreach my $ref (@{$cluster_ref}) {
					my ($chr,$str,$pos,$seq_id,$seq,$len,$abun,$n,$k,$hts,$best_pos,$pvalue) = @{$ref};
					print  OUT_CLU join("\t","o",$chr,$str,$pos,$seq_id,$seq,$len,$abun,$n,$k,$hts,$best_pos,$pvalue);
					print  OUT_CLU "\n";
				}
		}print CLUSTER_LITE "\n";
		close OUT_CLU;
	}print CLUSTER_LITE "\n";
}
close CLUSTER_LITE;


# _____       ___   _____    _____        _____   _____    _____   _____   _____   _____   _____   _   __   _   _____  
#|  _  \     /   | |  _  \  | ____|      |  _  \ |  _  \  /  _  \ /  ___| | ____| /  ___/ /  ___/ | | |  \ | | /  ___| 
#| |_| |    / /| | | |_| |  | |__        | |_| | | |_| |  | | | | | |     | |__   | |___  | |___  | | |   \| | | |     
#|  ___/   / / | | |  _  /  |  __|       |  ___/ |  _  /  | | | | | |     |  __|  \___  \ \___  \ | | | |\   | | |  _  
#| |      / /  | | | | \ \  | |___       | |     | | \ \  | |_| | | |___  | |___   ___| |  ___| | | | | | \  | | |_| | 
#|_|     /_/   |_| |_|  \_\ |_____|      |_|     |_|  \_\ \_____/ \_____| |_____| /_____/ /_____/ |_| |_|  \_| \_____/ 

my %PARE_alignment = (); # this hash will record abundance of each position
if ($pare_data) {# the pare format is tag only
	
	Ptime("starting to processing PARE data...");
	# --------------------------------PARE bowtie alignmen -------right now only support the 21nt register----- 
	my $pare_fas = generate_rand_filename(11);
	my %pare_summary = ();# structure $pare_summary{$cutoff}->{$siz} = $amount
	my ($pare_name,$path,$suffix) = fileparse($pare_fas,qr/\.[^.]*/);
	
	push(@temp_files,$pare_fas);# for the convenience of removing later
	
	open(IN, $pare_data);
	open(OUT,">$pare_fas");
	while (<IN>){
		chomp;
		$t_i++;
		my($seq,$abun) = split/\s+/,$_;
		print OUT ">seq_$t_i|$abun\n$seq\n";
	}
	close IN;
	close OUT;
	
	my @parameters = join(" ", "-f",  "$k -n $mm", "-p $cpu", "$options->{d}", "$pare_fas", " >$prefix.pare.output");
	my @bowtie_out = system("bowtie @parameters");
		
	# --------------------- parsing the PARE alignment, generate the data structure---------------------------
	# ----------------all the pare information will be merged and saved into the %PARE_alignment--------------
	my $total_pare_abun;
	open(PAREOUT,"$prefix.pare.output") or die "cant open the PARE out";
	while (<PAREOUT>) {
		#seq_50051|9     +       chr2    33727470        AAAAAAAAAAATCCCCCCAA    IIIIIIIIIIIIIIIIIIII    0
		chomp;
		
		my ($seqid, $strand, $target, $pos,  $seq, @other)=split/\t/, $_;
		my ($abun) = $seqid =~ /\|(\d+)/;
		die "not a number, abundance 4313\n" unless (looks_like_number($abun));
		my $total_pare_abun += $abun;
		my ($chr) = $target =~ /(\d+)/;
		die "not a number,chromosome, 14343\n"  unless (looks_like_number($chr)); $chr+= 0;
		my ($pos_new,$seq_new) = bowtie_process($strand,$pos,$seq);
		
		$PARE_alignment{$chr}->{$pos_new}->{$strand}->{'abun'} += $abun;
		$PARE_alignment{$chr}->{$pos_new}->{$strand}->{'seq'}   = pare_seq_combine($PARE_alignment{$chr}->{$pos_new}->{$strand}->{'seq'}, $seq_new);
	}
	close PAREOUT;
	
	# output the parsed PARE data with location, abundance , strand etc.
	Ptime("starting to output parsed result...");
	
	open(PARE_PARSE,">$pare_name.PARSED.txt");
	
	foreach my $chr (sort {$a <=> $b} keys %PARE_alignment) {
		foreach my $pos (sort {$a <=> $b} keys %{$PARE_alignment{$chr}}) {
			foreach my $strand (keys %{$PARE_alignment{$chr}->{$pos}}) {
				my $abun = ($PARE_alignment{$chr}->{$pos}->{$strand}->{'abun'})?($PARE_alignment{$chr}->{$pos}->{$strand}->{'abun'}):"NA";
				my $seq  = ($PARE_alignment{$chr}->{$pos}->{$strand}->{'seq'} )? $PARE_alignment{$chr}->{$pos}->{$strand}->{'seq'}  :"NA";
				print PARE_PARSE "$chr\t$pos\t$strand\t$abun\t$seq\n";
			}
		}
	}
	close PARE_PARSE;

	# ------------------------judge whether there is correlationship between the pare processing and phasiRNA prediciton, save the positive
	# because almost all of the PARE were loated outside the boundary, thus we only consider a relative range between 0 and 168bp to detect the most abundant PARE signatures
		
	# data structure: push(@{$boundary{$rl}->{$cut}->{$chr}}, [$start,$end]);
	foreach my $rl (sort {$a <=> $b} keys %boundary) {
			foreach my $cut (sort {$a <=> $b} keys %{$boundary{$rl}}) {
				foreach my $chr (sort {$a <=> $b} keys %{$boundary{$rl}->{$cut}}) {
					my @cluster_ref =@{$boundary{$rl}->{$cut}->{$chr}};
					while (my ($idx,$region_ref) = each @cluster_ref) {
						# retrieve the original boundary information without the assistant of PARE/degradome data.
						my ($cluster_start, $cluster_end) = @{$region_ref};
							
						my @top_three_pare_ref = &top_pare($chr,$cluster_start - $rl*4 - 2, $cluster_end + $rl*4 + 2, 3);
							
						# the stastics the 96.9% pare located within 3 ~ 4 sRNA size
						my $total_pare_abun_in_an_wider_range = pare_abun_total($chr, $cluster_start - $rl*4 - 2, $cluster_end + $rl*4 + 2); # incase of 2 nt overhang, 4 register length plus the 2nt overhang will be the maximal range to detect the total amount of the PARE signatures
						my $prime_5_outer_abun                = pare_abun_total($chr, $cluster_start - $rl*4 - 2, $cluster_start + $rl + 2);
						my $prime_3_outer_abun                = pare_abun_total($chr, $cluster_end - $rl -2, $cluster_end + $rl*4 + 2);
							
						#print "$chr,$cut,$cluster_start,$cluster_end, $total_pare_abun_in_an_wider_range, $prime_5_outer_abun, $prime_3_outer_abun\n";						
						if ($total_pare_abun_in_an_wider_range > 0) {
							my $motif_prime;
							if ($prime_5_outer_abun >= 6 or $prime_3_outer_abun>= 6) {
								
									if ($prime_5_outer_abun == $prime_3_outer_abun) {
										# unable to confirm which strand is the transcriptional strand and where the motif located
										push(@{$miscandidate{$rl}->{$cut}->{$chr}},[$cluster_start, $cluster_end,'UC','UC']);
									}
									else {
										
										my @a1 = top_pare($chr,  $cluster_start - $rl*4 - 2, $cluster_start + $rl   + 2, 1);
										my @a2 = top_pare($chr,  $cluster_end   - $rl   - 2, $cluster_end   + $rl*4 + 2, 1);
										
										my $top_1_pr5 = $a1[0];
										my $top_1_pr3 = $a2[0];
										
										#$j will indicate whether it's top3, $flag will indicate whether it's registered to the specified ends
										#initializing the $j5 $j3, $flag1 $flag2
										my $j5 = 0; my $j3 = 0; my $flag1 = 0; my $flag2 = 0;
										
										if ($top_1_pr5) {
											($j5, $flag1) = &top3_phase_judge(\@top_three_pare_ref, $top_1_pr5, 'ups', $cluster_start, $cluster_end, $rl); 
										}
										
										if ($top_1_pr3) {
											($j3, $flag2) = &top3_phase_judge(\@top_three_pare_ref, $top_1_pr3, 'dws', $cluster_start, $cluster_end, $rl);
										}
										
										if ($j5 + $j3 + $flag1 + $flag2 == 4) {
											#two?
											push(@{$miscandidate{$rl}->{$cut}->{$chr}},[$cluster_start, $cluster_end,'h2','h2']); # either top3 or phased or not
										}
										elsif ($j5 + $flag1 == 2 or $j3 + $flag2 == 2) {
											
											if($j5 + $flag1 == 2) { #push(@{$candidate{$rl}->{$cut}->{$top_1_pare_chr}},   [$cluster_start, $cluster_end,$motif_prime,$top_1_pare_strand]);
												$motif_prime = 'ups';
												candidate_judge($top_1_pr5, $motif_prime,  $cluster_start, $cluster_end, $rl, $cut);
											}
											elsif ($j3 + $flag2 == 2) {
												$motif_prime = 'dws';
												candidate_judge($top_1_pr3, $motif_prime,  $cluster_start, $cluster_end, $rl, $cut);
											}
											else {
												# at least one of the four return value is equal to 0
												die "weird 4th\n";
											}
										}
										else {
											# neither top3 nor phased
											push(@{$miscandidate{$rl}->{$cut}->{$chr}},[$cluster_start, $cluster_end,'tpon','tpon']); # either top3 or phased or not
										}
									}
							}
							else{
								#seems that there is no pare located in the two regions, perhaps this is not a PHAS locus
								push(@{$miscandidate{$rl}->{$cut}->{$chr}},[$cluster_start, $cluster_end,'lt6','lt6']);
							}						
						}
						else {
							push(@{$miscandidate{$rl}->{$cut}->{$chr}},[$cluster_start, $cluster_end,'tlt10','tlt10']); #total less than 10
						}
					}
				}
			}
	}
		
	# ----------------output PARE corrected phasiRNA cluster boundaries ---------------------------------
	#print Dumper(\%candidate);
	
	open(VALIDCLUSTER,">$prefix.valid.cluster.YES.by.PARE.txt");
	foreach my $rl (sort {$a<=>$b} keys %candidate) { 
		foreach my $cut (sort {$a <=> $b} keys %{$candidate{$rl}}) {
			foreach my $chr (sort {$a <=> $b} keys %{$candidate{$rl}->{$cut}}) {
				my @sorted_cluster_ref  = sort_cluster(@{$candidate{$rl}->{$cut}->{$chr}});
				foreach my $cluster (@sorted_cluster_ref){
					$pare_summary{'valid'}->{$cut}->{$rl}++;
					my ($cluster_start, $cluster_end,$motif_prime, $strand) = @{$cluster};
					print VALIDCLUSTER join("\t",$rl,$cut,$chr,$cluster_start, $cluster_end, $motif_prime, $strand);
					print VALIDCLUSTER "\n";
				}
			}print  VALIDCLUSTER "\n";
		}
	}	
	close 	VALIDCLUSTER;
	
	open(INVALIDCLUSTER,">$prefix.valid.cluster.NO.by.PARE.txt");
	foreach my $rl (sort {$a<=>$b} keys %miscandidate) { #push(@{$miscandidate{$rl}->{$cut}->{$chr}},[$cluster_start, $cluster_end]);
		foreach my $cut (sort {$a <=> $b} keys %{$miscandidate{$rl}}) {
			foreach my $chr (sort {$a <=> $b} keys %{$miscandidate{$rl}->{$cut}}) {
				my @sorted_cluster_ref  = sort_cluster(@{$miscandidate{$rl}->{$cut}->{$chr}});
				foreach my $cluster (@sorted_cluster_ref){
					$pare_summary{'invalid'}->{$cut}->{$rl}++;
					my ($cluster_start, $cluster_end, @other) = @{$cluster};
					eval {
						print INVALIDCLUSTER join("\t",$rl,$cut,$chr,$cluster_start, $cluster_end, @other);
						print INVALIDCLUSTER "\n";
					};
				}
			}print INVALIDCLUSTER "\n";
		}
	}	
	close 	INVALIDCLUSTER;
	
	# -------------------------------------- output the pare analysis summary ----------------------------
	open(PARE_VALID,">$prefix.summary.PARE.v$version.txt");
	foreach my $valid (sort {$b cmp $a } keys %pare_summary)  {
		print PARE_VALID "#The below is $valid clusters by PARE data\n";
		foreach my $cut (sort {$b <=> $a} keys %{$pare_summary{$valid}}) {
			print PARE_VALID "$cut";
			foreach my $rl (sort {$a <=> $b} keys %{$pare_summary{$valid}->{$cut}}) {
				my $amount = $pare_summary{$valid}->{$cut}->{$rl};
				print PARE_VALID "\t$rl\t$amount";
			}print PARE_VALID "\n";
		}print PARE_VALID "\n";
	}
	close PARE_VALID;
}


# -------------------------------- remove the temporary files  ------------------------------------------
foreach my $file (@temp_files) {
	unlink "$file";
}

Ptime("Gorgeous, you are done...");

# _____   _   _   __   _   _____   _____   _   _____   __   _   _____  
#|  ___| | | | | |  \ | | /  ___| |_   _| | | /  _  \ |  \ | | /  ___/ 
#| |__   | | | | |   \| | | |       | |   | | | | | | |   \| | | |___  
#|  __|  | | | | | |\   | | |       | |   | | | | | | | |\   | \___  \ 
#| |     | |_| | | | \  | | |___    | |   | | | |_| | | | \  |  ___| | 
#|_|     \_____/ |_|  \_| \_____|   |_|   |_| \_____/ |_|  \_| /_____/ 
#
# --------------------------------functions---------------------------------------------------------------
sub bowtie_process { # won't change the original strand information
	my ($strand,$pos,$seq) = @_;
	
	if ($strand eq '+') {
		$pos += 1;
	}

	elsif ($strand eq '-') {
		$pos = $pos + length($seq);
		$seq = reverse($seq);
		$seq =~ tr/ATCGN/TAGCN/;
	}
	else {die "something is wrong to check the bowtie output";}
	
	return ($pos,$seq);
}


sub Ptime{
	my $time = localtime;
	my ($msg)= @_;
	print "$time: $msg\n";
}

#-
sub n_calc {
	my ($chr,$str,$cor,$rl) = @_;
	my $n = 0;
	
	if ($str eq '+') {
		my $ss = $cor;
		my $ee = $cor + 11*$rl -1; #there would be total 11 position are in register to the cleavage site whatever size on any strand

		for (my $i = $ss; $i<= $ee ; $i++) { #count all the top strand siRNA
			if ($sh2{$chr}->{$i}->{'+'}) {
				$n++;
			}
		}
		for (my $i = $ss - 2;$i<= $ee - 2 ; $i++) {
			if ($sh2{$chr}->{$i}->{'-'}) {
				$n++;
			}
		}
	}
	elsif ($str eq '-') {
		my $ss = $cor - 11*$rl + 1; 
                
		if ($ss<=0) {
			$ss = 1;
		}

		my $ee = $cor;

		for (my $i = $ss;$i <= $ee ; $i++) {
			if ($sh2{$chr}->{$i}->{'+'}) {
				$n++;
			}
		}       
		for (my $i = $ss + 2; $i<=$ee + 2 ;$i++) {
			if ($sh2{$chr}->{$i}->{'-'}) {
				$n++;
			}
		}
	}

	else {
		die "wrong 1234\n";
	}
	return $n;
}

sub k_calc {
	my ($chr,$str,$cor,$rl) = @_;
	my $k = 0;
	
        if ($str eq '+') {
                my $ss = $cor;
                my $ee = $cor + 11*$rl -1;
                

                for (my $i = $ss; $i<= $ee ; $i += $rl) {#here the sRNA located in the top strand
                        if ($sh2{$chr}->{$i}->{'+'}) {
                                $k++;
                        }
                }
                
                for (my $i = $ss + $rl - 3; $i <= $ee - 2; $i += $rl) {
                        if ($sh2{$chr}->{$i}->{'-'}) {
                                $k++;
                        }
                }
        }
        elsif ($str eq '-') {
				my $ss = $cor;
				my $ee = $cor - 11*$rl + 1;
				for (my $i = $ss; $i >= $ee; $i -= $rl) {#here the sRNA located in the bottom strand
                        if ($sh2{$chr}->{$i}->{'-'}) {
                                $k++;
                        }
                }
                #for (my $i =  $ss + 2; $i <= $ee + 2; $i += $rl) {#here the sRNA located in the top strand
                 for (my $i =  $ss + 3 - $rl; $i >= $ee + 2; $i -= $rl) {#here the sRNA located in the top strand
						if ($sh2{$chr}->{$i}->{'+'}) {
                                $k++;
						}
                }
        }
        else {
                die "wrong 2345\n";
        }
        return $k;
}

sub return_phasiRNA {
	my ($chr,$cor,$str,$rl) = @_;
	my @valid_phased = ();
	
	if ($str eq '+') {
		for (my $i = $cor ; $i<= $cor + $rl*10; $i += $rl) {
			if ($sh2{$chr}->{$i}->{'+'}) {
				my @info_ref = @{$sh2{$chr}->{$i}->{'+'}};
				foreach my $pha (@info_ref){
					my ($seqid,$seq,$abun) = @{$pha};
					push(@valid_phased,[$chr,$i,'+', $seqid,$seq,$abun]);
				}
			}
		}
		
		for (my $i = $cor + $rl - 3 ; $i <= $cor+ $rl*11 - 3; $i += $rl ) {
			if ($sh2{$chr}->{$i}->{'-'}) {
				my @info_ref = @{$sh2{$chr}->{$i}->{'-'}};
				foreach my $pha (@info_ref){
					my ($seqid,$seq,$abun) = @{$pha};
					push(@valid_phased,[$chr,$i,'-', $seqid,$seq,$abun]);
				}
			}
		}
	}

	elsif ($str eq '-') {
		for (my $i = $cor ; $i >=  $cor - $rl*10; $i -= $rl) {
			if ($sh2{$chr}->{$i}->{'-'}) {
				my @info_ref = @{$sh2{$chr}->{$i}->{'-'}};
				foreach my $pha (@info_ref){
					my ($seqid,$seq,$abun) = @{$pha};
					push(@valid_phased,[$chr,$i,'-', $seqid,$seq,$abun]);
				}
			}
		}
		for (my $i = $cor - $rl + 3 ; $i >= $cor - $rl*11 + 3; $i -= $rl ) {
			if ($sh2{$chr}->{$i}->{'+'}) {
				my @info_ref = @{$sh2{$chr}->{$i}->{'+'}};
				foreach my $pha (@info_ref){
					my ($seqid,$seq,$abun) = @{$pha};
					push(@valid_phased,[$chr,$i,'+', $seqid,$seq,$abun]);
				}
			}
		}
	}
	
	else {
		die "you strand is not like + and -\n";
	}
	return @valid_phased;
}
#push(@{$sh2{$chr}->{$pos}->{$strand}},[$seqid,$seq,$abun]);

sub pha_score_replace {
	my ($phasiRNA_ref,$cor_best,$rl,$p,$n,$k) = @_;
	my @phasiRNA_ref = @{$phasiRNA_ref};
	
	foreach my $pha_ref(@phasiRNA_ref) {
		my ($chr,$cor,$str,$seqid,$seq,$abun) = @{$pha_ref};
		if ($pre_scored_sRNA{$rl}->{$chr}->{$cor}->{$str}->{$seq}) {
			
			my $p_old = $pre_scored_sRNA{$rl}->{$chr}->{$cor}->{$str}->{$seq}->[-1];
			
			if ($p_old>$p) {
				$pre_scored_sRNA{$rl}->{$chr}->{$cor}->{$str}->{$seq}->[-1] = $p;
				$pre_scored_sRNA{$rl}->{$chr}->{$cor}->{$str}->{$seq}->[-4] = $n;
				$pre_scored_sRNA{$rl}->{$chr}->{$cor}->{$str}->{$seq}->[-3] = $k;
				$pre_scored_sRNA{$rl}->{$chr}->{$cor}->{$str}->{$seq}->[-2] = $cor_best;
			}
		}else {
			$pre_scored_sRNA{$rl}->{$chr}->{$cor}->{$str}->{$seq} = [$seqid,$abun,$n,$k,$cor_best,$p];
		}
	}
}

sub phasing_analysis_mul {
	my ($chr,$rl) = @_;
	open(OUT, ">$prefix.$options->{o}\_all_sRNA_$rl\_out.chr$chr.txt"); # outfile of all small RNA clusters
	open(out2,">$prefix.$options->{s}\_p$p_cutoff\_sRNA_$rl\_out.chr$chr.txt"); # outfile of small RNA clusters with p<0.001
	
	foreach my $cor (sort {$a <=> $b} keys %{$sh2{$chr}}) {
		foreach my $str (sort {$a cmp $b} keys %{$sh2{$chr}->{$cor}}) {
					# 
					if ($str eq '+') {
						
							my $n = n_calc($chr,$str,$cor,$rl);
							my $k = k_calc($chr,$str,$cor,$rl);			
							my $p = hm_chen($n,$k,$rl);#code;
							
							#----------------calc the p value -------------------------------
							my @phasiRNA = return_phasiRNA($chr,$cor,$str,$rl);
							pha_score_replace(\@phasiRNA,$cor,$rl,$p,$n,$k);
					}
					elsif ($str eq '-') {

							my $n = n_calc($chr,$str,$cor,$rl);
							my $k = k_calc($chr,$str,$cor,$rl);		
							my $p = hm_chen($n,$k,$rl);#code;
		
							#----------------calc the p value -------------------------------							
							my @phasiRNA = return_phasiRNA($chr,$cor,$str,$rl);
							pha_score_replace(\@phasiRNA,$cor,$rl,$p,$n,$k);
					}
					else {
						die "wrong, no strand found\n";
					}
		}
	}

	#$pre_scored_sRNA{$rl}->{$chr}->{$cor}->{$str}->{$seq} = [$seqid,$seq,$abun,$n,$k,$cor_best,$p];
#	foreach my $chr  (sort {$a <=> $b} keys %{$pre_scored_sRNA{$rl}}) {
		foreach my $cor (sort {$a <=> $b} keys %{$pre_scored_sRNA{$rl}->{$chr}}) {
			foreach my $str (sort {$a cmp $b} keys %{$pre_scored_sRNA{$rl}->{$chr}->{$cor}}) {
				foreach my $seq (sort {$a cmp $b} keys %{$pre_scored_sRNA{$rl}->{$chr}->{$cor}->{$str}}) {
					my ($seqid,$abun,$n,$k,$cor_best,$p) = @{$pre_scored_sRNA{$rl}->{$chr}->{$cor}->{$str}->{$seq}};
					my $len = length($seq);
						print OUT join ("\t",$chr,$str,$cor,$seqid,$seq,$len,$abun,$n,$k,$cor_best,$p);
						print OUT "\n";
					                                           
						if ($p < $p_cutoff) {#select and output small RNA clusters with p<0.001
								print out2 join ("\t",$chr,$str,$cor,$seqid,$seq,$len,$abun,$n,$k,$cor_best,$p);
								print out2 "\n";
						}
				}
			}
		}
#	}
	close OUT;
	close out2;
}



sub generate_rand_filename {

     my $length_of_randomstring = shift;# the length of
     # the random string to generate

     my @chars=('a'..'z','A'..'Z','0'..'9','_');
     my $random_string;
     foreach (1..$length_of_randomstring)
     {
          # rand @chars will generate a random
          # number between 0 and scalar @chars
          $random_string .= $chars[rand @chars];
     }
     return $random_string; 

}

#v.2 2013-06-24, with proper modification to get it adapted to different register length
sub hm_chen { #v.2
	my ($n,$k,$rl) = @_;
	#$k = ($k>=22)?21:$k;
	
	my $p = 0;
	my $pr;

	for (my $w = $k; $w <= 22; $w++) {
	
		my $c=1;
		my $rr=1;
		my $rw=1;
        
		for (my $j = 0; $j <= $w-1; $j++){
			#$c=$c*($n-$j)/($j+1);
						$c=$c*($n-$j)/($j+1);
		}
        
		for (my $x = 0; $x <= $w-1; $x++){
			#$rr=$rr*(21-$x)/(461-$x);
			eval
					{
						$rr=$rr*(22 - $x)/($rl*11*2  - $x);
					}
		}
        
		for (my $y = 0; $y <= $n-$w-1; $y++){
			#$rw=$rw*(440-$y)/(461-$w-$y);			
			eval
					{
						$rw=$rw*($rl*22 - 22 - $y)/($rl*11*2 - $w - $y);
					}
						
		}
        
		$pr = $c*$rr*$rw;
		$p += $pr;
	}
	return $p;
}

sub merge_files {
	my ($file_ref,$output) = @_;
	my %file = ();
	foreach my $file(@{$file_ref}) {
		my ($chr) = $file =~/chr(\d+)?\.txt/;
		$file{$chr}= $file;
		#print "$chr\t$file\n";
	}
	my @files = map{$file{$_}} sort {$a<=>$b} keys %file;
	
	#print "cat @files > $output\n";
	system("cat @files > $output");
	
	foreach my $file(@files){
		unlink "$file";
	}
}

sub minimal {
	my @array = @_;
	my @array_sorted = map{$array[$_]} sort {$array[$a]->[0]<=>$array[$b]->[0]} 0..$#array;
	return $array_sorted[0];
}

# this function is mainly used to check the proportion of the interested small RNA abundance and average hits
sub cluster_abundance_check {
	my ($hash_ref, $rl, $ta_cf,$cut) = @_;
	my %cluster = %{$hash_ref};
	my $i = 0;
	my $phasedNumber = 0;
	
	foreach my $cluster_id (sort {$a <=> $b} keys %cluster) {
		my $total = 0;
		my $total_interested = 0;
		my $abun_HNA = 0;
		
		my @cluster_ref =  @{$cluster{$cluster_id}};
		my $j = 0;
		my $t_hits = 0;

		foreach my $record_ref (@cluster_ref) {
			$j++;
			my ($chr,$str,$pos,$seq_id,$seq,$len,$abun,$n,$k,$hts,$best_pos,$pvalue) = @{$record_ref};
			  $total += $abun;
			  $total_interested += $abun if ($len == $rl);
			  $abun_HNA += $abun/$hits{$seq_id};			
			  
			  my ($hits_value) = $hts =~ /=(\d+)/;
			  $t_hits += $hits_value;
		}
		# ----check the proportion of the interested small RNA-----
		if (100*$total_interested/$total >= $ta_cf and $t_hits/$j <= $avg_hit) {
			$i++;
			$sum{$cut}->{$rl}->{'phaseClus'}++;
			$phasedNumber += $abun_HNA;
			
			$qualified_cluster{$cut}->{$rl}->{$i} = $cluster{$cluster_id};
			
		}
	}
	$sum{$cut}->{$rl}->{'phaseNumb'} = int($phasedNumber*10_000_000/$total_r_abundance);
}

sub boundary {
	my ($cluster_ref,$rl) = @_;
	
	my @cluster = @{$cluster_ref};
	
	my $start = prime53($cluster[0][2],5,$cluster[0][1],$rl);
	my $end = prime53($cluster[-1][2],3,$cluster[-1][1],$rl);
	my $chr = $cluster[0][0];
	
	return($start,$end,$chr);
	
}

sub prime53 {
	my ($pos,$prime,$str,$rl) = @_;
	if ($prime == 5) {
		if ($str eq '+') {
			
		}
		elsif ($str eq '-') {
			$pos  = $pos - ($rl - 3);
		}
		else {
			die "writing 5' start calc\n";
		}
	}
	elsif ($prime == 3) {
		if ($str eq '+') {
			$pos = $pos + $rl - 3;
		}
		elsif ($str eq '-') {
			
		}
		else {
			die "writing 3' start calc\n";
		}		
	}
	else {
		die "writing edge start calc\n";
	}
	return $pos;
}

sub pare_seq_combine {
	my ($seq1,$seq2) = @_;
	
		if (defined $seq1) {
			if (length($seq1) >= length($seq2)) {
				return $seq1;
			}
			else {
				return $seq2;
			}
		}
		else {
			return $seq2;				
		}
}
#
sub pare_abun_total {
	my ($chr,$start,$end) = @_;
	
	my $region_total = 0;  # regional_total can be zero if there is no pare located in the specified range
	my $seq; 
	#call the varant %PARE_alignment
	foreach my $pos ($start..$end) {
		foreach my $strand ("+","-") {
			if ($PARE_alignment{$chr}->{$pos}->{$strand}->{'abun'}) {
				$region_total += $PARE_alignment{$chr}->{$pos}->{$strand}->{'abun'};
				#$seq = $PARE_alignment{$chr}->{$pos}->{$strand}->{'seq'};
			}
		}
	}
	return $region_total;
}

# similar to the above maxmial_PARE_scan, but it return top 3 pare information in a -$rl*4 and + $rl*4 range
sub top_pare {
	my ($chr,    $pos1,    $pos2, $return_no) = @_;
	my @tmp = ();
	
	foreach my $pos ($pos1..$pos2) {
		foreach my $strand ("+","-") {
			if ($PARE_alignment{$chr}->{$pos}->{$strand}->{'seq'}) {
				my $seq   = $PARE_alignment{$chr}->{$pos}->{$strand}->{'seq'};
				my $abun  = $PARE_alignment{$chr}->{$pos}->{$strand}->{'abun'};
				push(@tmp,[$abun,$chr,$pos,$strand,$seq]);
			}
		}
	}
	# order as desc
	my @sorted = map{$tmp[$_]} sort{$tmp[$b]->[0] <=> $tmp[$a]->[0]}  0..$#tmp;
    
	my @list;
	foreach my $idex (0..($return_no - 1)) {
		if ($sorted[$idex][1]) {
			push(@list,$sorted[$idex]);
		}
	}
	return(@list);
}


#
sub candidate_judge {
	my ($top_1_pr, $motif_prime, $cluster_start, $cluster_end,$rl,$cut) = @_;
	my ($top_1_pare_abun,$top_1_pare_chr,$top_1_pare_pos,$top_1_pare_strand,$top_1_pare_seq) = @{$top_1_pr};
	
		if ($motif_prime eq 'ups') {
			if ($top_1_pare_strand eq '+') {
				$cluster_start = $top_1_pare_pos;
			}
			elsif ($top_1_pare_strand eq '-') {
				$cluster_start = $top_1_pare_pos + 1;  #?need further check
			}
			else {
				die "wrong code 1231\n";
			}
		}
		elsif ($motif_prime eq 'dws') {
			if ($top_1_pare_strand eq '+') {
				#two hits model
				$cluster_end = $top_1_pare_pos - 1;  #?need further check
			}
			elsif ($top_1_pare_strand eq '-') {
				$cluster_end = $top_1_pare_pos;
			}
			else {
				die "wrong code 1232\n";
			}
		}
		else {
			die "wrong 3342\n";
		}
		
		# %candidate
		push(@{$candidate{$rl}->{$cut}->{$top_1_pare_chr}},   [$cluster_start, $cluster_end, $motif_prime, $top_1_pare_strand]);
}  
 
# 
sub  sort_cluster {
	my (@cluster_ref) = @_;
	my @sorted = map{$cluster_ref[$_]} sort{$cluster_ref[$a]->[0] <=> $cluster_ref[$b]->[0]}  0..$#cluster_ref;
	return @sorted;
}


#
sub top3_phase_judge {
	my ($top_3_pare, $top_1_prime_ref, $motif_prime, $old_cluster_start, $old_cluster_end, $rl) = @_;
	my $top  = 0;
	my $phased_or_not = 0;
	
	if (defined $top_1_prime_ref->[2]) {
	
		my ($abun,$chr,$pare_pos,$pare_strand,$seq) = @{$top_1_prime_ref};
	
	
		my $top_1_candidate = join("_",$abun,$chr,$pare_pos,$pare_strand,$seq);
		
		#print "$motif_prime : $old_cluster_start, $old_cluster_end, pare pos = $pare_pos, pare abun = $abun\n";
		
		foreach my $ref (@{$top_3_pare}) {
			
				my ($abun,$chr,$pos,$strand,$seq) = @{$ref};
				my $tmp = join("_",$abun,$chr,$pos,$strand,$seq);
				if ($tmp eq $top_1_candidate) {
					$top = 1; #match to the top1 
				}
			
			
			
			
		}
		
		
		# phased or not judge
		if    ($motif_prime eq 'ups') {
			
			if ($pare_strand eq '+') {
				
				if (abs($old_cluster_start - $pare_pos)%$rl == 0){
					$phased_or_not = 1;
				}
				else {
					$phased_or_not = 0;
				}
			}
			# example like TAS3 in rice
			elsif ($pare_strand eq '-') {
				
				if (abs($pare_pos + 3 - $old_cluster_start)%$rl == 0) {
					$phased_or_not = 1;
				}
				else {
					$phased_or_not = 0;
				}
			}
			else {
				die "wrong 34234\n";
			}		
		}
		elsif ($motif_prime eq 'dws') {
			
			if ($pare_strand eq '-') {
				
				
				if (abs($old_cluster_end - $pare_pos)%$rl == 0) {
					$phased_or_not = 1;	
				}
				else {
					$phased_or_not = 0;
				}
			}
			# example like TAS3 in rice
			elsif ($pare_strand eq '+') {
			
				
				#print "$pare_pos,$old_cluster_end","\t",abs($pare_pos - $old_cluster_end - 3)%$rl;
				if (abs($pare_pos - $old_cluster_end - 3)%$rl == 0 ) {
					$phased_or_not = 1;
				}
				else{
					$phased_or_not = 0;
				}
			}
			else {
				die "wrong 34234\n";
			}
		}
		
	}
	
	
	return ($top,  $phased_or_not);
}










