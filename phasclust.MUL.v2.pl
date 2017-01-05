#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Basename;

######################################################################
#  PhasiRNA Prediction Pipeline (PPP)
#  Copyright 2013 University of Delaware
#  Pingchuan Li: li@dbi.udel.edu
#  Blake Meyers: meyers@dbi.udel.edu
#
#  Department of Plant Soil Science
#  Delaware Biotechnology Insititute
#  University of Delaware
#  Newark DE, 19711
# 
#  This program is a FREE software; you can redistribute it and/or
#  modify as your wishes
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
######################################################################

my $version  = 1.0;
Ptime("start...");

# make sure Perl is over 5.12
no 5.12;
my %sum = ();
GetOptions(my $options = {},
              "-i=s","-f=s","-d=s","-v=i","-o=s","-rl=s","-s=s","-p=s","-cpu=i","-px=s","-n=i","-g=i","-t=i","-ht=i","-q=s","-k=i"#"-z=s",
);
########### USAGE ##############
my $USAGE = <<USAGE;
Scripts: cluster analysis for phased siRNA
The sRNA reads were mapped by bowtie by default.

Version: $version, Written by Pingchuan Li \@ Blake Meyers Lab
Usage and Arguments:

    Mandatory:
        -i        small RNA input file name
        -f        two file formats are supported,[t] for tagcount, [f] for fasta.
        -d        indexed genome by bowtie, indicating the prefix of the genome
        -px       prefix for each of the file, used to distinguish for eachother
        -rl       register len(rl), such as 21 or 24, separated by comma: 21,24 etc

    optional:
        -k        bowtie alignment hits report, default = all
        -q        degradome or PARE data in ONLY tagcount format: seq   abundance
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
$options->{a} = ($options->{a})?$options->{a}:'m';
my $pare_data = ($options->{q})?$options->{q}:0;
my $cpu = ($options->{cpu})?$options->{cpu}:2;
my $mm = (defined $options->{m})?$options->{m}:0;
my $input = $options->{i};
my $k = ($options->{k})?"-k $options->{k}" :"-a";  #bowtie options for hits reports limitation.

$options->{o} = ($options->{o})?$options->{o}:"output";
$options->{s} = ($options->{s})?$options->{s}:'score';

my %sh = ();#define the Small RNA Hash data. use the unique combination of information as the key.
my %sh2 = ();# a hash dataset for the qualified clusters
my %hits = (); # this will record the hits number for each of the small RNAs
my %qualified_cluster = ();
my %boundary = (); # this will record the boundary of each predicted phasiRNA cluster
my %pre_scored_sRNA = (); #this will pre-saving the score of each siRNA in a sliding window

# the below two hash will record validate cluster by pare data and unvalided cluster if possible;
my %candidate = ();
my %miscandidate = ();
push(@temp_files,$prefix.".".$bowtie_output);
# -----------tag count processing if the input is tagcount---------------
my $t_i = 0;
if ($options->{f} =~ /t/i) {
	my $fasta_file = generate_rand_filename(11);
	push(@temp_files,$fasta_file);
	
	open(IN, $input);
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
my @parameters = join(" ", "-f",  "-a -v $mm","-m 45", "-p $cpu", "$options->{d}", "$input", " >$prefix.$bowtie_output");

my @bowtie_out = system("bowtie @parameters");
Ptime("aligning is done");

#put all the signature and their information into the HASH %sh
#--------bowtie out put pre processing, organized as the position, deposited into the HASH, MYSQL processing will be better with lowere RAM requirement-------------------
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

	# calc the total mapped reads abundance
	$total_r_abundance += $abun unless ($hits{$seqid});
	$hits{$seqid}++;

	# in order to improve the sensitivity, run the below filter.
	next if ($abun <= $noise);
        
	($pos,$seq) = bowtie_process($strand,$pos,$seq);
	
	#if the same position of the same chromosome contains multiple small RNA in different length or strand, they will be put together
        push(@{$sh{$target}->{$pos}->{$strand}},[$seqid,$seq,$abun]);
}
Ptime("finished the loading of sRNA alignment data.");
# ------------ end of preprocessing ---------------
close DATA;

my @pre_splitted_cluster = ();   #@pre_splitted_cluster = ( [[p1],[p2.1],[p2.2],[p3]], [[p1]],  [ [p1],[p2] ] );
#-----------------------sorting, divided  the cluster prior to the calculation ------------------------------
foreach my $target (sort {$a cmp $b} keys %sh) {
	my $i = 0;
	my $temp_pos = 0;
	foreach my $pos (sort {$a <=> $b} keys %{$sh{$target}}) {
		foreach my $strand ( keys %{$sh{$target}->{$pos}}) {
			my @sRNA_ref = @{$sh{$target}->{$pos}->{$strand}};
			foreach my $sRNA_ref (@sRNA_ref) {
				my ($seqid,$seq,$abun) = @{$sRNA_ref};
				$i++;
				
				#initialization
				if ($i == 1) { #the first sRNA in a chr
				    $temp_pos = $pos;
				    push(@pre_splitted_cluster, [[$target,$pos,$strand,$seqid,$seq,$abun]]);
				}
				
				else {
					# judge the distance for two sRNA
					# if shorter than the pre-defined distance
					if ($pos - 400 <= $temp_pos) {
						
						push(@{$pre_splitted_cluster[-1]},[$target,$pos,$strand,$seqid,$seq,$abun]);					
						$temp_pos = $pos;
						
					} else {
						
						push(@pre_splitted_cluster, [[$target,$pos,$strand,$seqid,$seq,$abun]]);
						$temp_pos = $pos;
						
					}										
				}
			}
		}
	}
}

Ptime("finished the sorting of data.");

# push(@{$sh2{$chr}->{$pos}->{$strand}},[$seqid,$seq,$abun]);
my $total_cluster = 0;
my @qualified_cluster;
foreach my $cluster_ref (@pre_splitted_cluster) {
	my @record_ref = @{$cluster_ref};
	if ($#record_ref <= 1) {
		
	}else {
		$total_cluster++;
		foreach my $sRNA_ref (@record_ref) {
			my ($target,$pos,$strand,$seqid,$seq,$abun) = @{$sRNA_ref};
			
			#the sh2 is an image of qualified_cluster in Hash format
			push(@{$sh2{$target}->{$pos}->{$strand}},[$seqid,$seq,$abun]);
		}
		push(@qualified_cluster, $cluster_ref);
	}
}

Ptime("total $total_cluster presplitted sRNA clusters");

# -------------------------------- scoreing each of the cluster parallelly --------------------------------------------------
foreach my $rl (sort {$a <=> $b} @rl) {
	
	my @file_list = ();# save the temporary filenames
	# the below loops will spend a lot of time,  if this part can be optimized by parallel, then the efficiency woulbe be greatly improved-
	my $i = 0;
	my $j = 1;
	
	while (my ($cluster_idx, $cluster_ref) = each @qualified_cluster) { # the each for @array only works on the perl over 5.12
		$i++;
		if ($i/500 == $j ) {
			Ptime("processed $i clusters");
			$j++;
		}
		my $filename = phasing_analysis_v2($cluster_ref, $rl, $cluster_idx);		
		#push(@file_list,$filename);
	}
	# ----------------------------------end of score loop ---------------------------------------------------------------

	open(OUT,  ">$prefix.$options->{o}\_all_sRNA_$rl\_out.txt"); # outfile of all small RNA clusters
	open(out2, ">$prefix.$options->{s}\_p$p_cutoff\_sRNA_$rl\_out.txt") or die "cant create the scored file\n";

	#$pre_scored_sRNA{$rl}->{$chr}->{$cor}->{$str}->{$seq} = [$seqid,$seq,$abun,$n,$k,$cor_best,$p];
	foreach my $target  (sort {$a cmp $b} keys %{$pre_scored_sRNA{$rl}}) {
		foreach my $cor (sort {$a <=> $b} keys %{$pre_scored_sRNA{$rl}->{$target}}) {
			foreach my $str (sort {$a cmp $b} keys %{$pre_scored_sRNA{$rl}->{$target}->{$cor}}) {
				foreach my $seq (sort {$a cmp $b} keys %{$pre_scored_sRNA{$rl}->{$target}->{$cor}->{$str}}) {
					my ($seqid,$abun,$n,$k,$cor_best,$p) = @{$pre_scored_sRNA{$rl}->{$target}->{$cor}->{$str}->{$seq}};
					my $len = length($seq);
						print OUT join ("\t",$target,$str,$cor,$seqid,$seq,$len,$abun,$n,$k,$cor_best,$p);
						print OUT "\n";
						                                          
						if ($p < $p_cutoff) {#select and output small RNA clusters with p<0.001
								print out2 join ("\t",$target,$str,$cor,$seqid,$seq,$len,$abun,$n,$k,$cor_best,$p);
								print out2 "\n";
						}
				}
			}
		}
	}
	close OUT;
	close out2;
	
	Ptime("The scoring is done...");

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
				
				if ($temp_chr eq $array[0]) {
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
			print OUT_CLU ">cluster = $i | target = $chr and pos between $start and $end\n";
			print CLUSTER_LITE "$cut|$rl|target = $chr:$start..$end\n";
			
				foreach my $ref (@{$cluster_ref}) {
					my ($chr,$str,$pos,$seq_id,$seq,$len,$abun,$n,$k,$hts,$best_pos,$pvalue) = @{$ref};
					print  OUT_CLU join("\t","o",$chr,$str,$pos,$seq_id,$seq,$len,$abun,$n,$k,$hts,$best_pos,$pvalue);
					print  OUT_CLU "\n";
				}
		}print CLUSTER_LITE "\n";
		close OUT_CLU;
	}print CLUSTER_LITE "\n";
}
close CLUSTER_LITE ;


# -------------------------------- remove the temporary files  ------------------------------------------
foreach my $file (@temp_files) {
	unlink "$file";
}

Ptime("Gorgeous, you are done...");

#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------sub functions-----------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
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

sub phasing_analysis_v2 {
		my ($cluster_ref, $rl, $cluster_idx) = @_;
		my @record_ref = @{$cluster_ref};
		my $sRNA_last_index = $#record_ref;
		
		# the below array will deposited the qualified cluster
		my @cluster_verified = ();
			
		# process the cluster
		while (my ($curr_idx, $sRNA_ref) = each @record_ref) { # a cluster processing, loop for sRNA within the cluster on both strands
			
			# acquire the meta information for every sRNA in the cluster
			my ($target, $pos, $strand, $seqid, $seq, $abun) = @{$sRNA_ref};
			my $seq_len = length($seq);
			
			# the score array will contain all the possible scores
			my @score_array = ();
			
			# to check whether the other sRNA were phased with above sRNA in a 11 windows from upstream to downstream
			if ($strand eq '+') {
					my $n = n_calc($target,$strand,$pos,$rl);
					my $k = k_calc($target,$strand,$pos,$rl);			
					my $p = hm_chen($n,$k,$rl);#code;
					
				#----------------calc the p value -------------------------------
					my @phasiRNA = return_phasiRNA($target,$pos,$strand,$rl);
					pha_score_replace(\@phasiRNA,$pos,$rl,$p,$n,$k);
					
			}
			elsif ($strand eq '-') {
					my $n = n_calc($target,$strand,$pos,$rl);
					my $k = k_calc($target,$strand,$pos,$rl);		
					my $p = hm_chen($n,$k,$rl);#code;
	
				#----------------calc the p value -------------------------------							
					my @phasiRNA = return_phasiRNA($target,$pos,$strand,$rl);
					pha_score_replace(\@phasiRNA,$pos,$rl,$p,$n,$k);
			}
			else {
				die "wrong, no strand found\n";
			}
			#push(@cluster_verified,[$target,$strand,$pos,$seqid,$seq,$seq_len,$abun,$n,$k,$best_win_pos,$p]);
		}
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
#
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
			die "writing 5 start calc\n";
		}
	}
	elsif ($prime == 3) {
		if ($str eq '+') {
			$pos = $pos + $rl - 3;
		}
		elsif ($str eq '-') {
			
		}
		else {
			die "writing 3 start calc\n";
		}		
	}
	else {
		die "writing edge start calc\n";
	}
	return $pos;
}

# sub pare_seq_combine {
# 	my ($seq1,$seq2) = @_;
	
# 		if (defined $seq1) {
# 			if (length($seq1) >= length($seq2)) {
# 				return $seq1;
# 			}
# 			else {
# 				return $seq2;
# 			}
# 		}
# 		else {
# 			return $seq2;				
# 		}
# }
#
# sub pare_abun_total {
# 	my ($chr,$start,$end) = @_;
	
# 	my $region_total = 0;  # regional_total can be zero if there is no pare located in the specified range
# 	my $seq; 
# 	#call the varant %PARE_alignment
# 	foreach my $pos ($start..$end) {
# 		foreach my $strand ("+","-") {
# 			if ($PARE_alignment{$chr}->{$pos}->{$strand}->{'abun'}) {
# 				$region_total += $PARE_alignment{$chr}->{$pos}->{$strand}->{'abun'};
# 				#$seq = $PARE_alignment{$chr}->{$pos}->{$strand}->{'seq'};
# 			}
# 		}
# 	}
# 	return $region_total;
# }

# # similar to the above maxmial_PARE_scan, but it return top 3 pare information in a -$rl*4 and + $rl*4 range
# sub top_pare {
# 	my ($chr,    $pos1,    $pos2, $return_no) = @_;
# 	my @tmp = ();
	
# 	foreach my $pos ($pos1..$pos2) {
# 		foreach my $strand ("+","-") {
# 			if ($PARE_alignment{$chr}->{$pos}->{$strand}->{'seq'}) {
# 				my $seq   = $PARE_alignment{$chr}->{$pos}->{$strand}->{'seq'};
# 				my $abun  = $PARE_alignment{$chr}->{$pos}->{$strand}->{'abun'};
# 				push(@tmp,[$abun,$chr,$pos,$strand,$seq]);
# 			}
# 		}
# 	}
# 	# order as desc
# 	my @sorted = map{$tmp[$_]} sort{$tmp[$b]->[0] <=> $tmp[$a]->[0]}  0..$#tmp;
    
# 	my @list;
# 	foreach my $idex (0..($return_no - 1)) {
# 		if ($sorted[$idex][1]) {
# 			push(@list,$sorted[$idex]);
# 		}
# 	}
# 	return(@list);
# }


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
 
sub  sort_cluster {
	my (@cluster_ref) = @_;
	my @sorted = map{$cluster_ref[$_]} sort{$cluster_ref[$a]->[0] <=> $cluster_ref[$b]->[0]}  0..$#cluster_ref;
	return @sorted;
}

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



## Log Change
## Removed -k mode in bowtie and added -a, replaced -n mode with -v and added -m ceiling. 
## In case of transcriptome/scaffold version -m ceiling is kept high as these can match to multiple isoforms





