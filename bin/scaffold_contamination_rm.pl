#!/usr/bin/env perl

my $mod = "11/16/12 4:49 PM";
my $version = "0.2";
my $author = "Nick Youngblut";
#--------------------- version log ---------------------#
#
#
#-------------------------------------------------------#

### packages/perl_flags
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Temp;
use Pod::Usage;
use Set::IntervalTree;
use List::MoreUtils qw(pairwise);

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

### I/O
my ($verbose, $ref_file, $blast_file, $cov_cut, $seqID_cut);
my $keep_cut = 25;
my $len_cut = 500;
GetOptions(
	   "fasta=s" => \$ref_file,
	   "blast=s" => \$blast_file,
	   "keep=i" => \$keep_cut,				# drop contig if < $keep_cut is left (percentage length)
	   "length=i" => \$len_cut,				# minimum length of each contig written
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### Input error check
die " ERROR: provide a reference scaffold file (fasta format)\n" if ! $ref_file;
$ref_file = File::Spec->rel2abs($ref_file);
die $! if ! -e $ref_file;
die " ERROR: provide a blast_2_hsp_table file (.txt)\n" if ! $blast_file;
$blast_file = File::Spec->rel2abs($blast_file);
die $! if ! -e $blast_file;


### Routing main subroutines
my ($itree_r, $max_len_r) = load_parsed_blast($blast_file);					# loading blast file

#foreach my $q (keys %$itree_r){ print_interval_tree($$itree_r{$q}, $$max_len_r{$q}); } exit;		# for debugging

my $merged_hits_r = merge_hits($itree_r, $max_len_r);
my $fasta_r = load_fasta($ref_file);										# loading scaffolds
$fasta_r = scaffold_split($fasta_r, $merged_hits_r, $keep_cut, $len_cut);
write_fasta($fasta_r);


#----------------------Subroutines----------------------#
sub write_fasta{
# writing out fasta #
	my $fasta_r = shift;
	foreach my $seq (sort keys %$fasta_r){
		print join("\n", ">$seq", $$fasta_r{$seq}), "\n";
		}
	
	}

sub scaffold_split{
# removing sections of scaffold based on blast hits # 
	my ($fasta_r, $merged_hits_r, $keep_cut, $len_cut) = @_;
	
	#print Dumper $merged_hits_r; exit;
	
	# staring stats #
	my %filter_stats;
	$filter_stats{"split"} = 0;
	$filter_stats{"split_cnt"} = 0;
	$filter_stats{"length_cut_total"}  = 0;
	$filter_stats{"seq_rm"} = 0;
	
	$filter_stats{"seq_starting"} = scalar keys %$fasta_r;
	map{ $filter_stats{"length_starting_total"} += length $$fasta_r{$_} } keys %$fasta_r;
	
	# removing contaminating segments #
	my %filter_log;
	foreach my $q (keys %$merged_hits_r){
		die " ERROR: $q not found in fasta file!\n" if ! exists $$fasta_r{$q};
		
		# substringing sequence #
		my @frags;
		my $frag_len_total = 0;
		foreach my $frag (keys %{$$merged_hits_r{$q}}){
			push(@frags, substr($$fasta_r{$q}, $$merged_hits_r{$q}{$frag}{"start"}, 
				$$merged_hits_r{$q}{$frag}{"end"} - $$merged_hits_r{$q}{$frag}{"start"}) );
			$frag_len_total += length $frags[$#frags];
			}
		
		# delete original sequence #
		my $seq_len = length $$fasta_r{$q};
		delete $$fasta_r{$q};
		
		# filtering #
		if($frag_len_total / $seq_len * 100 < $keep_cut){
			$filter_stats{"seq_rm"}++;
			$filter_stats{"length_cut_total"} += $seq_len;
			next;
			};			# if total fragment length % is < starting 
		
		my $keep_cnt = 0;
		for my $i (0..$#frags){
			if (length $frags[$i] < $len_cut){
				$filter_stats{"length_cut_total"} += length $frags[$i];			# if fragment remaining is < lenght cutoff
				}
			else{
				$$fasta_r{ join("_", $q, $i+1, length $frags[$i]) } = $frags[$i];
				$keep_cnt++;
				}
			}

		if (! $keep_cnt){ $filter_stats{"seq_rm"}++;}
		else{ 
			$filter_stats{"split"}++; 
			$filter_stats{"split_cnt"} += $keep_cnt;
			}
		}
	write_filter_log(\%filter_stats);
		
		#print Dumper %filter_stats; exit;
		#print Dumper %$fasta_r; exit;		
	return $fasta_r;
	}

sub write_filter_log{
# writing filter log #
	my $log_r = shift;
	my $log_q_r = shift;

	# writing basic stats #
	print STDERR "### Sequences removed stats ###\n";
	print STDERR "Sequences starting: 		", $$log_r{"seq_starting"}, "\n";
	print STDERR "Sequences removed: 		", $$log_r{"seq_rm"}, "\n";
	print STDERR "Sequences split: 		", $$log_r{"split"}, "\n";
	print STDERR "Number of new fragments:	", $$log_r{"split_cnt"}, "\n";
	print STDERR "Sequences remaining: 		", $$log_r{"seq_starting"} - $$log_r{"seq_rm"} + $$log_r{"split_cnt"}, "\n";
	print STDERR "Percent remaining:		", 
		sprintf("%.1f", ( $$log_r{"seq_starting"} - $$log_r{"seq_rm"} + $$log_r{"split_cnt"} ) / $$log_r{"seq_starting"} * 100), 
		"%\n\n";
	
	print STDERR "### Length cut stats ###\n";
	print STDERR "Length starting:		", $$log_r{"length_starting_total"}, "\n";
	print STDERR "Length cut: 			", $$log_r{"length_cut_total"}, "\n";
	print STDERR "Length remaining:		", $$log_r{"length_starting_total"} - $$log_r{"length_cut_total"}, "\n";
	print STDERR "Percent remaining:		", 
		sprintf("%.1f", ($$log_r{"length_starting_total"} - $$log_r{"length_cut_total"}) / $$log_r{"length_starting_total"} * 100), 
		"%\n\n";
	
	return;
	
	# writing log #
	open OUT, ">scaffold_contam.log" or die $!;
	print OUT join("\t", qw/scaffold_id length_starting length_cut length_remaining kept_remove/), "\n";
	foreach my $q (sort keys %$log_q_r){
		print OUT join("\t", $q, $$log_q_r{$q}{"length_starting"}, $$log_q_r{$q}{"length_cut"}, 
			$$log_q_r{$q}{"length_remaining"}, $$log_q_r{$q}{"kept_remove"}), "\n";
		}
	close OUT;
	print STDERR "Log file written: scaffold_contam.log\n";
	}

sub merge_hits{
# merging overlapping hits in the same query #
# keeping locations not hit w/ contamination #
	my ($itree_r, $max_len_r) = @_;
	
	my %merged_hits;
	foreach my $q (keys %$itree_r){
		# finding all hit intervals # 
		my $cnt = 0;
		for my $i (0..$$max_len_r{$q}){
			my $res = $$itree_r{$q}->fetch($i, $i);
			if(@$res){
				$cnt++;
				}
			else{
				$merged_hits{$q}{$cnt}{"start"} = $i if ! exists $merged_hits{$q}{$cnt}{"start"};
				$merged_hits{$q}{$cnt}{"end"} = $i;
				}
			}
		}
		
		#print Dumper %merged_hits; exit;
	return \%merged_hits;		# %{query}{fragmentID}{start|end}=index
	}

sub load_parsed_blast{
# loading blast file parsed by blastxml_parse.pl #
# making a hash of interval trees to determine overlap of hits and merge #
	my ($blast_file) = @_;

	open IN, $blast_file or die $!;
	
	my %itree;
	my %max_len;
	while(<IN>){
		chomp;
		if($.==1){		# header 
			die " ERROR: the blast table header appears to be formatted incorrectly!\n"
				if ! /bit/i && ! /value/i;		# checking for 'bit' and 'value' in header
			}
		else{			# body
			my @tmp = split /\t/;
			
			# loading itree #
			if(! exists $itree{$tmp[0]}){
				$itree{$tmp[0]} = Set::IntervalTree->new();
				}
			die " ERROR: start position is > end position ($tmp[7] > $tmp[8])\n"
				if $tmp[7] > $tmp[8];
				
				#print Dumper join("-", @tmp[7..8]);
			$itree{$tmp[0]} -> insert(1, $tmp[7] -1 , $tmp[8] +1);
			
			# getting max index #
			if(! exists $max_len{$tmp[0]}){
				$max_len{$tmp[0]} = $tmp[8];
				}
			else{
				$max_len{$tmp[0]} = $tmp[8] if $tmp[8] > $max_len{$tmp[0]};		
				}
			
			}
		}

		#print Dumper scalar keys %itree; exit;
	close IN;
	return \%itree, \%max_len;
	}

sub print_interval_tree{
### Description: printing out interval tree

	my ($itree, $maxrange) = @_;
	
	print join("\t", qw/position value/), "\n";
	for my $i (1..$maxrange){
		my $res = $itree->fetch($i, $i);
		$res = join(":", @$res); 
		print join("\t", $i, $res), "\n";
		}	
	}

sub load_fasta{
	# version: 2.0
	# usage: load_fasta($fasta_file_name); returns a hash of sequences
	my $fasta_in = shift;
	open IN, $fasta_in or die $!;
	my (%fasta, $tmpkey);
	while(<IN>){
		chomp;
 		$_ =~ s/#.+//;
 		next if  $_ =~ /^\s*$/;	
 		if($_ =~ />.+/){
 			$_ =~ s/>//;
 			$fasta{$_} = "";
 			$tmpkey = $_;	# changing key
 			}
 		else{$fasta{$tmpkey} .= $_; }
		}
	close IN;
		#print Dumper(%fasta); exit;
	return \%fasta;
	} 



__END__

=pod

=head1 NAME

scaffold_contamination_rm.pl -- Remove contamination from scaffolds

=head1 SYNOPSIS

scaffold_contamination_rm.pl -r -b [options] > scaffolds.fna

=head2 options

=over

=item -f 	Fasta of scaffolds.

=item -b 	Blast_2_hsp_table output.

=item -l 	Fragment Length cutoff (bp). [100]

=item -k 	Cutoff for % of scaffold that must be left after contaminantion removal [25]
  		
=item -h 	This help message

=back

=head2 For more information:

perldoc scaffold_contamination_rm.pl

=head1 DESCRIPTION

Remove contamination from scaffolds while retain uncomtaminated portions.

Contamination is identified via blastn.

The blast results must be in xml format and then parsed by blastxml_parse.pl. 

Scaffolds containing contamination are split into multiple contigs/scaffolds.
These split scaffolds can be filtered out if they are too short or if only a small
percentage of the original scaffold remains. Split scaffolds that are kept will have
'_\d+_\d+' added to each split scaffold sequence name. The first number is the 
split number and the second is the length of the fragment (includes gaps).

=head1 EXAMPLES

=head2 blastn

blastn -db nt -query scaffolds.fna -out scaffolds.blast.xml -max_target_seqs 1 -max_hsps_per_subject 1 -evalue 0.00001 -outfmt 5

=head2 Parse blast output (only hsps >= 100bp)

blastxml_parse.pl < scaffolds.blast.xml -hsp_length 100 > contam.blast.txt

=head2 Remove contamination

scaffold_contamination_rm.pl -b contam.blast.txt -f scaffolds.fna > scaffolds_clean.fna

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut
