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
my $keep_cut = 50;
my $len_cut = 500;
GetOptions(
	   "fasta=s" => \$ref_file,
	   "blast=s" => \$blast_file,
	   "coverage=i" => \$cov_cut,
	   "sequence=i" => \$seqID_cut,
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

#foreach my $q (keys %$itree_r){ print_interval_tree($$itree_r{$q}, $$max_len_r{$q}); } exit;

my $merged_hits_r = merge_hits($itree_r, $max_len_r);
my $fasta_r = load_fasta($ref_file);										# loading scaffolds
$fasta_r = scaffold_split($fasta_r, $merged_hits_r, $keep_cut, $len_cut);
write_fasta($fasta_r);


#my $tmpfile = blast_2_hsp_table_2_fasta($blast_file);						# loading blast_2_hsp_table results
#my $res_ref = nucmer_map($tmpfile, $ref_file, $cov_cut, $seqID_cut);		# mapping sequences
#scaffold_split($fasta_ref, $res_ref, $keep_cut, $len_cut);

#----------------------Subroutines----------------------#
sub write_fasta{
# writing out fasta 
	}

sub scaffold_split{
# removing sections of scaffold based on blast hits # 
	my ($fasta_r, $merged_hits_r, $keep_cut, $len_cut) = @_;
	
	my %filter_stats;
	my %filter_log;
	$filter_stats{"seq_starting"} = scalar keys %$fasta_r;
	foreach my $q (keys %$merged_hits_r){
		die " ERROR: $q not found in fasta file!\n" if ! exists $$fasta_r{$q};
		
		# converting from string to array and pulling out fragments #		
		my $seq_len = length $$fasta_r{$q};
		$$fasta_r{$q} = [ (split //, $$fasta_r{$q} )[@{$$merged_hits_r{$q}}] ];
		$$fasta_r{$q} = join("", @{$$fasta_r{$q}});
		
		# filtering #
		my $cnt;
		$cnt++ while ($$fasta_r{$q} =~ /[ATCG]/g);
		
		my $seq_len_rm = length $$fasta_r{$q};
		$filter_log{$q}{"length_starting"} =  $seq_len;
		$filter_stats{"length_starting_total"} += $seq_len;
		$filter_log{$q}{"length_cut"} =  $seq_len - $seq_len_rm;
		$filter_log{$q}{"length_remaining"} =  $seq_len_rm;
		$filter_stats{"length_cut_total"} +=  $seq_len - $seq_len_rm;

		if( (length $$fasta_r{$q}) / $seq_len * 100 < $keep_cut ||
			$cnt < $len_cut ){
			delete $$fasta_r{$q};
			$filter_stats{"seq_rm"}++;
			$filter_log{$q}{"kept_remove"} = "removed";
			}
		else{
			$filter_log{$q}{"kept_remove"} = "kept";
			}
		
		}
	#$filter_stats{"seq_remaining"} = $filter_stats{"seq_rm"} / $filter_stats{"seq_starting"} * 100;

	write_filter_log(\%filter_stats, \%filter_log);
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
	print STDERR "Sequences starting: 	", $$log_r{"seq_starting"}, "\n";
	print STDERR "Sequences removed: 	", $$log_r{"seq_rm"}, "\n";
	print STDERR "Percent remaining:	", 
		sprintf("%.1f", ($$log_r{"seq_starting"} - $$log_r{"seq_rm"}) / $$log_r{"seq_starting"} * 100), 
		"%\n\n";
	
	print STDERR "### Length cut stats ###\n";
	print STDERR "Length starting:	", $$log_r{"length_starting_total"}, "\n";
	print STDERR "Length cut: 		", $$log_r{"length_cut_total"}, "\n";
	print STDERR "Percent remaining:	", 
		sprintf("%.1f", ($$log_r{"length_starting_total"} - $$log_r{"length_cut_total"}) / $$log_r{"length_starting_total"} * 100), 
		"%\n\n";

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
	my ($itree_r, $max_len_r) = @_;
	
	my %merged_hits;
	foreach my $q (keys %$itree_r){
		# finding all hit intervals # 
		my $cnt = 0;
		for my $i (0..$$max_len_r{$q}){
			my $res = $$itree_r{$q}->fetch($i, $i);
			push(@{$merged_hits{$q}},  $i) if $res;
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
			$itree{$tmp[0]} -> insert(1, $tmp[6] -1 , $tmp[7] +1);
			
			# getting max index #
			if(! exists $max_len{$tmp[0]}){
				$max_len{$tmp[0]} = $tmp[7];
				}
			else{
				$max_len{$tmp[0]} = $tmp[7] if $tmp[7] > $max_len{$tmp[0]};		
				}
			
			}
		}
	
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


################### LEGACY ####################
sub blast_2_hsp_table_2_fasta{
	my $blast_file = shift;
	open IN, $blast_file or die $!;
	
	my $tmpfile = File::Temp->new();
	open OUT, ">$tmpfile" or die $!;

	my %fasta; 		# needed to remove duplicates
	while(<IN>){
		my @line = split /\t/;
		$line[17] =~ s/-/N/g;
		$fasta{$line[17]} = 1;
		}
	close IN;
		
	my $cnt = 1;
	foreach(keys %fasta){
		print OUT join("\n", ">$cnt", $_), "\n";
		$cnt++;
		}
	close OUT;
		
	die " ERROR: no sequences found in $blast_file\n" if ! keys %fasta;
	return $tmpfile;
	}

sub nucmer_map{

	my ($tmpfile, $ref_file, $cov_cut, $seqID_cut) = @_;
	
	print STDERR "...Nucmer mapping\n";
	
	# getting file name for output #
	my @parts = File::Spec->splitpath($ref_file);
	
	# mapping reads with nucmer #
	if($verbose){ 
		print STDERR "nucmer -maxmatch -c 100 -p .nucmer $ref_file $tmpfile\n";
		`nucmer -maxmatch -c 100 -p .nucmer $ref_file $tmpfile`; 
		}
	else{ `nucmer -maxmatch -c 100 -p .nucmer $ref_file $tmpfile 2>/dev/null`; }
	`delta-filter -r .nucmer.delta > .nucmer.delta.filter`;
	open PIPE, "show-coords -r -c -l -T -H .nucmer.delta.filter |" or die $!;
	
	# making file for hits rejected & those that were the top hits #
	open REJECT, "> nucmer.reject.txt" or die $!;
	open TOPHIT, "> nucmer.tophits.txt" or die $!;
	
	my %res;					# checking for unique (or best) hit. Best = coverage & seqID
	my @queries;				# checking for any queries that had no hits above cutoff
	while (<PIPE>){
		chomp;
		# 12 = query name
		# 10 = query coverage
		# 6 = % ID
			#print $_;				# printing out nucmer output
		my @line = split /\t/;
		push @queries, $line[12];
		
		if ($line[10] < $cov_cut || $line[6] < $seqID_cut){
			print REJECT $_;
			next;
			}
		next if exists $res{$line[12]} && $res{$line[12]}{"coverage"} > $line[10];		# next if current hit is lower coverage
		next if exists $res{$line[12]} && $res{$line[12]}{"seqID"} > $line[6]; 			# next if current hit is lower SeqID

		$res{$line[12]}{"coverage"} = $line[10];
		$res{$line[12]}{"seqID"} = $line[6];
		$res{$line[12]}{"res"} = \@line;
		}
	
	close PIPE;
	close REJECT;
	
	# removing tmp files #
	my @rm_list = qw/.nucmer.delta .nucmer.delta.filter/;
	foreach(@rm_list){ unlink $_ or die $! if -e $_; }
	
	# writing top hits #
	foreach (sort keys %res){
		print TOPHIT join("\t", @{$res{$_}{"res"}}), "\n";
		}
	close TOPHIT;	
	
	# checking to see if any queries did not make cutoffs #
	foreach(@queries){
		print STDERR " WARNING: Contig $_ did not hit any scaffold above cutoffs\n" 
			if ! exists $res{$_};
		}
	
	# if no hits #
	die " NO HITS FOUND ABOVE CUTOFF!\n" if ! keys %res;
	
	return \%res;
	}


sub usage {
 my $usage = <<HERE;
Usage:
  scaffold_contamination_rm.pl -r -b [options] > scaffolds.fna
Options:
  -r 	fasta of scaffolds.
  -b 	blast_2_hsp_table output
  -c 	coverage cutoff for nucmer hits [95]
  -s 	sequence identity cutoff for nucmer hits [95]
  -k 	cutoff for percent of scaffold that must be left 
  		after contaminantion removal [75]
  -l 	length cutoff for keeping a fragment of a scaffold [100] bp
Description:
  Remove contamination from scaffolds.
  Provide a blast table (blast_2_hsp_table) output
  that identifies contamination.
  These regions are identified in the scaffolds by nucmer mapping.
  The regions are removed and the remaining fragments of the
  scaffold are kept if the fragment length or total length of 
  scaffold remaining meets the cutoffs (-k & -l).
  A warning is provided if a contaminating contig is not found 
  in the scaffold file.
Notes:
  Version: $version
  Last Modified: $mod
  Author: $author
Categories:
  Genome assembly

HERE
	print $usage;
    exit(1);
}

__END__

=pod

=head1 NAME

scaffold_contamination_rm.pl -- Remove contamination from scaffolds

=head1 SYNOPSIS

scaffold_contamination_rm.pl -r -b [options] > scaffolds.fna

=head2 options

=over

=item -r 	fasta of scaffolds

=item -b 	blast_2_hsp_table output

=item -c 	length cutoff [95]

=item -s 	sequence identity cutoff [95]

=item -k 	cutoff for percent of scaffold that must be left after contaminantion removal [75]
  		
=item -h 	This help message

=back

=head2 For more information:

perldoc scaffold_contamination_rm.pl

=head1 DESCRIPTION

Remove contamination from scaffolds while retain uncomtaminated portions.

Contamination is identified via BLAST.

Provide a parsed blast table (output parsed by blastxml_parse.pl; normal output).

These regions are identified in the scaffolds by nucmer mapping.
The regions are removed and the remaining fragments of the scaffold are kept if 
the fragment length or total length of scaffold remaining meets the cutoffs (-k & -l).
A warning is provided if a contaminating contig is not found in the scaffold file.

=head1 EXAMPLES

=head2 Parse blast output

scaffold_contamination_rm.pl < file.blast.xml > file.blast.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut
