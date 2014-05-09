#!/usr/bin/perl
my $mod = "10/10/12 2:06 PM";
my $version = "0.6";
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

### global variables
my ($error);

### I/O
my ($verbose, $method, $cutoff, $fasta_w);
GetOptions(
	   "method=i" => \$method,
	   "verbose" => \$verbose,
	   "fasta" => \$fasta_w,
	   "cutoff=i" => \$cutoff,
	   "help|?" => \&usage # Help
	   );

### Input error check
$cutoff = 0 if ! $cutoff;
$method = 1 if ! $method;

### Routing main subroutines
write_lengths($fasta_w, $cutoff, $method);

#----------------------Subroutines----------------------#
sub write_lengths{
	my ($fasta_w, $cutoff) = @_;
	
	my (%fasta, $tmpkey);
	while(<>){ 
 		chomp;
 		next if /^\s*$/ && ! eof;
 		if(eof){
 			die " ERROR: file seems to have no sequences!\n" if ! $tmpkey; 			
 			$fasta{$tmpkey} .= $_;
 			count_or_print(\%fasta, $tmpkey, $fasta_w, $cutoff, $method);
 			last;
 			}
 		if($_ =~ />.+/){		# if new header 
 			if($tmpkey){		# if sequence already loaded
 				count_or_print(\%fasta, $tmpkey, $fasta_w, $cutoff, $method);
 				}
 			$fasta{$_} = "";
 			$tmpkey = $_;	# changing key
 			}
 		else{$fasta{$tmpkey} .= $_; }
		}
	
	sub count_or_print{
		### counting or printing fasta ###
		my ($fasta_r, $tmpkey, $fasta_w, $cutoff, $method) = @_;
	 		#print Dumper $fasta_r; exit;
 		my $len = 0;
  		if ($method == 1){ $len++ while $$fasta_r{$tmpkey} =~ /[ATGCRYSWKMBDHVN]/gi; }
  		elsif ($method == 2){ $len++ while $$fasta_r{$tmpkey} =~ /[ATGC]/gi; }
  		else{ die " ERROR: method $method not recognized!\n"; }
  		
  		
  		print $len, "\n" if $len >= $cutoff and ! $fasta_w;
  		print "$tmpkey\n$$fasta_r{$tmpkey}\n" if $fasta_w and $len >= $cutoff;
  		
  		delete $$fasta_r{$tmpkey};
		}
	
	} #end load_fasta

sub error_routine{
	my $error = $_[0];
	my $exitcode = $_[1];
	print STDERR "ERROR: $error\nSee help: [-h]\n";
	exit($exitcode);
	}

sub usage {
 my $usage = <<HERE;
Usage:
 contig_lengths.pl [-f] [-c] [-m] < contig_file > output.txt
Options:
 -f 	Write fasta instead of lengths
 -c 	Length cutoff for writing output (>= cutoff)
    	 [default: 0]
 -m 	Method of counting:
    	 1 = all IUPAC letters 'ATGCRYSWKMBDHVN'
    	 2 = just 'ATGC'	(doesn't cound N-gaps)
    	 [1]
Description:
 The program produces a list of contig lengths
 for each contig in the fasta file.
 * only fasta formatted input allowed!
 * full IUPAC nucleotide code used.
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