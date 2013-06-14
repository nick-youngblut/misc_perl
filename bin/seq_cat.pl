#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
GetOptions(
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
my %cat;
my $taxon_cnt = 0;
foreach my $infile (@ARGV){
	die " ERROR: $infile not found!" unless -e $infile;
	
	# loading fasta #
	my $fasta_r = load_fasta($infile);
	check_lengths($fasta_r, $infile); 
	
	# checking number of taxa #
	my $file_taxa_cnt = scalar keys %$fasta_r;
	if($taxon_cnt && $taxon_cnt > $file_taxa_cnt){
		die " ERROR: $infile has too few taxa!\n";
		}
	elsif($taxon_cnt && $taxon_cnt < $file_taxa_cnt){
		die " ERROR: $infile has too many taxa!\n";
		}		
	
	# concatenating #
	foreach my $taxon (keys %$fasta_r){
		$cat{$taxon} .= $fasta_r->{$taxon};
		}
	}
	
# checking lengths #
check_lengths(\%cat, "the concatenated alignment");
write_fasta(\%cat);

### Subroutines
sub write_fasta{
	my $fasta_r = shift;
	foreach my $taxon (keys %$fasta_r){
		print "$taxon\n$fasta_r->{$taxon}\n";
		}
	}

sub check_lengths{
	my ($fasta_r, $infile) = @_;
	my %lens;
	my %taxa;
	foreach my $taxon (keys %$fasta_r){
		my $len = length $fasta_r->{$taxon};
		$lens{$len}++;
		$taxa{$taxon}{$len} = 1;
		}
	if( (scalar keys %lens) > 1 ){	# if multiple lengths	
		print STDERR " ERROR: $infile has sequences of varying lengths!\n";
		foreach my $taxon (keys %taxa){
			foreach my $len (sort{$b<=>$a} keys %{$taxa{$taxon}}){
				print STDERR "$taxon -> $len\n";
				}
			}
		exit(1);
		}

	}

sub load_fasta{
	my $fasta_in = shift;
	open IN, $fasta_in or die $!;
	my (%fasta, $tmpkey);
	while(<IN>){
		chomp;
 		 s/#.+//;
 		next if  /^\s*$/;	
 		if(/>.+/){
 			$fasta{$_} = "";
 			$tmpkey = $_;	# changing key
 			}
 		else{$fasta{$tmpkey} .= $_; }
		}
	close IN;
		#print Dumper %fasta; exit;
	return \%fasta;
	} #end load_fasta



__END__

=pod

=head1 NAME

seq_cat.pl -- concatenate fasta files of alignments

=head1 SYNOPSIS

seq_cat.pl *fasta > concat_aligned.fasta

=head2 options

=over

=item -h	This help message

=back

=head2 For more information:

perldoc seq_cat.pl

=head1 DESCRIPTION

Concatenate fasta files of alignments (e.g. aligned core gene clusters).
Sequence names must match!
The script will die otherwise!
Both constistent lengths of individual alignments & the final alignment
are checked.

=head1 EXAMPLES

=head2 Usage: 3 aligned core gene clusters

seq_cat.pl align1.fna align2.fna align3.fna > concat_aln.fna

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

