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
foreach my $infile (@ARGV){
	die " ERROR: $infile not found!\n" if ! -e $infile;
	sort_pairs($infile);
	}

## subroutines
sub sort_pairs{
# parsing out existing paired-end reads #
	my ($infile) = @_;

	open IN, $infile or die $!;

	my %pairs;
	my $pair_cnt = 0;			# counting number of pairs
	my $read_cnt = 0;			# total number of reads
	while(<IN>){	
		$read_cnt++;
		
		# load lines #
		my @line = split /\//;		# header, pair (& partition)
		die " ERROR: read names must be in old illumina format (>NAME/1 or >NAME/2)\n"
			if scalar @line != 2;
		my $nline = <IN>;

		#print Dumper @line; exit;

		# checking for pairs; counting #
		if( exists $pairs{$line[0]} ){			# if pair found
			delete $pairs{$line[0]};
			$pair_cnt += 2;
			}
		else{
			$pairs{$line[0]} = 1;
			}

		}
	close IN;

	# printing #
	print join(" ", $infile, $read_cnt, $pair_cnt), "\n";

	}

__END__

=pod

=head1 NAME

nseq-paired.pl -- Count the number of paired-end reads in a file

=head1 SYNOPSIS

nseq-paired.pl [options] file1.fasta file2.fasta ...

=head2 options

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc nseq-paired.pl

=head1 DESCRIPTION

Counts the number of paired-end reads in a fasta file containing reads.
The output is in 3 columns (space-delimited): 'file_name' 'total_reads' 'paired-end_reads'

=head2 WARNING

=over 

=item *

The file(s) must be in fasta format!

=item *

The read names must be in Illumina software format! Example: '@HWUSI-EAS100R:6:73:941:1973#0/1'

=back

=head1 EXAMPLES

=head2 Usage

nseq-paired.pl *fasta

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

