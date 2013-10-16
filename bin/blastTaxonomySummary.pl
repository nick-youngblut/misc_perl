#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use List::Util qw/max/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my @columns = (1);
my $NA = "";
GetOptions(
	   "columns=i{,}" => \@columns,		# columns for grouping
	   "NA=s" => \$NA, 					# to replace "N/A"
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
map{$_--} @columns;

### MAIN
load_blast_table();


### Subroutines
sub load_blast_table{
	
	my @vals = qw/qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids/;
	my %blast; 
	while(<>){
		chomp;
		next if /^\s*$/;
		
		# parsing line #
		my @l = split /\t/;
		die " ERROR: not enough columns in blast table!\n"
			unless scalar @l >= 12;
		
		# loading hash #
		#for my $i (0..$#vals){	
		for (my $i = $#vals-1; $i>=0; $i--){
			my $ii = $#vals - $i;
			$blast{ join("_::_", @l[@columns]) }{$vals[$i]} = $l[$#l-$ii]
			}
		
		}
	
	print Dumper %blast; exit;
	}


__END__

=pod

=head1 NAME

template.pl -- script template

=head1 SYNOPSIS

template.pl [options] < input > output

=head2 options

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc template.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2
   3) Step 3

=head1 EXAMPLES

=head2 Usage method 1

template.pl <read1.fastq> <read2.fastq> <output directory or basename>

=head2 Usage method 2

template.pl <library file> <output directory or basename>

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

