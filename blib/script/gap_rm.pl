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

my ($error, $regex);

GetOptions(
	   "regex=s" => \$regex,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$gap_cut = 1 if ! $gap_cut;
$regex = "[Nn-]{1,}" if ! $regex;

### MAIN
foreach(@infiles){
	$_ = File::Spec->rel2abs($_);
	die " ERROR: $_ not found\n" if ! -e $_;
	write_gaps($_, $cutoff, $gap_cut);
	}


### Subroutines
sub write_gaps{
	my ($infile, $cutoff, $gap_cut) = @_;

	my $gap_cnt = 0;
	while(<>){ 
 		next if /^\s*$/ && ! eof;
 		next if /^\s*/;
 		
 		s/[Nn]/
		}
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

