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

my ($verbose, $nlen);
GetOptions(
 	   "nlength=i" => \$nlen,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$nlen = 100 if ! $nlen;



### MAIN
foreach my $infile(@ARGV){
	my $infile = File::Spec->rel2abs($infile);
	die " ERROR: $infile not found!\n" if ! -e $infile;
	cat_scaffold($infile, $nlen);
	}

### Subroutines
sub cat_scaffold{
# concatenate scaffolds #
	my ($infile, $nlen) = @_;
	
	# I/O #
	open(IN, $infile) or die $!;
	(my $outfile = $infile) =~ s/\.[^\.]+$|$/_cat.fna/;
	open(OUT, ">$outfile") or die $!;

	# loading and concatenating #
	(my $cat_name = $infile) =~ s/.+\///;
	
	print OUT ">$cat_name\_cat-gap$nlen\n";
	my $first = 0;
	while(<IN>){
		chomp;
 		s/#.+//;
 		next if /^\s*$/;
 		if (/^\s*>/ && $first){ print OUT "N" x $nlen; next;}
 		elsif (/^\s*>/){ $first++; }
 		elsif(eof){ print OUT "$_\n"; }
 		else{ print OUT;}
		}
	close IN; close OUT;
	}


__END__

=pod

=head1 NAME

scaffold_cat.pl -- Concatenating scaffolds with a certain number of N's

=head1 SYNOPSIS

scaffold_cat.pl [options] file1 file2 ...

=head2 options

=over

=item -n 	Gap length [100]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc scaffold_cat.pl

=head1 DESCRIPTION

Concatenate scaffolds with a certain number of N's placed between each scaffold.
Input is one or more multi-fasta files. The file name is used to name the concatenated scaffolds.

=head1 EXAMPLES

=head2 Adding gaps of 1000 N's between each scaffold

scaffold_cat.pl -n 1000 input.fna 

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

