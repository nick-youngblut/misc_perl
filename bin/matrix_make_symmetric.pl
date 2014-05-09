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

my ($verbose, @line_skip, @mats);
GetOptions(
	   "line_skip=s{,}" => \@line_skip,
	   "matrix=s{,}" => \@mats,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide at least 1 matrix\n" if ! @mats;
my %line_skip;
map{$line_skip{$_} = 1} @line_skip if @line_skip;

### MAIN
foreach my $infile (@mats){
	die " ERROR: $infile not found!\n" if ! -e $infile;
	make_symmetric($infile, \%line_skip);
	}
	
### Subroutines
sub make_symmetric{
# making a lower triangle matric symmetric #
	my ($infile, $line_skip_ref) = @_;
	open IN, $infile or die $!;
	
	(my $outfile = $infile) =~ s/\.[^\.]+$|$/_sym.dist/;
	open OUT, ">$outfile" or die $!;
	
	# loading matrix #	
	my @mat;
	my $mat_line= 0;
	while (<IN>){
		next if exists $$line_skip_ref{$.};
		chomp;
		next if /^\s$/;
		$mat_line++;
		$mat[$mat_line] = [split /\t/, $_];
		}
	close IN;
	
	# loading upper triangle
	for my $i (0..$#mat){
		for my $ii (0..$#{$mat[$i]}){
			if ($ii == $#{$mat[$i]}){ 
				if($i == 0){
					$mat[$i][$ii + 1] = "\t";
					}
				else{
					$mat[$i][$ii + 1] = 0; 
					}
				}
			$mat[$ii][$i] = $mat[$i][$ii];
			}
		}
	#$mat[0][0] = "";	
	shift @{$mat[0]};
	
	# writing out matrix #
	foreach my $line (@mat){
		print OUT join("\t", @$line), "\n";
		}
	close OUT;
	
	}

__END__

=pod

=head1 NAME

matrix_make_symmetric.pl -- convert lower triangle matrix to symmetric matrix

=head1 SYNOPSIS

matrix_make_symmetric.pl [-l] -m

=head2 options

=over

=item -m 	Matrix file(s)

=item -l 	Line(s) to skip (Example: -l 1 2 3)

=item -h	This help message

=back

=head2 For more information:

perldoc matrix_make_symmetric.pl

=head1 DESCRIPTION

The script converts a lower triangle distance matrix (mothur phylip format) to
a symmetric matrix.

The symmetric matrix can then be loaded into R.

The output file name is modified from the input name ('_sym.dist').

=head1 EXAMPLES

=head2 Load matrix

matrix_make_symmetric.pl -m matrix1.dist

=head2 Load Mothur phylip matrix

matrix_make_symmetric.pl -m matrix1.dist -l 1

=head2 Load multiple matrices

matrix_make_symmetric.pl -m matrix1.dist matrix2.dist

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

