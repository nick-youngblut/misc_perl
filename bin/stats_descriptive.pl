#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Statistics::Descriptive;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my $sep = '\t';
my $round = 2;
GetOptions(
		"round=i" => \$round,
		"sep=s" => \$sep,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$round = join("", "%.", $round, "f");

### MAIN
my %cols;
while(<>){
	chomp;
	next if /^\s*$/;
	
	# loading each column #
	my @line = split /$sep/;
	for my $i (0..$#line){
		if(! exists $cols{$i} ){
			$cols{$i} = Statistics::Descriptive::Full->new();
			}
		$cols{$i}->add_data($line[$i]);
		}
	}

# writing out stats #
foreach my $col (sort{$a<=>$b} keys %cols){
	my $col_i = $col + 1;
	print join("\t", $col_i, "min", sprintf($round, $cols{$col}->min()) ), "\n";
	print join("\t", $col_i, "Q1", sprintf($round, $cols{$col}->percentile(25)) ), "\n";
	print join("\t", $col_i, "mean", sprintf($round, $cols{$col}->mean()) ), "\n";
	print join("\t", $col_i, "median", sprintf($round, $cols{$col}->median()) ), "\n";
	print join("\t", $col_i, "Q3", sprintf($round, $cols{$col}->percentile(75)) ), "\n";
	print join("\t", $col_i, "max", sprintf($round, $cols{$col}->max()) ), "\n";
	print join("\t", $col_i, "stdev", sprintf($round, $cols{$col}->standard_deviation()) ), "\n";
	}



__END__

=pod

=head1 NAME

stats_descriptive.pl -- summary stats for a table of values

=head1 SYNOPSIS

stats_descriptive.pl [options] < input > stats.txt

=head2 options

=over

=item -sep

Column seperator. ['\t']

=item -round

Number of decimals to round summary values. [2]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc stats_descriptive.pl

=head1 DESCRIPTION

Get summary descriptive statistics of a list or table of numbers.

=head2 Output format

=over

=item column 1

Column number (starting at 1)

=item column 2

Statistic (e.g. 'mean')

=item column 3

Value

=back

=head1 EXAMPLES

=head2 Usage: list of values

seq 1 1 1000 | stats_descriptive.pl

=head2 Usage: 2 column table

paste <(seq 1 1 1000) <(seq 1 1 1000) | stats_descriptive.pl

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

