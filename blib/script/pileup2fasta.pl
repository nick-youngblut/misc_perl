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

my ($verbose, $cov_cutoff, $len_cutoff);
GetOptions(
	   "coverage=s" => \$cov_cutoff,
	   "length=i" => \$len_cutoff,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$len_cutoff = 100 if ! $len_cutoff;
$cov_cutoff = "MS1" if ! $cov_cutoff;
die " ERROR: Coverage cutoff must be a number or 'MS#'\n" if $cov_cutoff !~ /^\d+$/ && $cov_cutoff !~ /^MS\d+$/i; 

### MAIN
$cov_cutoff = get_coverage_cutoff($cov_cutoff) if $cov_cutoff =~ /^MS\d+/i; 		# coverage cutoff >= mean+stdev
print STDERR "Coverage cutoff: $cov_cutoff\n";
pileup2fasta($cov_cutoff, $len_cutoff);

### Subroutines
sub get_coverage_cutoff{
# getting mean & stdev of genome coverage #
	my $cov_cutoff = shift;
	
	my @cov;
	while(<>){
		my @tmp = split /\t/;
		push(@cov, $tmp[3]);
		}
	my $mean = average(\@cov);
	my $stdev = stdev(\@cov);	
	
	(my $multiple = $cov_cutoff) =~ s/^MS//;
	return int ($mean + $stdev * $multiple);
	}

sub pileup2fasta{
# pulling out regions of high coverage #
	my ($cov_cutoff, $len_cutoff) = @_;
	
	die " ERROR: provide an mpileup file via STDIN\n" if -t STDIN;
	seek STDIN, 0, 0;
	
	my %scaffold;
	while(<>){
		my @tmp = split /\t/;

		if($tmp[3] >= $cov_cutoff){							# if change to high coverage region
			$scaffold{$tmp[0]}{"start"} = $tmp[1];			# start of region
			push(@{$scaffold{$tmp[0]}{"cov"}}, $tmp[3]);		# coverage 
			while(<>){
				my @tmp = split /\t/;
				if($tmp[3] < $cov_cutoff && $scaffold{$tmp[0]}{"seq"}){					# if end of high coverage region
					my $seq_name = join("_", $tmp[0], join("-", $scaffold{$tmp[0]}{"start"}, $tmp[1]), int average($scaffold{$tmp[0]}{"cov"}) );
					print join("\n", ">$seq_name", $scaffold{$tmp[0]}{"seq"}), "\n"
						if ($tmp[1] - $scaffold{$tmp[0]}{"start"}) >= $len_cutoff;
					delete $scaffold{$tmp[0]};
					last;
					}
				elsif(! exists $scaffold{$tmp[0]}){						# if high coverage, but transition to next scaffold
					foreach my $name (keys %scaffold){
						my $seq_name = join("_", $name, join("-", $scaffold{$tmp[0]}{"start"}, $tmp[1]), int average($scaffold{$tmp[0]}{"cov"}) );
						print join("\n", ">$seq_name", $scaffold{$name}{"seq"}), "\n" 
							if ($tmp[1] - $scaffold{$tmp[0]}{"start"}) >= $len_cutoff;
						delete $scaffold{$name};
						}
					last;
					}
				push(@{$scaffold{$tmp[0]}{"cov"}}, $tmp[3]);		# coverage 			
				$scaffold{$tmp[0]}{"seq"} .= $tmp[2];			# adding nucleotides
				}	
			}
		}
	print "\n";
	}

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

__END__

=pod

=head1 NAME

pileup2fasta.pl -- pull out regions of high coverage from a pileup file

=head1 SYNOPSIS

pileup2fasta.pl [options] < input > output

=head2 options

=over

=item -c 	Coverage cutoff [>= mean + stdev]

=item -l 	Length cutoff [>=100]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc pileup2fasta.pl

=head1 DESCRIPTION

Pull out regions of high coverage from an mpileup output file, which
can be useful for blasting the regions.

=head2 Coverage cutoff options:

=over

=item "MS#" 	Mean coverage + standard deviation * '#' [default = 1]

=item Number	User supplied cutoff value

=back

=head2 Output sequence name format = ">NAME_x-y_z"

=over

=item NAME 	= Scaffold/chromosome name

=item x 	= Start position in scaffold

=item y 	= End position in scaffold

=item z 	= Mean coverage across fragment

=back

=head1 EXAMPLES

=head2 Finding regions with a coverage >=100 reads & a length >=200bp

pileup2fasta.pl -c 100 -l 200 < file.pileup > file.fna

=head2 Finding regions with a coverage >= mean + 2 *stdev of the whole-genome coverage

pileup2fasta.pl -c MS2 < file.pileup > file.fna

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

