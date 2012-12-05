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
$cov_cutoff = "MS" if ! $cov_cutoff;
die " ERROR: Coverage cutoff must be a number, 'M', or 'MS'\n" if $cov_cutoff !~ /^\d+$/ && $cov_cutoff !~ /^MS*$/i; 

### MAIN
$cov_cutoff = get_coverage_cutoff($cov_cutoff) if $cov_cutoff =~ /^MS*$/i; 		# coverage cutoff >= mean+stdev
	#print Dumper $cov_cutoff; exit;
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
	
	if ($cov_cutoff =~ /^M$/i){ 
		return $mean; 							# mean
		}			
	elsif ($cov_cutoff =~ /^MS$/i){
		my $stdev = stdev(\@cov);	
		return $mean + $stdev;					# mean + stdev
		}
	else{ die " LOGIC ERROR! $!\n"; }
	}

sub pileup2fasta{
# pulling out regions of high coverage #
	my ($cov_cutoff, $len_cutoff) = @_;
	
	seek STDIN, 0, 0;
	
	my %scaffold;
	while(<>){
		my @tmp = split /\t/;

		if($tmp[3] >= $cov_cutoff){						# if change to high coverage region
			$scaffold{$tmp[0]}{"start"} = $tmp[1];
			while(<>){
				my @tmp = split /\t/;
					#print Dumper $tmp[3], $cov_cutoff;
				if($tmp[3] < $cov_cutoff){					# if end of high coverage region
					print join("\n", $tmp[0], $scaffold{$tmp[0]}{"seq"}), "\n";
					delete $scaffold{$tmp[0]};
					last;
					}
				elsif(! exists $scaffold{$tmp[0]}){						# if high coverage, but transition to next scaffold
					foreach my $name (keys %scaffold){
						print join("\n", $name, $scaffold{$name}{"seq"}), "\n";
						delete $scaffold{$name};
						}
					last;
					}
								
				$scaffold{$tmp[0]}{"seq"} .= $tmp[2];			# adding nucleotides
				}	
			}
		}
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

pileup2fasta.pl -- script pileup2fasta

=head1 SYNOPSIS

pileup2fasta.pl [options] < input > output

=head2 options

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc pileup2fasta.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2
   3) Step 3

=head1 EXAMPLES

=head2 Usage method 1

pileup2fasta.pl <read1.fastq> <read2.fastq> <output directory or basename>

=head2 Usage method 2

pileup2fasta.pl <library file> <output directory or basename>

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

