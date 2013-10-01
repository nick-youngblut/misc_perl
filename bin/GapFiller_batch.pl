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

my ($verbose, $contig_in, $read1_in, $read2_in);
my $overlap = 90;
my $coverage = 5;
my $insert_span = 0.2;
GetOptions(
	   "contig=s" => \$contig_in, 			# scaffold file
	   "r1=s" => \$read1_in,				# f-read
	   "r2=s" => \$read2_in,				# r-read
	   "overlap=i" => \$overlap,			# -m in GapFiller
	   "coverage=i" => \$coverage, 			# -o in GapFiller
	   "insert_range=f" => \$insert_span,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a contig file and read files\n" if ! -e $contig_in || ! -e $read1_in || ! -e $read2_in;
$contig_in = File::Spec->rel2abs($contig_in);
$read1_in = File::Spec->rel2abs($read1_in);
$read2_in = File::Spec->rel2abs($read2_in);

### MAIN

# chdir #
my @parts = File::Spec->splitpath($contig_in);
chdir $parts[1] or die $!;

# finding insert size #
my $insert_size = call_estimate_insert_size($contig_in, $read1_in, $read2_in);
make_lib_file($read1_in, $read2_in, $insert_size, $insert_span);
call_GapFiller($parts[2], $contig_in, $overlap, $coverage);

### subroutines ###
sub call_GapFiller{
	my ($base, $contig_in, $overlap, $coverage) = @_;

	(my $basename = $base) =~ s/\.[^.]+$|$/_fill/;
	
	my $cmd = "GapFiller.pl -l lib.txt -s $contig_in -m $overlap -o $coverage -b $basename";
	print STDERR "$cmd\n" unless $verbose;
	system($cmd); 
	}

sub make_lib_file{
	my ($read1_in, $read2_in, $insert_size, $insert_span) = @_;
	open OUT, ">lib.txt" or die $!;
	print OUT join(" ", "lib1", "bowtie", $read1_in, $read2_in, $insert_size, $insert_span, "FR"), "\n";
	close OUT;
	}

sub call_estimate_insert_size{
	my ($contig_in, $read1_in, $read2_in) = @_;
	my $cmd = "estimate_insert_size.pl $contig_in $read1_in $read2_in";
	print STDERR "$cmd\n" unless $verbose;
	my $out = `$cmd`;
	$out =~ s/.+median = (\d+).+/$1/s;
	die " ERROR: \"$out\" is not numberic. Cannot determine median insert size!\n"
		unless $out =~ /^\d+$/;
	
	return $out;
	}



__END__

=pod

=head1 NAME

GapFiller_batch.pl -- batch running of GapFiller.pl

=head1 SYNOPSIS

GapFiller_batch.pl [flags]

=head2 required flags

=item -contig

Contig/scaffold file.

=item -r1 	

Forward read file (Sanger fastq).

=item -r2 	

Reverse read file (Sanger fastq).

=back

=head2 optional flags

=over

=item -overlap

Minimum number of overlapping bases with the edge of the gap ('-m' in GapFiller). [90]

=item -coverage

Minimum number of reads needed to call a base during an extension ('-o' in GapFiller). [5]

=item -insert_range

Range around the provided insert size (provided in library file). [0.2]

=item -verbose

Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc GapFiller_batch.pl

=head1 DESCRIPTION

A wrapper around GapFiller. The median insert size estimated from
estimate_insert_size.pl is used for GapFiller run.


=head1 EXAMPLES

=head2 Usage

GapFiller_batch.pl -contig contigs.fna -r1 reads_F.fq -r2 reads_R.fq

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

