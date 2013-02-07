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
my $threads = 1;
my $kmer = 50;
GetOptions(
	   "contig=s" => \$contig_in, 			# scaffold file
	   "r1=s" => \$read1_in,			# f-read
	   "r2=s" => \$read2_in,			# r-read
	   "threads=i" => \$threads, 		# cpus
	   "kmer=i" => \$kmer,				# kmer for SEQuel
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

# running prep #
my $cmd = "prep.pl -c $contig_in -r1 $read1_in -r2 $read2_in -t $threads";
print STDERR "$cmd\n";
system("$cmd");

# finding insert size #
open IN, "prep/prep.log" or die $!;
my $insert;
while(<IN>){
	chomp;
	if(/insert-size:/){
		/\d+/;
		$insert = $&;
		}
	}
close IN;

if(! $insert){ die " ERRROR: insert size not found!\n"; }
else{ print STDERR " Insert size: $insert\n"; }

# running SEQuel #
$cmd = "SEQuel -A prep -i $insert -p $threads -k $kmer";
print STDERR "$cmd\n";
system("$cmd");

# end #
print STDERR "### complete ###";


__END__

=pod

=head1 NAME

SEQuel_batch.pl -- run prep.pl and SEQuel using insert size from prep.pl

=head1 SYNOPSIS

SEQuel_batch.pl -c -r1 -r2 [-t] [-k]

=head2 options

=over

=item -c 	Contig/scaffold file.

=item -r1 	Forward read file (Sanger fastq).

=item -r2 	Reverse read file (Sanger fastq).

=item -t 	Number of threads. [1]

=item -k 	Kmer for SEQuel run. [50]

=item -h	This help message

=back

=head2 For more information:

perldoc SEQuel_batch.pl

=head1 DESCRIPTION

A wrapper around prep.pl and SEQuel. The insert size estimated from
prep.pl is used for SEQuel run.

Files produced by both scripts are written to the same directory as
the contig/scaffold file.

=header 2 WARNING

Contig/scaffold names must not incluce '|'!

=head1 EXAMPLES

=head2 Usage (4 threads)

SEQuel_batch.pl -c scaffold.fna -r1 read_F.fq -r2 read_R.fq -t 4

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

