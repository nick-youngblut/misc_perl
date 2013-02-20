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

my ($verbose, $gff_in, $ref_in, $reads_in);
my $extra = 100;				# maybe it should be max read length?
GetOptions(
	   "gff=s" => \$gff_in,
	   "ref=s" => \$ref_in,
	   "reads=s" => \$reads_in,
	   "extra=i" => \$extra, 			# extra beyond genes
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a gff file\n" if ! $gff_in;
die " ERROR: provide a reference assembly file (genbank or fasta)\n" if ! $ref_in;
#die " ERROR: provide a fastq file of reads (searched for gene(s))\n" if ! $reads_in;

### MAIN
# workflow:
	# load gff, & either fasta or merged genbank
	# load reads
	# pull out genes + downstream & upstream regions
	# map reads onto genes using samtools
	# make pileup output
	# call pres/abs of gene based on pileup results

# loading files #
my $gff_r = load_gff($gff_in);
my $format = check_ref_format($ref_in);
my $ref_r;
#if($format eq "fasta"){ $ref_r = load_fasta($ref_in); }
#elsif($format eq "genbank"){ $ref_r = load_genbank($ref_in); }
else{ die " ERROR: reference assembly format no recognized!\n"; }

# parsing out genes #
get_gene_region($gff_r, $ref_r, $extra);

### Subroutines

sub get_gene_region{
# Usage:
# Description:
# 	Getting nucleotide sequence of gene +/- extra sequence around gene
	my ($gff_r, $ref_r, $extra) = @_;
	
	my %gene_list;
	foreach my $gene (@$gff_r){
		$gene_list{$$gene[8]} = substr($ref_r, $$gene[3] - $extra, $$gene[4] + $extra - $$gene[3]);
		}

		print Dumper %gene_list; exit;
	}
	
sub load_fasta{
# Usage:
# Description:
# 	Loading fasta file sequence as string
	my $fasta_in = shift;
	
	open IN, $fasta_in or die $!;
	my $fasta;
	while(<IN>){
		chomp;
		next if />/;
		$fasta .= $_;
		}
	close IN;
		#print Dumper $fasta; exit;
	return $fasta;
	} #end load_fasta

sub check_ref_format{
# Usage:
# Description:
# 	checking format of refence assembly (fasta or genbank)
	my ($ref_in) = @_;

	open IN, $ref_in or die $!;
	
	my $format;
	while(<IN>){
		next if $format && $. > 1;
		if(/LOCUS/){ $format = "genbank"; }
		elsif(/>/){ $format = "fasta"; }
		}
	close IN;
	
	return $format;
	}

sub load_gff{
# Usage:
# 	load_gff($gff_in)
# Description:
# 	loads a gff file (RAST export: gff3)
	my ($gff_in) = @_;
	open IN, $gff_in or die $!;
	
	my @gff;
	while(<IN>){
		chomp;
		push(@gff, [split /\t/]);
		}
	close IN;
		#print Dumper @gff; exit;
	return \@gff;
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

