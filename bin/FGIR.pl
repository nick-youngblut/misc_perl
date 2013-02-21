#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::SeqIO;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $gff_in, $genbank_in, $read1_in, $read2_in, $regex);
my $extra = 90;				# maybe it should be max read length?
GetOptions(			
	   "gff=s" => \$gff_in,
	   "genbank=s" => \$genbank_in,
	   "r1=s" => \$read1_in,
	   "r2=s" => \$read2_in,
	   "regex=s" => \$regex, 			# regex for finding CDS of interest
	   "extra=i" => \$extra, 			# extra beyond genes
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
#die " ERROR: provide a gff file\n" if ! $gff_in;
die " ERROR: provide a reference assembly file (genbank or fasta)\n" if ! $genbank_in;
die " ERROR: provide read files in fastq format!\n" if ! $read1_in || ! $read2_in;
die " ERROR: read file not found!\n" if ! -e $read1_in || ! -e $read2_in;

$regex = qr/$regex/ if $regex;


### MAIN
# workflow:
	# load genbank
	# pull out genes + downstream & upstream regions
	# map reads onto genes using samtools
	# make pileup output
	# call pres/abs of gene based on pileup results

# loading files #
#my $gff_r = load_gff($gff_in);
my $gbk_fna_r = load_genbank_fna($genbank_in);
my $gene_se_r = get_cds_locations($genbank_in, $regex);

# parsing out genes #
my $region_list_r = get_gene_region($gene_se_r, $gbk_fna_r, $extra);
#write_region_list($region_list_r, $gene_se_r, $extra);

# mapping reads #
foreach my $id (keys %$region_list_r){
	foreach my $prod (keys %{$region_list_r->{$id}}){ 
		
		# status #
		print STDERR " Mapping reads onto: 'CDS$id; $prod'\n";

		write_seq($region_list_r, $id, $prod);
		call_bowtie($id, $read1_in, $read2_in);
		call_samtools($id);
		
		exit;
		}
	}



### Subroutines
sub call_samtools{
# Description:
#	Calling samtools to produce an mpileup file
	my ($id) = @_;
	
	# convert to bam #
	my $cmd = "samtools view -bS .CDS$id.sam > .CDS$id.bam";
	print_run($cmd);
	
	# sort #
	$cmd = "samtools sort .CDS$id.bam .CDS$id.sorted";
	print_run($cmd);

	# faidx #
	$cmd = "samtools faidx .CDS$id.fna";
	print_run($cmd);
	
	# mpileup #
	$cmd = "samtools mpileup -B -d 10000 -f .CDS$id.fna .CDS$id.sorted.bam > CDS$id.pile";
	print_run($cmd);
	

	}

sub call_bowtie{
# Description:
# 	Calling bowtie to map reads onto gene region
	my ($id, $read1_in, $read2_in) = @_;
	
	# calling #
	my $cmd = "bowtie2-build .CDS$id.fna .CDS$id";
	print_run($cmd);
	
	$cmd = "bowtie2 -x .CDS$id -1 $read1_in -2 $read2_in -S .CDS$id.sam";
	print_run($cmd);
	}
	
sub print_run{
	my $cmd = shift;
	print STDERR "\n\t", $cmd, "\n";
	system("$cmd");		
	}

sub write_seq{
# Description:
# 	Writing out individual sequence for read mapping #
	my ($region_list_r, $id, $prod) = @_;
	
	open OUT, ">.CDS$id.fna" or die $!;
	print OUT join("\n", ">$prod\_$id", $region_list_r->{$id}{$prod}), "\n";
	close OUT;
	}

sub write_region_list{
# Usage:
# Description:
# 	Writing out fasta of genes. Name = product_start_end
	my ($region_list_r, $gene_se_r, $extra) = @_;

	foreach my $id (keys %$region_list_r){
		foreach my $prod (keys %{$region_list_r->{$id}}){
			my $start = ${$gene_se_r->{$id}{$prod}}[0] - $extra;
			my $end = ${$gene_se_r->{$id}{$prod}}[1] + $extra;
			print join("\n", ">$prod\__$start\_$end",
			 	$region_list_r->{$id}{$prod}), "\n";
			}
		}
	}

sub get_gene_region{
# Usage:
# Description:
# 	Getting nucleotide sequence of gene +/- extra sequence around gene
	my ($gene_se_r, $gbk_fna_r, $extra) = @_;
	
	# status #
	print STDERR " Getting gene region (gene +/- $extra nucs)\n";
	
	my %region_list;
	foreach my $id (keys %$gene_se_r){
		foreach my $prod (keys %{$gene_se_r->{$id}} ){
			$region_list{$id}{$prod} = substr($gbk_fna_r, 
				${$gene_se_r->{$id}{$prod}}[0] - $extra, 
				${$gene_se_r->{$id}{$prod}}[1] + $extra - ${$gene_se_r->{$id}{$prod}}[0]);
			}
		}

		#print Dumper %region_list; exit;
	return \%region_list;
	}
	

sub get_cds_locations{
	my ($infile, $regex) = @_;

	# status #
	print STDERR " Parsing locations for CDS(s) of interest\n";

	# getting cds features #
	my @cds_features = grep { $_->primary_tag eq 'CDS' } Bio::SeqIO->new(-file => $infile)->next_seq->get_SeqFeatures;

	# making hash: cdsID => product => sequence #
	my %gene_sequences;
	for my $i (0..$#cds_features){
		my @product = $cds_features[$i]->get_tag_values('product');
		my @trans = $cds_features[$i]->get_tag_values('translation');
		
		# selection for certain annotations #
		next if $regex && $product[0] !~ $regex;
		
		# warning #
		print STDERR " WARNING: multiple product entries found for CDS$i\n" if (scalar @product) > 1;
		print STDERR " WARNING: multiple translation entries found for CDS$i\n" if (scalar @trans) > 1;
		
		# loading CDS location #
		push(@{$gene_sequences{$i + 1}{ $product[0] }}, 
			$cds_features[$i]->location->start,
			$cds_features[$i]->location->end);


			#print Dumper %gene_sequences; exit;
		}

		#print Dumper %gene_sequences; exit;
	return \%gene_sequences;			# returning %%@ (uniqueID, product, start-end)
	}




sub load_genbank_fna{
# Usage:
# Description:
#	Loading genbank file using seqIO; returning string of nuc assembly
	my ($genbank_in) = @_;

	# status #
	print STDERR " Loading in genbank as fasta\n";
	
	# loading using bioperl #
	my $seqin;
	eval{
		$seqin = Bio::SeqIO->new(-file => $genbank_in,
							-format => "genbank");
		}; 
	
	my $fasta;
	while (my $seq = $seqin->next_seq() ){
		$fasta .= $seq->seq;
		}
		
		#print Dumper $fasta; exit;
	return $fasta;
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

FGIR.pl -- Find Genes in Reads

=head1 SYNOPSIS

FGIR.pl [options] < input > output

=head2 options

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc FGIR.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2
   3) Step 3

=head1 EXAMPLES

=head2 Usage method 1

FGIR.pl <read1.fastq> <read2.fastq> <output directory or basename>

=head2 Usage method 2

FGIR.pl <library file> <output directory or basename>

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

