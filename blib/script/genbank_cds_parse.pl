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

my ($verbose, $nuc_bool, $regex, $invar);
GetOptions(
	   "nucleotide" => \$nuc_bool, 			# nucleotide or protein sequences? [protein]
	   "regex=s" => \$regex,				# regex for selecting annotations 
	   "invariant" => \$invar,				# captialization invariant for regex? [TRUE]
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
if($regex && $invar){ $regex = qr/$regex/; }
elsif($regex){ $regex = qr/$regex/i; }


### MAIN
foreach my $infile (@ARGV){
	my $mfile = merge_genbank($infile);
	my $seqs_r = get_cds_features($mfile, $nuc_bool, $regex);
	
	if(scalar @ARGV > 1){ write_seqs_fasta($seqs_r, $infile); }
	else{ write_seqs_fasta($seqs_r); }
	}


### Subroutines
sub write_seqs_fasta{
	my ($seqs_r, $infile) = @_;
	
	$infile =~ s/\.[^\.]+$// if $infile;
	
	foreach my $gid (sort {$a<=>$b} keys %$seqs_r){
		foreach my $product (keys %{$$seqs_r{$gid}}){
			#print join("\n", join("", ">", join("__", $infile, "CDS$gid", "$product")), $$seqs_r{$gid}{$product}), "\n";
			if($infile){
				print join("\n", ">$infile\__CDS$gid\__'$product'", $$seqs_r{$gid}{$product}), "\n";
				}
			else{
				print join("\n", ">CDS$gid\__'$product'", $$seqs_r{$gid}{$product}), "\n";
				}
			}
		}
	}

sub get_cds_features{
	my ($infile, $nuc_bool, $regex) = @_;

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
		
		# nuc or AA #
		if($nuc_bool){		# nucleotide
			$gene_sequences{$i + 1}{ $product[0] } = $cds_features[$i]->spliced_seq->seq;
			}
		else{
			$gene_sequences{$i + 1}{ $product[0] } = $trans[0];
			}
			#print Dumper %gene_sequences; exit;
		}

		#print Dumper %gene_sequences; exit;
	return \%gene_sequences;
	}

sub merge_genbank{
# merging a multi-contig genbank file using 'union' from EMBOSS #
	my $infile = shift;
	
	(my $outfile = $infile) =~ s/\.[^\.]+$|$/_merged.gbk/;
	my $cmd = "union -sequence $infile -sformat genbank -outseq $outfile -osformat genbank -auto -feature";
	
	print STDERR "### Merging genbank file ###\n";
	print STDERR $cmd, "\n";
	system($cmd);
	
	return $outfile;
	}

sub example{
foreach(@ARGV){
	die " ERROR: $_ not found\n" if ! -e $_;
	my $seqioo = Bio::SeqIO->new(-file => $_);
	for my $seq_o ($seqioo->next_seq){
		
		my $anno_collection = $seq_o->annotation;
		print Dumper $anno_collection; exit;
		
		for my $feat_o ($seq_o->get_SeqFeatures){
			next if $feat_o->primary_tag ne "CDS";

     		if ($feat_o->has_tag('translation') && $feat_o->has_tag('product')) {
        		for my $val ($feat_o->get_tag_values('translation')){
            		print "trana: ",$val,"\n";
    	        	}
	    	     }		
			}
		}
	}
	}


__END__

=pod

=head1 NAME

genbank_cds_parse -- parse out CDS records from a genbank file

=head1 SYNOPSIS

genbank_cds_parse [options] file(s).gbk > out.fasta

=head2 options

=over

=item -n

AA or Nucleotide output? [AA]

=item -r 

Regular expression for parsing by 'product.'
Example: '^Methyl.+alpha'

=item -i

Capitalization invariant for regex? [TRUE]

=item -v	

Verbose output

=item -h	

This help message

=back

=head2 Requires:

=over 

=item *

Bioperl

=item *

EMBOSS

=back

=head2 For more information:

perldoc genbank_cds_parse

=head1 DESCRIPTION

Parse out all or a seletion of CDS features from a genbank file.
The CDS features are converted to fasta and can be written in AA or
nucleotide format.

EMBOSS is used to merge multi-record genbank files.

=head2 Output fasta name format:

=over

=item 1 input file:

">cds#__'product'"


=item >1 input file:

">file__cds#__'product'"

=back

=head1 EXAMPLES

=head2 Parse all CDS features

genbank_cds_parse.pl file1.gbk > feats.fa

=head2 Regex search (2 files)

genbank_cds_parse.pl file1.gbk file2.gbk -r 'coenzyme M reductase alpha' > mcrA.fa

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

