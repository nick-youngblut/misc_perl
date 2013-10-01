#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::DB::GenBank;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
GetOptions(
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults


### MAIN
my $acc_r = load_accessions();

my $gb = new Bio::DB::GenBank;
my $seqio = $gb->get_Stream_by_acc( $acc_r );

# header #
print join("\t", qw/accession acc_classification acc_authors acc_description man_title man_authors man_location man_pubmed/), "\n";

# getting info from seq objects #
my @ref_q = qw/title authors location pubmed/;
while( my $seq =  $seqio->next_seq ) {
	# seqfeature #
	my @line = (
		$seq->accession_number,
		join(";", reverse($seq->species->classification()) ),
		$seq->authority,
		$seq->desc
		);

	# annotations #
	my @annotations = $seq->annotation->get_Annotations('reference');
	foreach my $val (@annotations){
		my $hash_r = $val->hash_tree;
		foreach my $q (@ref_q){
			if(exists $hash_r->{$q}){
				push @line, $hash_r->{$q};
				}
			else{
				push @line, "";
				}
			}
		}
	
	# checking for undefined variables #
	map{$_ = "" unless $_} @line;
	
	# writing data for accession #
	print join("\t", @line), "\n";
    }

### Subroutines
sub load_accessions{
	my @acc;
	while(<>){
		chomp;
		next if /^\s*$/;
		
		my @line = split /\t/;
		push @acc, $line[0];
		}
	return \@acc;
	}


__END__

=pod

=head1 NAME

accession2genbankInfo.pl -- Get Genbank info associated with accession numbers

=head1 SYNOPSIS

accession2genbankInfo.pl [options] < accessions.txt > genbank_info.txt

=head2 options

=item -help			

This help message

=back

=head2 For more information:

perldoc accession2genbankInfo.pl

=head1 DESCRIPTION

Provide a list of accession numbers
(1 number per line; can be 1st column of table)
and get metadata info associated with
the accession number.

=head2 Key for column names

=over

=item * 'sub_' = info directly related to the accession number

=item * 'man_' = info related to the manuscript associated with 
the accession number

=back

=head2 Reqirements:

Bioperl -> use Bio::DB::GenBank

=head1 EXAMPLES

=head2 Basic usage

accession2genbankInfo.pl < accessions.txt > genbank_info.txt

=head2 Using output from blastxml_parse.pl

cat blast_result.xml | blastxml_parse.pl | cut -f 2 | accession2genbankInfo.pl

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

