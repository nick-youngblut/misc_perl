#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::DB::EUtilities;
use Bio::DB::Taxonomy;
use Forks::Super;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my $sep = "\t";
my $taxid_col = 1;
my ($div_bool, $sci_bool, $genus_bool, $sp_bool, $subsp_bool, $lin_bool);
my $fork = 0;
GetOptions(
	   "column=i" => \$taxid_col, 
	   "delimiter=s" => \$sep,
	   "division" => \$div_bool,
	   "scientific" => \$sci_bool,
	   "genus" => \$genus_bool,
	   "species" => \$sp_bool,
	   "subspecies" => \$subsp_bool,
	   "lineage" => \$lin_bool,
	   "fork=i" => \$fork,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide file name!\n" unless $ARGV[0];
die " ERROR: cannot find $ARGV[0]!\n" unless -e $ARGV[0];
$taxid_col--;

### MAIN
my $taxids_r = get_taxids($ARGV[0]);
my $tax_info_r = get_tax_info($taxids_r);
taxid2taxinfo($ARGV[0], $tax_info_r);


### Subroutines
sub taxid2taxinfo{
	my ($infile, $tax_info_r) = @_;
	
	# status #
	print STDERR " Converting taxids to seleted taxonomy info\n" unless $verbose;
	
	open IN, $infile or die $!;
	while(<IN>){
		chomp;
		if(/^\s*$/){ print; next; }
		 
		# parsing line #
	 	my @l = split /$sep/;
 		my @taxids = split / *; */, $l[$taxid_col];
 		
 		# assembling output #
    	my @ids;
  		foreach my $id (@taxids){
    		next unless exists $tax_info_r->{$id};
    		my @ll;
    		foreach my $doc_cnt (keys %{$tax_info_r->{$id}}){
				push @ll, ${$tax_info_r->{$id}{$doc_cnt}{"div"}}[0]
					if $div_bool && ${$tax_info_r->{$id}{$doc_cnt}{"div"}}[0];
				push @ll, ${$tax_info_r->{$id}{$doc_cnt}{"lineage"}}[0]
					if $lin_bool && ${$tax_info_r->{$id}{$doc_cnt}{"lineage"}}[0];
				}
			push @ids, join("|", @ll) if @ll;
			}
			
		# writing new line #
	    $l[$taxid_col] = join(";", @ids) if @ids;
    	print join($sep, @l), "\n";
    	}
    
	close IN;
	}

sub get_taxids{
	my ($infile) = @_;

 	# status #
 	print STDERR " Getting all taxids\n" unless $verbose;

	open IN, $infile or die $!;

	my %taxids;
	while(<IN>){
	 	chomp;
 
	 	# parsing line #
 		my @l = split /$sep/;
	 	my @tmp = split / *; */, $l[$taxid_col];	
	 	map{$taxids{$_} = 1 unless $_ !~ /^\d+$/} @tmp;
	 	}
 	seek(IN, 0, 0);
 	close IN;
 	
 	# status #
 	print STDERR " Number of unique taxids found: ", scalar keys %taxids, "\n"
 		unless $verbose;
 		
 	return \%taxids;
 	}

sub get_tax_info{
	my ($taxids_r) = @_;
	
	# status #
	print STDERR " Using Bio::DB::EUtilities to get seleted info for each unique taxid\n"
		unless $verbose;
	
	my %tax_info;
	foreach my $taxid (keys %$taxids_r){
		my $job = fork {
			max_proc => $fork,
			share => [ \%tax_info ],
			sub => sub{
				my $factory = Bio::DB::EUtilities->new(-eutil => 'esummary',
                                   -email => 'mymail@foo.bar',
                                      -db    => 'taxonomy',
                                      -id    => $taxid );
				
				# names #
				my $doc_cnt = 0;
				for my $doc ($factory->next_DocSum()){
					$tax_info{$taxid}{$doc_cnt}{"div"} = [$doc->get_contents_by_name('Division')]
						if $div_bool;
					$tax_info{$taxid}{$doc_cnt}{"sci"} = [$doc->get_contents_by_name('ScientificName')]
						if $sci_bool;
					$tax_info{$taxid}{$doc_cnt}{"genus"} = [$doc->get_contents_by_name('Genus')]
						if $genus_bool;
					$tax_info{$taxid}{$doc_cnt}{"sp"} = [$doc->get_contents_by_name('Species')]
						if $sp_bool;
					$tax_info{$taxid}{$doc_cnt}{"subsp"} = [$doc->get_contents_by_name('Subsp')]
						if $subsp_bool;
					$doc_cnt++;
					}
			
				$tax_info{$taxid}{0}{"lineage"} = get_lineage($taxid) if $lin_bool;
				}
			};
		}
	waitall;
		
		#print Dumper %tax_info; exit;
	return \%tax_info;
	}

sub get_lineage{
# if lineage needed #
	my $id = shift;

	my @levels = qw/superkingdom phylum class order family genus species/;

	my %lineage;	
	my $db = new Bio::DB::Taxonomy(-source => 'entrez');
	my $node= $db->get_Taxonomy_Node(-taxonid => $id);

	# get basal node #
	my $anc = $node;
	while(1){
		last unless $anc->ancestor();
		$anc = $anc->ancestor();
		my $rank = $anc->rank;
		if( grep(/$rank/i, @levels) ){
			$lineage{$rank} = $anc->scientific_name;
			}
		}
		
	my @desc = $db->get_all_Descendents($node);

	for my $child ( @desc ) {
		my $rank = $anc->rank;
		if( grep(/$rank/i, @levels) ){
			$lineage{$rank} = $anc->scientific_name;
			}
		}
	
	# checking for each level #
	map{ $lineage{$_} = "NA" unless $lineage{$_} } @levels;

		#print join(":", @lineage{@levels}), "\n"; exit;
	return [join(":", @lineage{@levels})];
	}

__END__

=pod

=head1 NAME

taxid2name.pl -- get taxonomy names for a list of ncbi taxonomy IDs

=head1 SYNOPSIS

taxid2name.pl [options] input.txt > output.txt

=head2 options

=over

=item -column  <int>

Column number containing taxIDs (indexed by 1). [1]

=item -delimiter  <char>

Delimiter separating columns in table. ['\t']

=item -lineage  <bool>

Include entire taxonomic lineage? [FALSE]

=item -division  <bool>

Include division name? [FALSE]

=item -scientific  <bool>

Include scientific name? [FALSE]

=item -genus  <bool>

Include genus name? [FALSE]

=item -species  <bool>

Include species name? [FALSE]

=item -subspecies  <bool>

Include subspecies name? [FALSE]

=item -fork  <int>

Number of parallel EUtilities queries. [1]

=item -verbose  <bool>

Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc taxid2name.pl

=head1 DESCRIPTION

Bio::DB::EUtilities is used to get taxonomy info
for ncbi taxIDs. The input should be a delimited table
(tab-delimited by default).

Multiple taxonomy IDs in the same column should be
separated with a ';'

The output is the same table with the taxonomy IDs 
replaced by taxonomy names. Names are separated by "|" and
taxonomy ID entries are separated by ";"

A row in the table will not be altered unless a valid
taxonomy ID can be found.

=head1 EXAMPLES

=head2 Just taxonomic lineage

taxid2name.pl -lin taxids.txt > taxids_lineage.txt

=head2 Just scientific names

taxid2name.pl -sci taxids.txt > taxids_sci.txt

=head2 Just division, genus, & species

taxid2name.pl -div -sp -sub taxids.txt > taxids_div-sp-subsp.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

