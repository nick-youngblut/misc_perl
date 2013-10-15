#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::DB::EUtilities;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my $sep = "\t";
my $taxid_col = 1;
my ($div_bool, $sci_bool, $genus_bool, $sp_bool, $subsp_bool);
GetOptions(
	   "column=i" => \$taxid_col, 
	   "delimiter=s" => \$sep,
	   "division" => \$div_bool,
	   "scientific" => \$sci_bool,
	   "genus" => \$genus_bool,
	   "species" => \$sp_bool,
	   "subspecies" => \$subsp_bool,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$taxid_col--;

### MAIN
while(<>){
 	chomp;
 	if(/^\s*$/){ print; next; }
 	
 	# parsing line #
 	my @l = split /$sep/;
 	my @taxids = split / *; */, $l[$taxid_col];
	
	# getting names from taxIDs #
	my @out;
	foreach my $id (@taxids){
		next unless $id =~ /^\d+$/;
		my $factory = Bio::DB::EUtilities->new(-eutil => 'esummary',
                                       -email => 'mymail@foo.bar',
                                       -db    => 'taxonomy',
                                       -id    => $id );
    	
		# names #
		my %names;
		my $doc_cnt = 0;
		for my $doc ($factory->next_DocSum()){
			$names{$doc_cnt}{"div"} = [$doc->get_contents_by_name('Division')]
				if $div_bool;
			$names{$doc_cnt}{"sci"} = [$doc->get_contents_by_name('ScientificName')]
				if $sci_bool;
			$names{$doc_cnt}{"genus"} = [$doc->get_contents_by_name('Genus')]
				if $genus_bool;
			$names{$doc_cnt}{"sp"} = [$doc->get_contents_by_name('Species')]
				if $sp_bool;
			$names{$doc_cnt}{"subsp"} = [$doc->get_contents_by_name('Subsp')]
				if $subsp_bool;
			
			$doc_cnt++;
			}

		push @out, \%names;
        }
        
    # assembling output #
    my @ids;
    foreach my $id (@out){
    	my @ll;
    	foreach my $doc_cnt (keys %{$id}){
			push @ll, ${$id->{$doc_cnt}{"div"}}[0] 
				if $div_bool && ${$id->{$doc_cnt}{"div"}}[0];
			push @ll, ${$id->{$doc_cnt}{"sci"}}[0] 
				if $sci_bool && ${$id->{$doc_cnt}{"sci"}}[0];
			push @ll, ${$id->{$doc_cnt}{"genus"}}[0] 
				if $genus_bool && ${$id->{$doc_cnt}{"genus"}}[0];
			push @ll, ${$id->{$doc_cnt}{"sp"}}[0] 
				if $sp_bool && ${$id->{$doc_cnt}{"sp"}}[0];
			push @ll, ${$id->{$doc_cnt}{"subsp"}}[0] 
				if $subsp_bool && ${$id->{$doc_cnt}{"subsp"}}[0];						
			}
		push @ids, join("|", @ll) if @ll;
    	}
    
    # writing new line #
    $l[$taxid_col] = join(";", @ids) if @ids;
    print join($sep, @l), "\n";
 	}


__END__

=pod

=head1 NAME

taxid2name.pl -- get taxonomy names for a list of ncbi taxonomy IDs

=head1 SYNOPSIS

taxid2name.pl [options] < input.txt > output.txt

=head2 options

=over

=item -column  <int>

Column number containing taxIDs (indexed by 1). [1]

=item -delimiter  <char>

Delimiter separating columns in table. ['\t']

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

=item -v	Verbose output

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

=head2 Just scientific names

taxid2name.pl -sci < taxids.txt > taxids_sci.txt

=head2 Just division, genus, & species

taxid2name.pl -div -sp -sub < taxids.txt > taxids_div-sp-subsp.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

