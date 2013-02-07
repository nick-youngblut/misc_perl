#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::SearchIO;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my $evalue = 10;
my $hsp_length = 0;
my $query_length = 0;
my $percentID = 0;
my $bit_score = 0;
my ($verbose, $taxa_summary, $header);
GetOptions(
	   "header" => \$header,								# [TRUE]
	   "evalue=s" => \$evalue,
	   "hsp_length=i" => \$hsp_length,
	   "query_length=i" => \$query_length,
	   "percentID=i" => \$percentID,
	   "bit_score=i" => \$bit_score,
	   "taxa_summary=i" => \$taxa_summary,				# make a summary of taxa hit [FALSE]
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
my $in = Bio::SearchIO->new(-format => "blastxml", -fh => \*STDIN);

# header #
if(! $taxa_summary && ! $header){
	print join("\t", qw/query_desc hit_accession hit_desc perc_id evalue bit_score hsp_length q_start q_end hit_start hit_end/), "\n";
	}

my %tax_sum;									# Summary table of taxa hit (query, hit-gi, hit-desc, count) 
while( my $result = $in->next_result ) {		# result object
	while( my $hit = $result->next_hit ) {		# hit object
		while( my $hsp = $hit->next_hsp ) {			# hsp object
     		 # filtering #
			if( $hsp->length('query') >= $query_length &&
    			$hsp->hsp_length >= $hsp_length &&
    			$hsp->percent_identity >= $percentID &&
    			$hsp->bits > $bit_score &&
    			$hsp->evalue < $evalue
    		 	) {
				if($taxa_summary && $taxa_summary == 1){
					$tax_sum{$hit->accession}{$hit->description}{$result->query_description}++;
					}
				elsif($taxa_summary && $taxa_summary ==2){
					# summing:
						# number hits
						# mean percent id
						# mean length
						# mean evalue
					$tax_sum{$hit->accession}{$hit->description}{"hsp_cnt"}++;
					$tax_sum{$hit->accession}{$hit->description}{"percentID"} += $hsp->percent_identity;
					$tax_sum{$hit->accession}{$hit->description}{"length"} += $hsp->hsp_length;
					$tax_sum{$hit->accession}{$hit->description}{"evalue"} += $hsp->evalue;
					$tax_sum{$hit->accession}{$hit->description}{"bits"} += $hsp->bits;
					}
				elsif($taxa_summary && $taxa_summary ==3){
					# summing:
						# number hits
						# mean percent id
						# mean length
						# mean evalue
					$tax_sum{$hit->accession}{$hit->description}{$result->query_description}{"hsp_cnt"}++;
					$tax_sum{$hit->accession}{$hit->description}{$result->query_description}{"percentID"} += $hsp->percent_identity;
					$tax_sum{$hit->accession}{$hit->description}{$result->query_description}{"length"} += $hsp->hsp_length;
					$tax_sum{$hit->accession}{$hit->description}{$result->query_description}{"evalue"} += $hsp->evalue;
					$tax_sum{$hit->accession}{$hit->description}{$result->query_description}{"bits"} += $hsp->bits;
					}
				else{
					print join("\t", 
						$result->query_description, 			# qseqid
						$hit->accession,						# accession
						$hit->description,					# qdescription
						$hsp->percent_identity,				# pidnet (percentID)
						$hsp->evalue,						# evalue
						$hsp->bits,							# bit score
						$hsp->hsp_length,					# length
						$hsp->start('query'),				# start and end locations
						$hsp->end('query'),
						$hsp->start('hit'),
						$hsp->end('hit')
						), "\n";
        			}
       			} 
			}
		}  
	}


# printing taxon summary table #
exit if ! $taxa_summary;
if($taxa_summary == 1){			# basic summary
	print join("\t", qw/hit_accession hit_desc query_name hit_count/), "\n" if ! $header;
	foreach my $name (keys %tax_sum){
		foreach my $desc (keys %{$tax_sum{$name}}){
			foreach my $q (keys %{$tax_sum{$name}{$desc}}){
				print join("\t", $name, $desc, $q, $tax_sum{$name}{$desc}{$q}), "\n";
				}
			}
		}
	}
elsif($taxa_summary == 2){	# if more involved summary 	
	print join("\t", qw/hit_accession hit_desc summary_category value/), "\n" if ! $header;
	foreach my $name (keys %tax_sum){
		foreach my $desc (keys %{$tax_sum{$name}}){	
			foreach my $cat (keys %{$tax_sum{$name}{$desc}}){ 		# summary category
				if ($cat eq "hsp_cnt"){
					print join("\t", $name, $desc, $cat, $tax_sum{$name}{$desc}{$cat}), "\n";
					}
				else{
					print join("\t", $name, $desc, $cat, 
						$tax_sum{$name}{$desc}{$cat} / $tax_sum{$name}{$desc}{"hsp_cnt"}), "\n";
					}
				}
			}
		}
	}
elsif($taxa_summary == 3){	# if more involved summary 	
	print join("\t", qw/hit_accession hit_desc query_desc summary_category value/), "\n" if ! $header;
	foreach my $name (keys %tax_sum){
		foreach my $desc (keys %{$tax_sum{$name}}){	
			foreach my $q (keys %{$tax_sum{$name}{$desc}}){
				foreach my $cat (keys %{$tax_sum{$name}{$desc}{$q}}){ 		# summary category
					if ($cat eq "hsp_cnt"){
						print join("\t", $name, $desc,$q, $cat, $tax_sum{$name}{$desc}{$q}{$cat}), "\n";
						}
					else{
						print join("\t", $name, $desc, $q, $cat, 
							$tax_sum{$name}{$desc}{$q}{$cat} / $tax_sum{$name}{$desc}{$q}{"hsp_cnt"}), "\n";
						}
					}
				}
			}
		}
	}

__END__

=pod

=head1 NAME

blastxml_parse.pl -- Parse blast output (xml format)

=head1 SYNOPSIS

blastxml_parse.pl [options] < file.blast.xml > file.blast.txt

=head2 options

=over

=item -taxa_summary 		

Print a summary of taxa hit

=item -evalue 			

e-value cutoff (examples: '10' or '1e-5'). [10]

=item -hsp_length 		

Hit length (gaps included)

=item -query_length 		

Query length (gaps included)

=item -percent_identity 	

Percent identity (values from 0-100)

=item -bit_score 		

Bit score

=item -header 			

Print header? [TRUE]

=item -help			

This help message

=back

=head2 For more information:

perldoc blastxml_parse.pl

=head1 DESCRIPTION

Parse the xml output from blast+ (and probably blast).

A summary of taxa hit can also be produced. The summary can be grouped by different categories.

=head1 EXAMPLES

=head2 Parse blast output

blastxml_parse.pl < file.blast.xml > file.blast.txt

=head2 Make a taxa-hit summary count (grouping by hitID & queryID)

blastxml_parse.pl -taxa 1 < file.blast.xml > file.taxa-hit.txt

=head2 Make a taxa-hit summary stats (grouping by hitID)

blastxml_parse.pl -taxa 2 < file.blast.xml > file.taxa-hit.txt

=head2 Make a taxa-hit summary stats (grouping by hitID & queryID)

blastxml_parse.pl -taxa 3 < file.blast.xml > file.taxa-hit.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

