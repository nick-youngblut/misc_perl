#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use List::Util qw/max/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my $columns = join(" ", qw/qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids/);
my $groups = "qseqid";
my $NA = "";
GetOptions(
	   "columns=s" => \$columns,		# column names
	   "groups=s" => \$groups, 		# grouping columns
	   "NA=s" => \$NA, 					# to replace "N/A"
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
my @columns = split / /, $columns;
my @groups = split / /, $groups;

### MAIN
my $col_index_r = make_column_index(\@columns, \@groups);
my $blast_r = load_blast_table($col_index_r);
my ($tax_sum_r, $tax_total_r) = summarize_taxonomy($blast_r);
normalize_tax_sum($tax_sum_r, $tax_total_r);

# writing summary #
write_summary($tax_sum_r);

### Subroutines
sub write_summary{
	my ($tax_sum_r) = @_;
	
	my @levels = qw/superkingdom phylum class order family genus species/;
	
	foreach my $group (keys %$tax_sum_r){
		foreach my $level (@levels){
			foreach my $cls (sort keys %{$tax_sum_r->{$group}{$level}}){
				print join("\t", split(/_::_/, $group), 
							$level, $cls, $tax_sum_r->{$group}{$level}{$cls}), "\n";
				}
			}
		}
	
	}
	
sub normalize_tax_sum{
	my ($tax_sum_r, $tax_total_r) = @_;
	
	foreach my $group (keys %$tax_sum_r){
		foreach my $level (keys %{$tax_sum_r->{$group}}){
			foreach my $cls (keys %{$tax_sum_r->{$group}{$level}}){
				$tax_sum_r->{$group}{$level}{$cls} = 
					$tax_sum_r->{$group}{$level}{$cls} /
					$tax_total_r->{$group}{$level};
				}
			}
		}
	
		#print Dumper $tax_sum_r; exit;
	}

sub summarize_taxonomy{
	my ($blast_r) = @_;
	
	my @levels = qw/superkingdom phylum class order family genus species/;

	my %tax_sum;
	my %tax_total;
	foreach my $group (keys %$blast_r){							# each group #
		# counting by group #
		foreach my $l (keys %{$blast_r->{$group}}){				# each hit #
			foreach my $ll (@{$blast_r->{$group}{$l}{"staxids"}}){
				# getting the levels shared by the taxids #
				my %shared; my %total;
				my @taxids = split /;/, $ll;		# each taxonomy
				for my $taxid (@taxids){
					if($taxid eq "N/A"){
						for my $i (0..$#levels){					# each taxonomic level
							$shared{$levels[$i]}{"N/A"}++;
							$total{$levels[$i]}++;
							#$tax_total{$group}{$levels[$i]}++;							
							}
						}
					else{
						my @tmp = split /:/, $taxid;
						for my $i (0..$#tmp){		# each taxonomic level
							$shared{$levels[$i]}{$tmp[$i]}++;
							$total{$levels[$i]}++;
							#$tax_total{$group}{$levels[$i]}++;
							}
						}
					}	
				
				# summing for the blast hit #
				for my $level (keys %shared){
					for my $cls (keys %{$shared{$level}}){
						# shared; tax_sum++
						if($shared{$level}{$cls} == $total{$level}){
							$tax_sum{$group}{$level}{$cls}++;
							$tax_total{$group}{$level}++;
							}
						}
					# if none in level shared, level = "N/A" #
					unless(keys %{$tax_sum{$group}{$level}}){
						$tax_sum{$group}{$level}{"NOT_SHARED"}++;
						$tax_total{$group}{$level}++;	
						}
					}
				}
			}
		}
	
		#print Dumper %tax_sum; exit;
	return \%tax_sum, \%tax_total;
	}

sub make_column_index{
	my ($columns_r, $groups_r) = @_;

	# lower case #
	map{tr/A-Z/a-z/} @$columns_r;
	map{tr/A-Z/a-z/} @$groups_r;	
	
	
	# blast columns # 
	my @req = qw/qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore/;
	#my %req;
	#for my $i (0..$#req){ $req{$req[$i]} = $i; }
	
	# setting index #
	my %col_index;
	for my $i (0..$#$columns_r){
		if(grep(/$$columns_r[$i]/, @req)){				# blast column
			$col_index{"blast"}{$$columns_r[$i]} = $i;

			push @{$col_index{"blast_array"}}, $i;
			}
		elsif(grep(/$$columns_r[$i]/, "staxids")){		# staxid
			$col_index{"staxids"}{$$columns_r[$i]} = $i;

			push @{$col_index{"staxids_array"}}, $i;			
			}
		else{
			$col_index{"extra"}{$$columns_r[$i]} = $i;
			
			push @{$col_index{"extra_array"}}, $i;
			}
		
		if(grep(/$$columns_r[$i]/, @$groups_r)){		# group column
			push @{$col_index{"group"}}, $i;
			}
		}

		#print Dumper %col_index; exit;
	return \%col_index;
	}

sub load_blast_table{
	my ($col_index_r) = @_;
	
	#my @vals = qw/qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids/;	
	
	# parsing input #
	my %blast; 
	while(<>){
		chomp;
		next if /^\s*$/;
		
		# parsing line #
		my @l = split /\t/;
		die " ERROR: not enough columns in blast table!\n"
			unless scalar @l >= 12;
		
		# loading hash #
		## group key ##
		my $group = join("_::_", @l[@{$col_index_r->{"group"}}]); 
		
		## taxonomy values ##
		$blast{$group}{$.}{"staxids"} = [@l[@{$col_index_r->{"staxids_array"}}]];
		
		## blast values ##
		$blast{$group}{$.}{"blast"} = [@l[@{$col_index_r->{"blast_array"}}]];
		
		## extra values ##
		$blast{$group}{$.}{"extra"} = [@l[@{$col_index_r->{"extra_array"}}]];
		}
	
		#print Dumper %blast; exit;
	return \%blast;
	}


__END__

=pod

=head1 NAME

blastTaxonomySummary.pl -- summarize blast hits to various lineages

=head1 SYNOPSIS

blastTaxonomySummary.pl [options] < blast.txt > tax_summary.txt

=head2 options

=over

=item -columns

Column names in blast table. ["qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid"]

=item -groups

Column names used for grouping taxonomies. ["qseqid"]


=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc blastTaxonomySummary.pl

=head1 DESCRIPTION

Make a summary table of which lineages were hit by blast query(s).

Output is normalized by total number of hits for the group.

If multiple taxids for a hit, the taxonomic level is only kept
if it is found in all taxids. Otherwise: 'NOT_SHARED';

"N/A" = no taxid

=head1 EXAMPLES

$ cat TMA-WC_blastp_tax.txt | blastTaxonomySummary.pl -c "clade cluster qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" -g "clade"

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

