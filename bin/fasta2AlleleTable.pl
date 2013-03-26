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

my ($verbose);
GetOptions(
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
my $alleles_r = find_alleles();
get_multi_alleles($alleles_r);
write_positions($alleles_r);
fasta2allele($alleles_r);

### Subroutines
sub fasta2allele{
# converting fasta to allele table 
	print STDERR "...writing out allele table\n";
	my ($allele_r) = @_;
	
		#print Dumper $allele_r; exit;
	
	my (@line, $taxon);
	my $pos = 0;
	while(<>){
		chomp;
 		$_ =~ s/#.+//;
 		next if  $_ =~ /^\s*$/;	
 		if($_ =~ />.+/){			# taxon name
 			print "\n" if $pos;
 			$_ =~ s/^>//;
 			print "$_\t";
 			$pos = 0;
 			}
 		else{ 						# nucleotide seq
 			$_ = uc $_;
 			my @line = split //;
 			my @tmp; 
 			for my $i (0..$#line){
 					#print Dumper $pos, $line[$i] if exists $allele_r->{$pos}{$line[$i]};
 				push(@tmp, $allele_r->{$pos + 1}{$line[$i]})
 					if exists $allele_r->{$pos + 1}{$line[$i]};
 				$pos++;
 				}
 			print join("\t", @tmp);
 			}
		}
	print "\n";
	} #end load_fasta

sub write_positions{
# writing out SNP positions #
	my $alleles_r = shift;
	print join("\t", "Taxon", sort {$a<=>$b} keys %$alleles_r), "\n";
	}

sub get_multi_alleles{
# removing homo-allelic positions (i.e. not a snp) #
	print STDERR "...Purging all homo-alleles\n";
	my ($alleles_r) = @_;
	
	my %summary;
	foreach my $pos (keys %$alleles_r){
		if(scalar values %{$alleles_r->{$pos}} == 1){
			$summary{1}++;
			delete $alleles_r->{$pos};
			}
		else{
			my $cnt = 0;
			my $gap = 0;
			foreach my $nuc (keys %{$alleles_r->{$pos}}){
				if($nuc =~ /[ATCG]/i){			# not counting any SNPs caused by gaps or other characters
					$cnt++;
					$alleles_r->{$pos}{$nuc} = $cnt;
					}
				else{
					delete $alleles_r->{$pos};
					$gap++;
					last;
					}
				}
			$gap ? $summary{"gap or odd nucleotide"}++  : 
				$summary{scalar keys %{$alleles_r->{$pos}}}++;
			}
		}

	# printing Summary #
	print STDERR "### allele summary ###\n";
	foreach my $key (sort keys %summary){
		print STDERR join("\t", $key, $summary{$key}), "\n";
		}
	print STDERR "\n";
	}

sub find_alleles{
# finding alleles from fasta alignment #
	print STDERR "...Making allele table\n";
	
	my %alleles;
	my @nucs = ("A", "T", "G", "C");
	my $taxon;
	while(<>){
		chomp;
 		$_ =~ s/#.+//;
 		next if  $_ =~ /^\s*$/;	
 		if($_ =~ />.+/){
 			$taxon = $_;
 			}
 		else{ 
 			foreach my $nuc (@nucs){
 				print STDERR "$taxon -> $nuc\n" unless $verbose;
 				while ($_ =~ /$nuc/ig){
 					$alleles{$-[0] + 1}{$nuc} = 1;
 					}
 				}
 			}
		}
	
	seek(STDIN, 0, 0);
	
	return \%alleles;
	} #end load_fasta


__END__

=pod

=head1 NAME

fasta2AlleleTable.pl -- Convert fasta alignment to SNP-allele table

=head1 SYNOPSIS

fasta2AlleleTable.pl [options] < align.fna > align.snp

=head2 options

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc fasta2AlleleTable.pl

=head1 DESCRIPTION

Make a table of SNP (multi-allele) positions. Each allele (different
nucleotide) is represented as a number in the table.

A heatmap of the table can be plotted in R.

Positions with gaps and other characters other than [ACGTacgt]
will be skipped.

=head1 EXAMPLES

=head2 Usage:

fasta2AlleleTable.pl < align.fna > snp.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

