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

my ($verbose, @names);
my $cutoff = 0;
my $method = 1;
GetOptions(
	   "method=i" => \$method,		# count all characters or not
	   "names=s{,}" => \@names,		# file names
	   "verbose" => \$verbose,	   
	   "cutoff=i" => \$cutoff,		# length cutoff for counting
	   "help|?" => \&usage # Help
	   );

### I/O error & defaults
my @infiles;
if(! @ARGV){
	while(<>){
		chomp;
		push(@infiles, $_);
		}
	}
else{ @infiles = @ARGV; }

if(@names && $#names != $#infiles){
	die " ERROR: number of names does not match number of files!\n";
	}


### MAIN
my %fasta_cnt;
my $max_length = 0;
my %names;
for my $i (0..$#infiles){
	my $key;
	if(@names){ $key = $names[$i]; }
	else{ $key = $infiles[$i]; }

	print STDERR " Processing '$key'\n" if $verbose;
	($fasta_cnt{$key}, $max_length) = count_lengths($infiles[$i], $method, $cutoff, $max_length);
	
	$names{$infiles[$i]} = $names[$i] if @names;
	}

write_table(\%fasta_cnt, $max_length, \%names, \@infiles);


### Subroutines
sub write_table{
# writting table of contig lengths #
	my ($fasta_cnt_r, $max_length, $names_r, $infiles_r) = @_;

	my @keys = sort keys %$fasta_cnt_r;

	# writing header #
	print join("\t", @keys), "\n";	
	
	for my $i (0..$max_length){
		my @line;
		foreach my $key (@keys){
			if(${$$fasta_cnt_r{$key}}[$i]){
			 	push(@line, ${$$fasta_cnt_r{$key}}[$i]);
				}
			else{ push(@line, "NA"); }
			}
		print join("\t", @line), "\n";
		}
	}

sub count_lengths{
# counting length of each contigs #
	my ($infile, $method, $cutoff, $max_length) = @_;
	
	open IN, $infile or die $!;	

	my (%fasta, $tmpkey);
	my @counts;
	while(<IN>){ 
 		chomp;
 		next if /^\s*$/ && ! eof;
 		if(eof){
 			die " ERROR: file seems to have no sequences!\n" if ! $tmpkey; 			
 			$fasta{$tmpkey} .= $_;
			push(@counts, count_seq(\%fasta, $tmpkey, $cutoff, $method) );
 			last;
 			}
 		if($_ =~ />.+/){		# if new header 
 			if($tmpkey){		# if sequence already loaded
 				push(@counts, count_seq(\%fasta, $tmpkey, $cutoff, $method));
 				}
 			$fasta{$_} = "";
 			$tmpkey = $_;	# changing key
 			}
 		else{$fasta{$tmpkey} .= $_; }
		}
	close IN;
	
	@counts = sort{$b <=> $a} @counts;
	
	$max_length = $#counts if $#counts > $max_length;
	return \@counts, $max_length;
	
	sub count_seq{
		### counting or printing fasta ###
		my ($fasta_r, $tmpkey, $cutoff, $method) = @_;
	 		#print Dumper $fasta_r; exit;
 		my $len = 0;
  		if ($method == 1){ $len++ while $$fasta_r{$tmpkey} =~ /[ATGCRYSWKMBDHVN]/gi; }
  		elsif ($method == 2){ $len++ while $$fasta_r{$tmpkey} =~ /[ATGC]/gi; }
  		else{ die " ERROR: method $method not recognized!\n"; }
  		
  		delete $$fasta_r{$tmpkey};
  		
  		return $len if $len >= $cutoff;
		}
	
	}




__END__

=pod

=head1 NAME

contig_length_table.pl -- make a table of contig lengths from multiple contig files

=head1 SYNOPSIS

contig_length_table.pl [options] < file(s).fna > output

=head2 options

=over

=item -n 	

Names for each (must be in same order as fasta files)

=item -c

Length cutoff for counting contigs (>= cutoff)

=item -m 	

Method of counting [1]

=over 

=item 1

all IUPAC letters 'ATGCRYSWKMBDHVN'

=item 2

just 'ATGC'	(doesn't cound N-gaps)

=back

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc contig_length_table.pl

=head1 DESCRIPTION

Make a table of contig/scaffold lengths from a set of contig/scaffold files
(fasta format).

=head1 EXAMPLES

=head2 Usage: basic

contig_length_table.pl *.fna > lengths.txt

=head2 Usage: adding names

contig_length_table.pl file1.fna file2.fna -n name1 name2 > lengths.txt

=head2 Usage: using find

find . -name "*fna" | contig_length_table.pl -n name1 name2 > lengths.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

