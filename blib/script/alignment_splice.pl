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

my ($verbose, $table_in, $range);
my $cutoff = 0;
GetOptions(
	   "cutoff=i" => \$cutoff,
	   "table=s" => \$table_in,			# table of values 
	   "range=s" => \$range,			# range of values to keep
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
my %splice_tbl;
if($table_in){
	load_splice_table(\%splice_tbl, $table_in, $cutoff);
	}
else{
	load_range(\%splice_tbl, $range);
	}

my $fasta_r = load_fasta();
splice_fasta($fasta_r, \%splice_tbl);

### Subroutines
sub splice_fasta{
# splicing fasta based on splice table #
	my ($fasta_r, $splice_tbl_r) = @_;
	
	foreach my $taxon (sort keys %$fasta_r){
		my $seq;
		foreach my $start (sort {$a<=>$b} keys %$splice_tbl_r){
			$splice_tbl_r->{$start} = length($fasta_r->{$taxon}) 
				if $splice_tbl_r->{$start} eq "end";
 
			$seq .= substr($fasta_r->{$taxon}, 
				$start -1, $splice_tbl_r->{$start} - $start + 1);
			}
		print join("\n", "$taxon", $seq), "\n";
		}
	}

sub load_fasta{
# version: 2.0
# usage: load_fasta($fasta_file_name); returns a hash of sequences
	my $fasta_in = shift;
	my (%fasta, $tmpkey);
	while(<>){
		chomp;
 		$_ =~ s/#.+//;
 		next if  $_ =~ /^\s*$/;	
 		if($_ =~ />.+/){
 			$fasta{$_} = "";
 			$tmpkey = $_;	# changing key
 			}
 		else{$fasta{$tmpkey} .= $_; }
		}
		#print Dumper(%fasta); exit;
	return \%fasta;
	} #end load_fasta

sub load_splice_table{
# loading table for splicing #
# 3 column: start, end, value (optional) #
	my ($splice_tbl_r, $table_in, $cutoff) = @_;

	open IN, $table_in or die $!;
	while(<IN>){
		chomp;
		my @line = split /[\t ]+/;
		if($line[2] && $cutoff){
			next unless $line[2] >= $cutoff;
			}
		
		$splice_tbl_r->{$line[0]} = $line[1];
		}
	close IN;
	
	}
	
sub load_range{
# loading range of numbers to splice #
	my ($splice_tbl_r, $range) = @_;
	my @ranges = split /,/, $range;
	
	foreach my $r (@ranges){
		my @line = split /-/, $r;
		$line[0] = 0 if ! $line[0];
		$line[1] = "end" if ! $line[1];
		$splice_tbl_r->{$line[0]} = $line[1];
		}
		
		#print Dumper $splice_tbl_r; exit;
	}



__END__

=pod

=head1 NAME

alignment_splice.pl -- splice regions of an alignment (fasta format)

=head1 SYNOPSIS

alignment_splice.pl [-r | -t] [-c] < align.fna > align_splice.fna

=head2 options

=over

=item -range 	

List start-end values for regions to keep. If the end value is not given,
the max alignment length is used. Example: "1-10,15-20,25-"

=item -table

2 or 3 column table (*txt) with start, end, and percent sequence ID (optional).
These regions will be kept. This table can be produced by sliding_percent_id.pl

=item -cutoff

Cutoff value for percent ID column (regions with values < cutoff will not
be kept. 

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc alignment_splice.pl

=head1 DESCRIPTION

Remove regions of a fasta alignment. The regions kept can be given via a list
(-r) or with a table (-t). Percent sequence ID and a cutoff (-c) can be used
to determine which regions are kept.

=head2 WARNING:

Overlapping ranges will produce repeating regions in the alignment!

=head1 EXAMPLES

=head2 Usage: provided range

alignment_splice.pl -r 1-10,15- < file.fna > splice.fna

=head2 Usage: table from sliding_percent_id.pl (97% ID cutoff)

alignment_splice.pl -t file_pID.txt -c 97 < file.fna > splice.fna

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

