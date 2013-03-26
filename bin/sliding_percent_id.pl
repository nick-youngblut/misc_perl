#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::AlignIO;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($format_in, $verbose, $total_bool, $pop_in);
my $window = 100;
my $jump = 100;
GetOptions(
	   "in=s" => \$format_in,
	   "window=i" => \$window, 		# window size
	   "jump=i" => \$jump,			# jump size
	   "total" => \$total_bool,		# total percent ID? [FALSE]
	   "population=s" => \$pop_in, 	# population table file
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
if(! $format_in){
	print STDERR "SeqIO is guessing input file format.\n (See http://www.bioperl.org/wiki/HOWTO:AlignIO_and_SimpleAlign for formats)\n";
	}


### MAIN
# loading alignment #
my ($alnin, $alnout);
print STDERR "...loading alignment\n";
eval{
	$alnin = Bio::AlignIO->new(-fh => \*STDIN,
							-format => $format_in);
	};
if($@){		# if error
  	print STDERR "Was not able to open files, sorry!\n";
   	print STDERR "Full error is\n\n$@\n";
   	exit(-1);
	}

### calculating percent ID ###
# if population table #
if($pop_in){				
	# loading pop table #
	my ($pops_r, $pop_list_r) = load_pop($pop_in);
	my ($name_ord_r);
	while (my $aln = $alnin->next_aln){			# alignment
		# name order for parsing alignment #
		$name_ord_r = get_name_order($pops_r, $aln);
		add_name_order($pops_r, $name_ord_r);

		# pairwise comparisons of populations #	
		for my $i (0..$#$pop_list_r){
			my $pop_aln1 = $aln->select_noncont(@{$pops_r->{$$pop_list_r[$i]}});
		
			for my $ii (0..$#$pop_list_r){
				next if $i <= $ii; 		# lower triangle
				my $comp = join("__", $$pop_list_r[$i], $$pop_list_r[$ii]);
				print STDERR "...calculating percent ID: Pop $$pop_list_r[$i] vs Pop $$pop_list_r[$ii]\n";
	
				my $pop_aln2 = $aln->select_noncont(@{$pops_r->{$$pop_list_r[$ii]}});
				my $pop_aln12 = $aln->select_noncont(@{$pops_r->{$$pop_list_r[$i]}},
									@{$pops_r->{$$pop_list_r[$ii]}});
				
				# calculating #
				window_percent_id($pop_aln1, $window, $jump, $comp, "within_pop1");
				window_percent_id($pop_aln2, $window, $jump, $comp, "within_pop2");
				window_percent_id($pop_aln12, $window, $jump, $comp, "between");
				}
			}
		}
	}
	
# just all alignment #
else{	
	while (my $aln = $alnin->next_aln){
		print STDERR "...calculating percent ID\n";
	
		write_total_percent_id($aln) if $total_bool;

		window_percent_id($aln, $window, $jump);
		}
	}


### Subroutines
sub add_name_order{
# adding name order to populations #
	my ($pops_r, $name_ord_r) = @_;
	
	foreach my $pop (keys %$pops_r){
		foreach my $i (0..$#{$pops_r->{$pop}}){
			if(exists $name_ord_r->{${$pops_r->{$pop}}[$i]} ){
				${$pops_r->{$pop}}[$i] = $name_ord_r->{${$pops_r->{$pop}}[$i]}; 
				}
			else{
				die " ERROR: ${$pops_r->{$pop}}[$i] not found in alignment!\n";
				}
			}
		}
		#print Dumper %$pops_r; exit;
	}

sub get_name_order{
# getting position of each taxon in the alignment #
	my ($pops_r, $aln) = @_;
	
	my %name_ord;
	my $cnt = 0;
	#while (my $name = $aln->each_seq->display_id()){
	foreach my $seq ($aln->each_seq){
		$cnt++;
		$name_ord{$seq->display_id()} = $cnt;
		}
		#print Dumper %name_ord; exit;
	return \%name_ord; 
	}

sub load_pop{
# loading population table #
	my ($pop_in) = @_;
	open IN, $pop_in or die $!;
	my %pops;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @line = split /\t/;
		die " ERROR: population table is not formatted correctly\n"
			unless scalar @line == 2;
		push(@{$pops{$line[1]}}, $line[0]);
			#$pops{$line[1]}{$line[0]} = 1;
		}
	close IN;
		#print Dumper %pops; exit;
	return \%pops, [keys %pops];
	}

sub write_total_percent_id{
# writing percent ID for total alignment #
	my ($aln) = @_;
	print join("\t", 0, $aln->length, $aln->percentage_identity), "\n";
	}

sub window_percent_id{
# getting percent ID for a window #
	my ($aln, $window, $jump, $comp, $group) = @_;
	
	for (my $i=1; $i<=($aln->length -1); $i+=$jump){
		print STDERR "$i," if $verbose;
		
		my $sub_aln = $aln->slice($i, $i + $window - 1);		# subalignment
		my $pID = $sub_aln->percentage_identity;				# percent ID
		$pID = "NA" if ! $pID;
		print join("\t", $i, $i+$window -1, $pID, $comp, $group), "\n";		# writing
		}
	}


__END__

=pod

=head1 NAME

sliding_percent_id.pl -- Sliding window of percent sequence identity

=head1 SYNOPSIS

sliding_percent_id.pl [options] < alignment > percentID.txt

=head2 options

=over

=item -i 	Alignment file format. [bioperl will guess]

=item -w 	Window length (bp). [100]

=item -j 	Jump length (bp). [100]

=item -p	Population table (2 column: taxon	population).

=item -t 	Get percent ID of total alignment length? [FALSE]

=item -v 	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc sliding_percent_id.pl

=head1 DESCRIPTION

Get percent sequence identity for a sequence alignment over a sliding scale.

Just using Bio::AlignIO subroutines.

If a population table is provided (*txt; 2 column), then each population is
compared via within vs between for each window. Two extra columns will be
added to the output: the comparison, & which populations are being compared 
(within each population or between).

=head1 EXAMPLES

=head2 Usage: basic

sliding_percent_id.pl -i fasta < file.fna > file_percent-ID.txt

=head2 Usage: comparing populations

sliding_percent_id.pl -i fasta -p pops.txt < file.fna > file_percent-ID-pop.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

