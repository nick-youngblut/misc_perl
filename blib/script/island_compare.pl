#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path qw/rmtree/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $maf_in, $ntaxa, $gap_rm);
my $max_len = 80000;
my $min_len = 8000;
my $LCB_sort = "descending";
my $fasta_dir = "./GI_fasta/";
my $repeat_min = 30;
my $end_region = 300;
GetOptions(
	   "maf=s" => \$maf_in,
	   "max=i" => \$max_len,
	   "min=i" => \$min_len,
	   "n=i" => \$ntaxa,
	   "sort=s" => \$LCB_sort,
	   "directory=s" => \$fasta_dir,
	   "gap" => \$gap_rm,
	   "repeat=i" => \$repeat_min, 			# '-n' in nucmber
	   "end=i" => \$end_region,				# region of LCBs that must have direct repeat
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a maf file!\n" unless $maf_in;
die " ERROR: $maf_in not found!\n" unless -e $maf_in;
die " ERROR: provide the number of taxa in the genome alignment!\n" unless $ntaxa;
$fasta_dir = File::Spec->rel2abs($fasta_dir);

### MAIN
my ($islands_r, $LCBs_r) = parse_n_filter_maf($maf_in, $min_len, $max_len, $ntaxa);

# PA table #
my $LCB_sort_r = sort_LCB_by_prevalence($islands_r, $LCB_sort);
make_PA_table($islands_r, $LCB_sort_r);

# table of start - end for genomic islands #

exit;

# writing GI sequences #
make_fasta_dir($fasta_dir);
my $file_list_r = write_LCB_seqs($islands_r, $fasta_dir);

# getting direct repeats #
get_direct_repeats($file_list_r, $islands_r, $repeat_min);


### Subroutines
sub get_direct_repeats{
	my ($file_list_r, $islands_r, $repeat_min) = @_;

	foreach my $infile (@$file_list_r){
		open PIPE, "repeat-match -n $repeat_min $infile |";

		# getting LCB - taxon #
		my @parts = File::Spec->splitpath($infile);
		$parts[2] =~ s/\.[^.]+$//;
		my ($LCB, $taxon) = split /__/, $parts[2], 2; 
		$LCB =~ s/^LCB//;
		die " LOGIC ERROR: LCB$LCB -> $taxon not found in GI hash!\n"
			unless exists $islands_r->{$LCB}{$taxon};
		
		# region end #
		my $region_end;
		if($islands_r->{$LCB}{$taxon}{"start"} < 
			$islands_r->{$LCB}{$taxon}{"end"} ){
			$region_end = $islands_r->{$LCB}{$taxon}{"end"} - $end_region;
			}
		else{
			$region_end = $islands_r->{$LCB}{$taxon}{"start"} - $end_region;
			}
						
		while(<PIPE>){
			chomp;
			next if /^\s*$/;
			my @line = split / +/;		# s1=1, s2=2, len=3
			next if $line[0];
			next if $line[1] =~ /Start/;

			#print Dumper @line; exit;
			
			if($line[1] <= $end_region	&&	#checking 1st of repeat pair; must be within $end_region of start#
				$line[2] >= $region_end){	#checking 1st of repeat pair; must be >= length of LCB - $end_region
				$islands_r->{$LCB}{$taxon}{"DR"} = 1;
				}
			}
		close PIPE;
		}
	
		print Dumper %$islands_r; exit;
	}

sub write_LCB_seqs{
	my ($islands_r, $fasta_dir) = @_;

	my @file_list;
	foreach my $LCB (keys %$islands_r){
		foreach my $taxon ( keys %{$islands_r->{$LCB}} ){
			# making output file #
			my $outfile = "$fasta_dir/LCB$LCB\__$taxon.fna";
			open OUT, ">$outfile" or die $!;
			push(@file_list, $outfile);
			

			# removing gaps #
			$islands_r->{$LCB}{$taxon}{"sequence"} =~ s/-//g unless $gap_rm;
				
			# writing sequence #
			my $name = join("__", $taxon, $islands_r->{$LCB}{$taxon}{"scaffold_id"});
			print OUT join("\n", ">$name", $islands_r->{$LCB}{$taxon}{"sequence"}), "\n";
			}
		close OUT;
		}
	return \@file_list;
	}

sub make_fasta_dir{
	my $fasta_dir = shift;
	rmtree($fasta_dir) if -d $fasta_dir;
	mkdir $fasta_dir;
	}

sub sort_LCB_by_prevalence{
	my ($islands_r, $direction) = @_;
	
	# number of taxa per LCB #
	my %LCB_prev;
	foreach my $LCB (keys %$islands_r){
		$LCB_prev{$LCB} = scalar keys %{$islands_r->{$LCB}};
		}
	
	# sort & return #
	if(! $direction || $direction =~ /ascending/i){
		return [sort{ $LCB_prev{$a} <=> $LCB_prev{$b}} keys %LCB_prev];
		}
	elsif($direction =~ /descending/i){
		return [sort{ $LCB_prev{$b} <=> $LCB_prev{$a}} keys %LCB_prev];
		}
	}

sub make_PA_table{
	my ($islands_r, $LCB_sort_r) = @_;
	
	# getting taxa #
	my %taxa;
	foreach my $LCB (keys %$islands_r){
		map{$taxa{$_} = 1} keys %{$islands_r->{$LCB}};
		}
	my @taxa = keys %taxa;
	
	# LABELS #
	my @LCBs;
	if($LCB_sort_r){
		@LCBs = @$LCB_sort_r;
		}
	else{ 
		@LCBs = keys %$islands_r;
		}
	print join("\t", "LABELS", @LCBs), "\n";
	
	# COLORS #
	print join("\t", "COLORS", ("#000000") x scalar @LCBs), "\n";
	
	# PA #
	my @PA;
	foreach my $taxon (@taxa){
		my @row;
		foreach my $LCB (@LCBs){
			if(exists $islands_r->{$LCB}{$taxon}){
				# adding length #
				push(@row, 
					abs( $islands_r->{$LCB}{$taxon}{"end"} -
						 $islands_r->{$LCB}{$taxon}{"start"})
						 );
				}
			else{
				push(@row, 0);
				}
			}
		print join("\t", $taxon, @row), "\n";
		}
	
	}

sub parse_n_filter_maf{
# parsing & filtering maf #
	my ($maf_in, $min_len, $max_len, $ntaxa) = @_;
	
	open IN, $maf_in or die $!;
	my %islands;
	my $LCB_skip = 0;
	my $label;
	while(<IN>){
		chomp;
		next if /^#/; 		# skipping comment lines

		if(/^a/){
			my @line = split /[= ]/;			# mult = [7]
			$label = $line[4]; 
			$LCB_skip = 1 if $line[6] >= $ntaxa; 		# skip LCB if core
			}
		elsif(! $LCB_skip && /^s/){
			# parse line #
			my @line = split /[\t ]+/;
			my @name = split /[. ]/, $line[1];
			
			# filter by lcb length #
			my $LCB_len = $line[3];

			if($LCB_len > $max_len || $LCB_len < $min_len){
				$LCB_skip = 1;
				next;
				}
			else{
				# loading hash #
				$islands{$label}{$name[0]}{"scaffold_id"} = $name[1];		# lcb=>taxa=>cat=value
				$islands{$label}{$name[0]}{"sequence"} = $line[$#line];
				$islands{$label}{$name[0]}{"start"} = $line[2];
				$islands{$label}{$name[0]}{"end"} = $line[2] + $line[3];
				$islands{$label}{$name[0]}{"strand"} = $line[4];
				$islands{$label}{$name[0]}{"scaffold_len"} = $line[5];
				$islands{$label}{$name[0]}{"DR"} = 0;
				}
			
			}
		elsif(/^\s*$/ ){
			$LCB_skip = 0; 		# reset
			}
		}
	close IN;
	
	
		#print Dumper %islands; exit;
	return \%islands;
	}


__END__

=pod

=head1 NAME

template.pl -- script template

=head1 SYNOPSIS

template.pl [options] < input > output

=head2 options

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc template.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2
   3) Step 3

=head1 EXAMPLES

=head2 Usage method 1

template.pl <read1.fastq> <read2.fastq> <output directory or basename>

=head2 Usage method 2

template.pl <library file> <output directory or basename>

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

