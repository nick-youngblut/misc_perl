#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::TreeIO;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $tree_in, $maf_in, $format, $group_in);
GetOptions(
	"maf=s" => \$maf_in, 			# maf file input
	"tree=s" => \$tree_in,
	"format=s" => \$format,
	"groups=s" => \$group_in, 		# defining groups in the tree & maf
	"verbose" => \$verbose,
	"help|?" => \&pod2usage # Help
	);

### I/O error & defaults
die " ERROR: provide a tree file!\n" unless $tree_in;
die " ERROR: provide a maf file!\n" unless $maf_in;
$format = check_tree_format($format);

### MAIN
# loading group file #
my $group_r = load_group($group_in) if $group_in;

# getting branch length distance matrix #
my $treeio = tree_io($tree_in, $format); 
my ($brlen_dist_r, $leaves_r) = branch_dist($treeio);

# getting nucleotide shared distance matrix #
my ($nuc_dist_r, $len_dist_r) = percent_alignment_shared($maf_in, $leaves_r);

# writing table #
write_table($brlen_dist_r, $nuc_dist_r, $len_dist_r, $group_r);


### Subroutines
sub write_table{
# writing table #
	my ($brlen_dist_r, $nuc_dist_r, $len_dist_r, $groups_r) = @_;
	
	foreach my $t1 (keys %$brlen_dist_r){
		foreach my $t2 (keys %{$brlen_dist_r->{$t1}}){
			
			# nucleotides shared #
			my $nuc_dist;
			if(exists $nuc_dist_r->{$t1}{$t2}){ 
				$nuc_dist = $nuc_dist_r->{$t1}{$t2};
				}
			elsif(exists $nuc_dist_r->{$t2}{$t1}){ 
				$nuc_dist = $nuc_dist_r->{$t2}{$t1};
				}
			else{ $nuc_dist = "NA";}
			
			# length shared #
			my $len_dist;
			if(exists $len_dist_r->{$t1}{$t2}){ 
				$len_dist = $len_dist_r->{$t1}{$t2};
				}
			elsif(exists $len_dist_r->{$t2}{$t1}){ 
				$len_dist = $len_dist_r->{$t2}{$t1};
				}
			else{ $len_dist = 0;}
			
			# groups #
			my $group;
			$group = join("__", $groups_r->{$t1}, $groups_r->{$t2})
				if exists $groups_r->{$t1} and exists $groups_r->{$t2}; 
			
			# writing out table #
			print join("\t", $t1, $t2, $brlen_dist_r->{$t1}{$t2},
				$nuc_dist, $len_dist, $nuc_dist/$len_dist * 100, $group), "\n";
			}
		}
	}

sub load_group{
# loading group file (2 column: taxon\tgroup) #
	my ($group_in) = @_;

	open IN, $group_in or die $!;
	
	my %group;
	while(<IN>){
		chomp;
		s/#.+//;
		next if /^\s*$/;
		
		my @line = split /\t/;
		die " ERROR: group file should be 2 columns (*txt format)\n"
			unless scalar @line == 2;

		$group{$line[0]} = $line[1];		
		}
	close IN;
		
		#print Dumper %group; exit;
	return \%group;
	}

sub branch_dist{
# branch length distance matrix between leaves on tree #
	my ($treeio) = @_;
	
	my @leaves = $treeio->get_leaf_nodes;

	my %brlen_dist;
	for my $i (0..$#leaves){
		for my $ii (0..$#leaves){
			next if $ii <= $i; 		# lower triangle matrix
			$brlen_dist{$leaves[$i]->id}{$leaves[$ii]->id} = $treeio->distance(-nodes => [$leaves[$i],$leaves[$ii]] );
			}
		}
	
	# getting leaf order # (pairwise comparisons should be the same)
	my @leaf_ids;
	map { push @leaf_ids, $_->id } @leaves;
		
		#print Dumper %brlen_dist; exit;
	return \%brlen_dist, \@leaf_ids;
	}

sub percent_alignment_shared{
	my ($maf_in, $leaves_r) = @_;
	
	open IN, $maf_in or die $!;

	my %seqs;
	my $lcb_mult;
	my %nuc_shared;
	my %len_shared;
	while(<IN>){
		chomp;
		next if /^#/;				# skipping any comment lines
		if(/^a/){					# lcb header line
			my @line = split / /;
			print STDERR "Comparing sequences in LCB: ", join(" ", @line[(2,3,1)]), "\n";
			($lcb_mult = $line[3]) =~ s/mult=//;		# getting number of taxa in LCB
			}
		elsif(/^s/ && $lcb_mult > 1){					# sequence line
			# parsing line #
			my @line = split /\t+/;
			my @name = split /[. ]/, $line[0];		# [1] = name
			my @seq = split / /, $line[1];			# [4] = sequence
			
			$seqs{$name[1]} = $seq[4];
			#my @nucs = split //, $seq[4];
			}
		elsif(/^\s*$/ || eof(IN)){			# blank line between LCBs
			next if $lcb_mult == 1; 
			
			get_pairwise_id(\%seqs, \%nuc_shared, \%len_shared);
			
			map{ delete $seqs{$_} } keys %seqs;		# deleting keys for next round
			}
		}
	close IN;

		#print Dumper \%nuc_shared; exit;
	return \%nuc_shared, \%len_shared;
	}
	
sub get_pairwise_id{
# pairwise number of nucleotides shared #
	my ($seqs_r, $nuc_shared_r, $len_shared_r) = @_;
	
	my @taxa = keys %$seqs_r;
	foreach my $i (0..$#taxa){
		foreach my $ii (0..$#taxa){
			next if $ii >= $i;			# lower triangle
			die " ERROR: $taxa[$i] not found in the maf file\n"
				unless exists $seqs_r->{$taxa[$i]};
			die " ERROR: $taxa[$ii] not found in the maf file\n"
				unless exists $seqs_r->{$taxa[$ii]};
			
			# converting to uppercase #
			$seqs_r->{$taxa[$i]} = uc $seqs_r->{$taxa[$i]};
			$seqs_r->{$taxa[$ii]} = uc $seqs_r->{$taxa[$ii]};			
			
			# parsing sequences #			
			my @nuc1 = split //, $seqs_r->{$taxa[$i]};
			my @nuc2 = split //, $seqs_r->{$taxa[$ii]};
			die " ERROR sequences are not the same length\n!"
				if $#nuc1 != $#nuc2;
			
			# comparing nucleotides #
			my $shared_nuc = 0;
			my $shared_len = 0;
			for my $iii (0..$#nuc1){
				next if $nuc1[$iii] eq "-";
				next if $nuc2[$iii] eq "-";
				$shared_nuc++ if $nuc1[$iii] eq $nuc2[$iii]; 
				$shared_len++;
				}
			$nuc_shared_r->{$taxa[$i]}{$taxa[$ii]} += $shared_nuc;		# summing number of shared nucleotides among taxa
			$len_shared_r->{$taxa[$i]}{$taxa[$ii]} += $shared_len;		# summing length shared (no gaps)
			}
		}
	
	#print Dumper %nuc_shared; exit;
	}

sub tree_io{
	# loading tree object #
	my ($tree_in, $format) = @_;
	my $input = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format);
	my $treeio = $input->next_tree;		

	return $treeio;
	}

sub check_tree_format{
	my $format = shift;
	$format = "newick" if ! $format;
	$format =~ s/^new$/newick/i;
	$format =~ s/^nex$/nexus/i;
	die " Designated tree format ($format) not recognized.\n" if $format !~ /newick|nexus/;
	return $format;
	}



__END__

=pod

=head1 NAME

genome_align_conservation.pl -- comparing shared (alignment length or nucleotides)
vs branch length among aligned genomes

=head1 SYNOPSIS

genome_align_conservation.pl -m -t [-f] [-g] > nuc_v_brlen.txt

=head2 options

=over

=item -maf 

maf file (mugsy format).

=item -tree

Tree file (newick or nexus format).

=item -format

Tree file format (newick or nexus). [newick]

=item -group

Group file for defining comparisons.

=item -h	This help message

=back

=head2 For more information:

perldoc genome_align_conservation.pl

=head1 DESCRIPTION

Compare the number of nucleotides shared 
   
=head2 Group file format

=over 

=item * 2 column; tab-delimited.

=item * 1st column: Taxon names

=item * 2nd column: GroupID

=back 

=head2 Output columns:

=over 

=item 1) Taxon 1

=item 2) Taxon 2

=item 3) Branch length

=item 4) Nucleotides shared

=item 5) Length shared

=item 6) Percent ID

=item 7) Comparison 

=back 

=head1 EXAMPLES

=head2 Usage:

genome_align_conservation.pl -maf file.maf -t tree.nwk -g groups.txt > nuc_vs_brlen.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

