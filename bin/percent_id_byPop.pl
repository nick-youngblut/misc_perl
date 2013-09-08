#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::AlignIO;
use Parallel::ForkManager;
use List::Util qw/sum/;

### args/flags
#pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $pop_in);
my $fork = 0;
my $procs = 1;
my $delim = " ";
my $mothur_cmd = "dist.seqs(fasta=?, calc=onegap, countends=T)";
GetOptions(
	   "population=s" => \$pop_in, 	# population table file
	   "fork=i" => \$fork,
	   "delimiter=s" => \$delim,
	   "processors=i" => \$procs,
	   "mothur=s" => \$mothur_cmd,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$delim = qr/$delim/;

### MAIN
# getting working directory #
my $cwd = File::Spec->rel2abs(File::Spec->curdir());

# loading files if provided as a list #
my $infiles_r;
if(@ARGV){ $infiles_r = \@ARGV; }
else{ $infiles_r = load_files(); }

# loading population file (if provided) #
my $pops_r = load_pop($pop_in) if $pop_in;		# taxon => population

# loop for forking  #
my $pm = new Parallel::ForkManager($fork);

foreach my $infile (@$infiles_r){
	$infile = File::Spec->rel2abs($infile);
	
	# forking #
	my $pid = $pm->start and next;
	
	# making temp directory #
	use File::Temp qw/tempdir/;
	my $tmpdir = File::Temp->newdir(); 		# temp directory
	my $dirname = $tmpdir->dirname;
	chdir($dirname) or die $!;
	
	# symlinking #
	my @parts = File::Spec->splitpath($infile);
	my $link_name = "$dirname/$parts[2]";
	symlink($infile, $link_name) or die $!;
	
	# calling mother #
	my $dist_file = call_mother($link_name, $mothur_cmd, $procs);
	
	# loading distance file & summing distances per population comparison #
	sum_distances($dist_file, $pops_r, $infile);
	
	# moving back to original cwd #
	chdir($cwd) or die $!;
		
	$pm->finish;
	}

$pm->wait_all_children;


### Subroutines
sub sum_distances{
	my ($dist_file, $pops_r, $infile) = @_;
	
	my %pop_pdist;
	open IN, $dist_file or die $!;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @line = split / /;
		
		for my $i (0..1){
			my @tmp = split /$delim/, $line[$i], 2;
			$line[$i] = $tmp[0];
			}
		
		# checking for pop of each taxon (if pops designated) #
		if($pops_r){
			die " ERROR: $dist_file -> $line[0] not found in population file!\n"
				unless exists $pops_r->{$line[0]};
			die " ERROR: $dist_file -> $line[1] not found in population file!\n"
				unless exists $pops_r->{$line[1]};
		
			# loading pdistance into correct population comparison #
			my $pops = join("__", sort{$a cmp $b} ($pops_r->{$line[0]}, $pops_r->{$line[1]}) );
			push @{$pop_pdist{$pops}}, $line[2];
			}
		
		# loading pdistance for total pdistance #
		push @{$pop_pdist{"total"}}, $line[2];
		}
	close IN;
	
	# summing distances #
	foreach my $pop (keys %pop_pdist){
		#$pop_pdist_sum{$pop} = ( sum(@{$pop_pdist{$pop}}) / scalar @{$pop_pdist{$pop}} ) * 100;
		print join("\t", $infile, $pop, 
				( sum(@{$pop_pdist{$pop}}) / scalar @{$pop_pdist{$pop}} ) * 100 ), "\n";
		}
		
	}

sub load_files{
	my @infiles;
	while(<>){
		chomp;
		next if /^\s*$/;
		push @infiles, $_;
		}
	return \@infiles;
	}

sub call_mother{
# calling mothur to get pairwise distances among all sequences #
	my ($infile, $mothur_cmd, $procs) = @_;
	$mothur_cmd =~ s/\?/$infile/g;
	$mothur_cmd =~ s/\)/, processors=$procs)/;
		#print Dumper $mothur_cmd; exit;
	`mothur "#$mothur_cmd"`;
	(my $dist_file = $infile) =~ s/\.[^.]+$|$/.dist/;
	return $dist_file;
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
		$pops{$line[0]} = $line[1];
		}
	close IN;
		#print Dumper %pops; exit;
	return \%pops;
	}
	

__END__

=pod

=head1 NAME

percent_id_byPop.pl -- Calculating average percent ID within & between populations

=head1 SYNOPSIS

=head2 Piping in a list of files

find . -name "*fasta" | percent_id_byPop.pl [options] > percentID.txt

=head2 Providing files as arguments

percent_id_byPop.pl [options] file(s).fasta > percentID.txt

=head2 options

=over

=item -population

Population file (2-column, tab-delimited, no header). 1st=taxon_name 2nd=population_ID

=item -delimiter

Delimiter separating taxon_name from annotation in sequence files. [" "]

=item -fork

Number of files to process in parallel. [1]

=item -processors 

Number of processors used by Mothur. [1]

=item -mothur

Mothur cmd used to produce a pdistance matrix. ["dist.seqs(fasta=?, calc=onegap, countends=T)"]

=item -h	This help message

=back

=head2 For more information:

perldoc percent_id_byPop.pl

=head1 DESCRIPTION

Get average percent sequence identity of all comparisons of sequences.
Average percent identity comparisons are parsed by population
(both within and between) if a population file is provided. 

=head2 Required population table columns:

=over

=item 1) taxon name (must match sequence names!)

=item 2) population number (1,2,3,etc.)

=back

Not all taxa in the population file need to be found in each alignment.
The '-delimiter' option can be used is annotations or other info is 
appended on the taxon names in the sequence files.

=head1 EXAMPLES

=head2 Just average percent identity among all sequences

percent_id_byPop.pl < *.fasta > percent-ID.txt

=head2 Comparing populations

percent_id_byPop.pl -p pops.txt < *.fasta > percent-ID-pop.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

