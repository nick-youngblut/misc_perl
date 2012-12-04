#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Parallel::ForkManager;
use File::Path qw/remove_tree/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $fork, @ass_list, @fwd_list, @rev_list, @mpile_list, $bin);
GetOptions(
	   "mpile=s{,}" => \@mpile_list,		# list of mpileup files
	   "bin=i" => \$bin, 					# bin size for summarizing
	   "assembly=s{,}" => \@ass_list,
	   "forward=s{,}" => \@fwd_list,
	   "reverse=s{,}" => \@rev_list,
	   "processes=s" => \$fork,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
# defaults #	
$fork = 1 if ! $fork;
$bin = 5000 if ! $bin;

# skipping bowtie2 & mpilup if summarizing #
if(@mpile_list){
	print "File\t" if scalar @mpile_list > 1;
	print join("\t", qw/Bin_end_position Average_cov Average_mapping_qual/), "\n";
	foreach my $infile (@mpile_list){
		die " ERROR: $infile not found\n" if ! -e $infile;
		$infile = File::Spec->rel2abs($infile);
		make_mpileup_summary($infile, scalar @mpile_list, $bin);
		}
	exit;
	}

# I/O error
die " ERROR: provide 1 or more assemblies (fasta format), a foward read file (fastq), & a reverse file file (fastq)\n"
	if ! @ass_list || ! @fwd_list || ! @rev_list;
die " ERROR: the number of read files does not match the number of assembly files\n"
	if $#ass_list != $#fwd_list || $#ass_list != $#rev_list;
@ass_list = check_files(@ass_list);
@fwd_list = check_files(@fwd_list);
@rev_list = check_files(@rev_list);


### MAIN
# calling bowtie2 & mpileup #
my $pm = new Parallel::ForkManager($fork);
for my $i (0..$#ass_list){
	my $pid = $pm ->start and next;
	
	# making directory #
	my @ass_parts = File::Spec->splitpath($ass_list[$i]);
	my $newdir = make_dir( \@ass_parts );
	
	# making symlinks of the necessary files #
	symlink($ass_list[$i], "$newdir/$ass_parts[2]") or die $!;
	my @fwd_parts = File::Spec->splitpath($fwd_list[$i]);
	my @rev_parts = File::Spec->splitpath($rev_list[$i]);
	symlink($fwd_list[$i], "$newdir/$fwd_parts[2]") or die $!;
	symlink($rev_list[$i], "$newdir/$rev_parts[2]") or die $!;
	
	# line_wrapping genome file and calling faidx #
	run_command("perl -pe 's/([^>]{100})/\$1\\n/g; s/\\n\\s*\\n/\\n/g' $ass_parts[2] > $ass_parts[2]\_lw.fna");
	run_command("samtools faidx $ass_parts[2]\_lw.fna");
	
	# running bowtie #
	run_command("bowtie2-build $ass_parts[2] $ass_parts[2]\_mapped");
	run_command("bowtie2 -x $ass_parts[2]\_mapped -1 $fwd_parts[2] -2 $rev_parts[2] -S $ass_parts[2]\_mapped.sam");
	
	# making sorted bam file #
	run_command("samtools view -bS $ass_parts[2]\_mapped.sam > $ass_parts[2]\_mapped.bam");
	run_command("samtools sort $ass_parts[2]\_mapped.bam $ass_parts[2]\_mapped.sorted");
	
	# running mpileup #
	run_command("samtools mpileup -B -d 10000 -f $ass_parts[2] $ass_parts[2]\_mapped.sorted.bam > $ass_parts[2].mpile");
	
	$pm->finish;
	}

$pm->wait_all_children;


### Subroutines
sub run_command{
	my $cmd = shift;
	print STDERR "\$ ", $cmd, "\n";
	print STDERR `$cmd`;
	}

sub check_files{
# checking existence of files #
	my @file_list = @_;
	foreach(@file_list){
		die " ERROR: cannot find $_\n" if ! -e $_;
		$_ = File::Spec->rel2abs($_);
		}
	return @file_list;
	}

sub make_dir{
# making directory and changing directory #
	my $ass_parts_ref = shift;

	my $cdir = File::Spec->rel2abs(File::Spec->curdir());
	(my $newdir = $$ass_parts_ref[2]) =~ s/\.[^\.]+$|$/_mpile/;

	$newdir = "$cdir/$newdir";
	remove_tree($newdir) if -d $newdir;			# removing old version of directory
	mkdir $newdir or die $!; 
	chdir $newdir or die $!;

	return $newdir;
	}

sub make_mpileup_summary{
# making summary of mpileup output #
# binning mean coverage and mapping quality #
	my ($infile, $file_cnt, $bin) = @_;
	
	open IN, $infile or die $!;
	
	my(@cov, @qual);
	while(<IN>){
		my @line = split /\t/;
		if($line[1] % $bin == 0 || eof){
			print $infile, "\t" if $file_cnt > 1;
			print join("\t", $line[1], average(\@cov), average(\@qual)), "\n";
			@cov = ();
			@qual = ();
			}
		push(@cov, $line[3]);			# pushing coverage
		
		add_quals($line[4], \@qual)		# pushing mapping quality scores
		}
	}

sub add_quals{
# adding mapping quality scores to summary array #
	my ($mapline, $qual_ref) = @_;
	while($mapline =~ /\^./g){
		(my $tmp = $&) =~ s/\^//;
		$tmp = ord($tmp) -33;
		die " ERROR: qual score->$tmp is outside of range (1-42)\n" if $tmp > 42 || $tmp < 0;
		push(@$qual_ref, $tmp);			# phred+33; highest score appears to be 42
		}
		#print Dumper @$qual_ref;
	}

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

__END__

=pod

=head1 NAME

multi_mpileup.pl -- Parallel runs of bowtie2 and mpileup

=head1 SYNOPSIS

=head2 Parallel runs of bowtie2 and mpileup

multi_mpileup.pl -a *assembly.fna -f *read_F.fastq -r *read_R.fastq

=head2 Summarizing mpileup files

multi_mpileup.pl -m *.mpile > mpile_summary.txt

=head2 options

=over

=item -a 	Assembly files (fasta)

=item -f 	Forward read files (fastq)

=item -r 	Reverse read files (fastq)

=item -p 	Number of forked processes [1].

=item -m 	Mpileup files

=item -b 	Bin size [5000]

=item -v 	Verbose output

=item -h 	This help message

=back

=head2 For more information:

perldoc multi_mpileup.pl

=head1 DESCRIPTION

Run bowtie2 followed by mpileup in parallel with multiple assemblies.
The input is a list of assembly files (fasta format) and their corresponding
read files (forward and reverse files in fastq format). The final output is *.mpile

=head2 Warning

Each list of files is assumed to be in the same order!
Check the commands as they are run (printed to stderr).

=head2 Workflow:

=over

=item 1)	make a directory for each assembly and reads

=item 2)	samtools faidx on assembly file

=item 3)	read mapping with bowtie2

=item 4)	making and sorting a bam file

=item 5)	running mpileup

=back

=head1 EXAMPLES

=head2 Single assembly

multi_mpileup.pl -a Msar_barkeri_fusaro.fna -f Msar_barkeri_fusaro_2m_F.fq -r Msar_barkeri_fusaro_2m_R.fq 

=head2 Multiple assemblies (4 parallel processes)

multi_mpileup.pl -a *.fna -f *_F.fq -r *_R.fq -p 4

=head2 Summarizing mpileup files (bin size of 10000bp)

multi_mplieup.pl -b 10000 -m *.mpile > mpile_summary.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

