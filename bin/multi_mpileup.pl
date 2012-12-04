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

my ($verbose, $fork, @ass_list, @fwd_list, @rev_list);
GetOptions(
	   "assembly=s{,}" => \@ass_list,
	   "forward=s{,}" => \@fwd_list,
	   "reverse=s{,}" => \@rev_list,
	   "processes=s" => \$fork,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide 1 or more assemblies (fasta format), a foward read file (fastq), & a reverse file file (fastq)\n"
	if ! @ass_list || ! @fwd_list || ! @rev_list;
die " ERROR: the number of read files does not match the number of assembly files\n"
	if $#ass_list != $#fwd_list || $#ass_list != $#rev_list;
@ass_list = check_files(@ass_list);
@fwd_list = check_files(@fwd_list);
@rev_list = check_files(@rev_list);

$fork = 1 if ! $fork;


### MAIN
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

__END__

=pod

=head1 NAME

a5_pipeline.pl -- Assemble isolate genomes from Illumina data with ease

=head1 SYNOPSIS

multi_mpileup.pl -a *assembly.fna -f *read_F.fastq -r *read_R.fastq

=head2 options

=over

=item -a 	Assembly files (fasta)

=item -f 	Forward read files (fastq)

=item -r 	Reverse read files (fastq)

=item -p 	Number of forked processes [1].

=item -v 	Verbose output

=item -h 	This help message

=back

=head2 For more information:

perldoc multi_mpileup.pl

=head1 DESCRIPTION

Run bowtie2 followed by mpileup in parallel with multiple assemblies.
The input is a list of assembly files (fasta format) and their corresponding
read files (forward and reverse files in fastq format).

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

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

