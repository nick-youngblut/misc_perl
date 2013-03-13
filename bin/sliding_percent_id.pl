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

my ($format_in, $verbose, $total_bool);
my $window = 100;
my $jump = 100;
GetOptions(
	   "in=s" => \$format_in,
	   "window=i" => \$window, 		# window size
	   "jump=i" => \$jump,			# jump size
	   "total" => \$total_bool,		# total percent ID? [FALSE]
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
if(! $format_in){
	print STDERR "SeqIO is guessing input file format.\n (See http://www.bioperl.org/wiki/HOWTO:AlignIO_and_SimpleAlign for formats)\n";
	}


### MAIN
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

while (my $aln = $alnin->next_aln){
	print STDERR "...calculating percent ID\n";
	
	write_total_percent_id($aln) if $total_bool;

	window_percent_id($aln, $window, $jump);
	}

### Subroutines
sub write_total_percent_id{
# writing percent ID for total alignment #
	my ($aln) = @_;
	print join("\t", 0, $aln->length, $aln->percentage_identity), "\n";
	}

sub window_percent_id{
# getting percent ID for a window #
	my ($aln, $window, $jump) = @_;
	
	for (my $i=1; $i<=($aln->length -1); $i+=$jump){
		my $sub_aln = $aln->slice($i, $i + $window - 1);		# subalignment
		my $pID = $sub_aln->percentage_identity;				# percent ID
		$pID = "NA" if ! $pID;
		print join("\t", $i, $i+$window -1, $pID), "\n";		# writing
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

=item -t 	Get percent ID of total alignment length? [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc sliding_percent_id.pl

=head1 DESCRIPTION

Get percent sequence identity for a sequence alignment over a sliding scale.

Just using Bio::AlignIO subroutines.

=head1 EXAMPLES

=head2 Usage

sliding_percent_id.pl -i fasta < file.fna > file_percent-ID.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

