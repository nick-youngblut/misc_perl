#!/usr/bin/perl


### packages/perl_flags
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::SeqIO;
use Pod::Usage;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


### global variables
my ($error);

### I/O
my ($format_in, $format_out, $verbose);
GetOptions(
	   "in=s" => \$format_in,
	   "out=s" => \$format_out,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error check
if(! $format_in){
	die "ERROR: Provide the sequence input format.\n (See http://www.bioperl.org/wiki/HOWTO:SeqIO for formats)\n (";
	}
if(! $format_out){
	die "ERROR: Provide the sequence ouput format.\n (See http://www.bioperl.org/wiki/HOWTO:SeqIO for formats)\n (";
	}	


### data processing
my ($seqin, $seqout);
eval{
	$seqin = Bio::SeqIO->new(-fh => \*STDIN,
							-format => $format_in);
	$seqout = Bio::SeqIO->new(-fh => \*STDOUT,
							-format => $format_out);
	};
if($@){		# if error
  	print STDERR "Was not able to open files, sorry!\n";
   	print STDERR "Full error is\n\n$@\n";
   	exit(-1);
	}
while (my $seq = $seqin->next_seq){
	if($format_in =~ /genbank/i){			# genbank->fasta: using LOCUS name
		$seq->desc($seq->display_id);
		}
	$seqout->write_seq($seq);
	}


__END__

=pod

=head1 NAME

seqIO.pl -- use Bio::SeqIO to convert among sequence file formats

=head1 SYNOPSIS

seqIO.pl [flags] < input > output

=head2 Required flags

=over

=item -i 	Input file format.

=item -o 	Output file format.

=back

=head2 Optional flags

=over

=item -h	This help message

=back

=head2 For more information:

perldoc seqIO.pl

=head1 DESCRIPTION

Use Bioperl's Bio::SeqIO module to convert
among various sequence file formats.

 If converting genbank => fasta, LOCUS name used.

=head1 EXAMPLES

=head2 Convert genbank to fasta

seqIO.pl -i genbank -o fasta < file.gbk > file.fasta

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

