#!/usr/bin/perl 

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell


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
my ($format_in, $format_out, $verbose, $names_in);
GetOptions(
	   "in=s" => \$format_in,
	   "out=s" => \$format_out,
	   "name=s" => \$names_in, 
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
	
# loading names file if provided #
my $names_r = load_names($names_in) if defined $names_in;
	
# writing sequences #
while (my $seq = $seqin->next_seq){	
	# changing description if genbank #
	if($format_in =~ /genbank/i){			# genbank->fasta: using LOCUS name
		$seq->desc($seq->display_id);
		}
	# filtering by names list #
	next if defined $names_r && ! exists $names_r->{$seq->primary_id};

	# writing out sequence #	
	$seqout->write_seq($seq);
	}


#--- Subroutines ---#
sub load_names{
	my ($names_in) = @_;
	die "ERROR: cannot find $names_in!\n" 
		unless -e $names_in || -l $names_in;
	
	open IN, $names_in or die $!;
	my %names;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		$names{$_}++;	
		}
	close IN;
	
	# checking for > 1 #
	foreach my $name (keys %names){
		print STDERR "WARNING: $name found $names{$name} times in the names file!\n"
			if $names{$name} > 1;
		}
	
	return \%names;
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

=item -n 	File of sequence names to keep (1 name per line)

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

=head2 Parse/keep just some sequences in a fasta

seqIO.pl -i fasta -o fasta -n seqs2keep.txt < file.fasta > file_subset.fasta

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

