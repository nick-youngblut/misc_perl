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

my ($verbose);
my $datatype = "DNA";
GetOptions(
	   "datatype=s" => \$datatype,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: datatype '$datatype' not recognized (DNA|RNA|PROTEIN)!\n"
	unless $datatype =~ /^(DNA|RNA|PROTEIN)$/i;
	

### MAIN
my $fasta_r = load_fasta();
my $len = get_seq_len($fasta_r);
fasta2newnexus($fasta_r, $datatype, $len);

### Subroutines
sub fasta2newnexus{
# writing 'new' nexus file (includes taxa block #
	my ($fasta_r, $datatype, $len) = @_;
	
	# writing header #
	print "#NEXUS\n\n";
	
	# taxa block #
	my @taxa = keys %$fasta_r;
	print "BEGIN taxa;\n";
	print "\tDIMENSIONS ntax=", scalar @taxa, ";\n";
	print "TAXLABELS\n";
	for my $i (0..$#taxa){
		print "[$i]\t$taxa[$i]\n";
		}
	print ";\nEND;\n";
	
	# characters block #
	print "BEGIN characters;\n";
	print "\tDIMENSIONS nchar=$len;\n";
	print "\tFORMAT\n";
	print "\t\tdatatype=$datatype\n";
	print "\t\tmissing=?\n";
	print "\t\tgap=-\n";
	if($datatype =~ /DNA/i){
		print "\t\tsymbols=\"A C G T\"\n";
		}
	elsif($datatype =~ /RNA/i){
		print "\t\tsymbols=\"A C G U\"\n";	
		}
	elsif($datatype =~ /PROTEIN/i){
		print "\t\tsymbols=\"A R N D C Q E G H I L K M F P S T W Y V Z\"\n";
		}
	else{ die " LOGIC ERROR! $!\n"; }
	print "\t\tlabels\n";
	print "\t\tno transpose\n";
	print "\t\tno interleave\n";
	print "\t;\n";
	
	# Matrix #
	print "\tMATRIX\n";
	foreach my $taxon (keys %$fasta_r){
		print join("\t\t", $taxon, $fasta_r->{$taxon}), "\n";
		}
	print ";\nEND;\n\n";
	}

sub get_seq_len{
	my ($fasta_r) = @_;
	my $len = 0;
	foreach my $taxon (keys %$fasta_r){
		my $new_len = length $fasta_r->{$taxon};
		die " ERROR: sequence lengths are not equal!\n"
			if $len && $len != $new_len;
		$len = $new_len;
		}
	return $len;
	}

sub load_fasta{
	my (%fasta, $tmpkey);
	while(<>){
		chomp;
 		 s/#.+//;
 		next if  /^\s*$/;	
 		if(/>.+/){
 			s/^>//;
 			$fasta{$_} = "";
 			$tmpkey = $_;	# changing key
 			}
 		else{$fasta{$tmpkey} .= $_; }
		}
		#print Dumper %fasta; exit;
	return \%fasta;
	} #end load_fasta



__END__

=pod

=head1 NAME

fasta2NewNexus.pl -- convert a fasta to the 'new' nexus format (readable by SplitsTree)

=head1 SYNOPSIS

fasta2NewNexus.pl [options] < input > output

=head2 options

=over

=item -datatype

'datatype' setting in nexus (DNA|RNA|PROTEIN). [DNA]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc fasta2NewNexus.pl

=head1 DESCRIPTION

Convert a fasta to a nexus format that SplitsTree will accept.
The fasta must contain aligned sequences!

=head1 EXAMPLES

=head2 Usage: DNA alignment

fasta2NewNexus.pl < file.fasta > file.nex

=head2 Usage: Amino Acid alignment

fasta2NewNexus.pl -d PROTEIN < file.fasta > file.nex

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

