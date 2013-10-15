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

my ($verbose, $header_bool, $keep_meth);
my $length_cutoff = 0.75;
my $identity_cutoff = 50;
GetOptions(
		"length=f" => \$length_cutoff,
		"keep" => \$keep_meth,
		"x" => \$header_bool,
		"identity=f" => \$identity_cutoff,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
foreach my $infile (@ARGV){
	die " ERROR: file not found!\n" unless -e $infile;

	my $blast_res_r = load_blast_table($infile, $length_cutoff, $identity_cutoff);	
	find_top_hit($infile, $blast_res_r, $keep_meth);
	}

sub find_top_hit{
	my ($infile, $blast_res_r, $keep_meth) = @_;
	
	my $top_hit;
	my $top_meth_hit;
	foreach my $row (@$blast_res_r){
		if(! $top_hit){				# top hit
			$top_hit = $row;
			}
		if(! $keep_meth && ! $top_meth_hit && $$row[11] =~ /^methano/i){		# top methanogen hit
			$top_meth_hit = $row;
			}
		
		next if ($top_hit && $top_meth_hit) || ($top_hit && ! $keep_meth);
		}
	
	# no hits #
	return unless $top_hit;
	
	# checking if meth hit is same as top hit #
	my $same = 0;												# if top hit not methano
	$same = 1 if $$top_hit[$#$top_hit] =~ /^methano/i;			# if top hit methano
	$same = 2 if $same == 1 && $$top_meth_hit[$#$top_meth_hit] =~ /^methano/i;		# if both are methano
	
	print join("\t", $infile, "top_hit", $same, pop @$top_hit, @$top_hit), "\n";
	print join("\t", $infile, "top_meth_hit", $same, pop @$top_meth_hit, @$top_meth_hit), "\n" if $top_meth_hit;
	
	}

sub load_blast_table{
# loading blast table from blastxml parse #
	my ($infile, $length_cutoff, $identity_cutoff) = @_;

	open IN, $infile or die $!;	
	my @blast_res;
	while(<IN>){
		next if $. == 1 && ! $header_bool;
		chomp;
		
		# filtering #
		my @line = split /\t/;
		my $q_len = abs($line[8] - $line[7]);
		my $hit_len = abs($line[10] - $line[9]);
		next if $hit_len / $q_len < $length_cutoff;
		next if $line[3] < $identity_cutoff;
		
		# getting taxa #
		my @taxa = split /^.*?\[|\].+\[|\].*$/, $line[2];
		
		if($taxa[2]){ push @line, $taxa[2]; }
		elsif($taxa[1]){ push @line, $taxa[1]; }
		else{ push @line, "NA"; }
		
		push(@blast_res, \@line);
		}
	close IN;

	# sorting by percent ID #
	@blast_res = sort{$$b[3] <=> $$a[3]} @blast_res;
		
		#print Dumper @blast_res; exit;
	return \@blast_res;
	}



__END__

=pod

=head1 NAME

blastxml_parse_top_hit.pl -- get top blast hits (& taxon hit) from blastxml_parse.pl output

=head1 SYNOPSIS

blastxml_parse_top_hit.pl [options] *blast.txt

=head2 options

=over

=item -l 	Hit length cutoff (hit length / query length). [0.75]

=item -i 	Sequence identity cutoff. [50]

=item -k	Keep top /methano/ hit? [TRUE]

=item -x 	Header in blast input table? [TRUE]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc blastxml_parse_top_hit.pl

=head1 DESCRIPTION

Summarizing output from blastxml_parse.pl.
Just finding the top hit and the top methanogen
hit. 

The top methanogen hit is just found with regex '/methano/i',
so it could it other taxa!

=head1 EXAMPLES

=head2 Usage method 1

blastxml_parse_top_hit.pl *blast.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

