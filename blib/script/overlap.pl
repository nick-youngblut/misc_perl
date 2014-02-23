#!/usr/bin/env perl

=pod

=head1 NAME

overlap.pl -- get amount of overlap between (scaffold)-start-end entries

=head1 SYNOPSIS

overlap.pl [flags] query.txt subject.txt

=head2 Required flags

=over

NONE

=back

=head2 Optional flags

=over

=item -pos1  <int>

start, stop, [scaffold] columns in query table (2 or 3 arguments). [1 2]

=item -pos2  <int>

start, stop, [scaffold] columns in subject table (2 or 3 arguments). [1 2]

=item -col1  <char>

output columns from query table. [1-]

=item -col2  <char>

output columns from subject table. [1-]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc overlap.pl

=head1 DESCRIPTION

Find any overlap between a query table of start-stop-[scaffold]
positions and a subject table of start-stop-[scaffold] positions.


=head1 EXAMPLES

=head2 Basic usage:

overlap.pl gtlEnvA5udCFSMzCC_core_m-C16_rDTL.gff ../genomic_islands/Methanosarcina_mazei_C16_IsView_HL.txt -pos1 4 5 1 -pos2 2 3 1

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


#--- modules ---#
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Set::IntervalTree;
use List::Util qw/min max/;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b, @pos1, @pos2);
my ($col1, $col2) = ("1-", "1-");
GetOptions(
	"pos1=i{2,3}" => \@pos1,
	"pos2=i{2,3}" => \@pos2,
	"col1=s" => \$col1,
	"col2=s" => \$col2,
	"verbose" => \$verbose_b,
	"help|?" => \&pod2usage # Help
	);

#--- I/O error ---#
die "ERROR: provide 2 input tables (tab-delimited)\n"
	unless defined @ARGV[0..1];
die "ERROR: cannot find '$ARGV[0]'!\n"
	unless -e $ARGV[0];
die "ERROR: cannot find '$ARGV[1]'!\n"
	unless -e $ARGV[1];

@pos1 = set_default_pos(@pos1);
@pos2 = set_default_pos(@pos2);

#--- MAIN ---#
my $itrees_r = load_itrees($ARGV[1], \@pos2);			# subject
query_itrees($itrees_r, $ARGV[0], \@pos1, \@pos2, $col1, $col2);		# query 

#--- Subroutines ---#
sub query_itrees{
	my ($itrees_r, $infile, $pos1_r, $pos2_r, $col1, $col2) = @_;
	
	open IN, $infile or die $!;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @l = split /\t/;
		die "ERROR: $infile must have >= 2 columns!\n"
			unless scalar @l >= 2;
			
		# checking for start-end columns #
		die "ERROR: cannot find start column in $infile -> line $.!\n"
			unless defined $l[$$pos1_r[0]];
		die "ERROR: cannot find end column in $infile -> line $.!\n"
			unless defined $l[$$pos1_r[1]];
			
		# flipping coords if needed #
		@l[ $$pos1_r[0], $$pos1_r[1] ] = 
			coord_flip($l[$$pos1_r[0]], $l[$$pos1_r[1]]);
		
		# scaffold column #
		my $scaf = "";
		if (defined $$pos1_r[2]){
			die "ERROR: cannot find column '$$pos1_r[2]'\n"
				unless defined $l[$$pos1_r[2]];
			$scaf = $l[$$pos1_r[2]];
			}
			
		# querying itree #
		next unless exists $itrees_r->{$scaf};
		my $ret = $itrees_r->{$scaf}->fetch( $l[$$pos1_r[0]], $l[$$pos1_r[1]] );
		next unless @$ret;
		
		foreach my $r (@$ret){
			# determining amount of overlap #
			## query ##
			my ($qs, $qe) = ($l[$$pos1_r[0]], $l[$$pos1_r[1]]);
			## subject ##
			my ($ss, $se) = ($$r[$$pos2_r[0]], $$r[$$pos2_r[1]]);
			
			my $overlap = min($qe, $se) - max($qs, $ss);
			my $query_frac = $overlap / ($qe - $qs);
			my $sub_frac = $overlap / ($se - $ss);
			
			# col #
			my @col1 = expand_ranges($col1, scalar @l);		# query 
			my @col2 = expand_ranges($col2, scalar @$r);	# subject

			
			# writing output #
			print join("\t", $overlap, $query_frac, $sub_frac,
					@l[ @col1 ],
					@$r[ @col2 ]
					), "\n";
			}
		
		}
	close IN;
	}

sub expand_ranges{
	my ($col, $nrow) = @_;
	
	my @cols;
	foreach my $x ( split /,/, $col ){
		if($x =~ /^\d+-$/){
			my @tmp = split /-/, $x;
			for my $i ($tmp[0]..$nrow){
				push @cols, $i -1;
				}
			}
		elsif($x =~ /^\d+-\d+$/){
			my @tmp = split /-/, $x;
			for my $i ($tmp[0]..$tmp[1]){
				push @cols, $i - 1;
				}
			}
		elsif($x =~ /^\d+$/){
			push @cols, $x - 1;
			}
		else{
			die "ERROR: '$col' not formated correctly!\n";
			}
		}
	
	return sort{$a<=>$b} @cols;
	}

sub load_itrees{
	my ($infile, $pos_r) = @_;
	
	open IN, $infile or die $!;
	my %itrees;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @l = split /\t/;
		die "ERROR: $infile must have >= 2 columns!\n"
			unless scalar @l >= 2;

		# checking for start-end columns #
		die "ERROR: cannot find start column in $infile -> line $.!\n"
			unless defined $l[$$pos_r[0]];
		die "ERROR: cannot find end column in $infile -> line $.!\n"
			unless defined $l[$$pos_r[1]];
			
		# flipping coords if needed #
		@l[ $$pos_r[0], $$pos_r[1] ] = 
			coord_flip($l[$$pos_r[0]], $l[$$pos_r[1]]);
		
		# scaffold column #
		my $scaf = "";
		if (defined $$pos_r[2]){
			die "ERROR: cannot find column '$$pos_r[2]'\n"
				unless defined $l[$$pos_r[2]];
			$scaf = $l[$$pos_r[2]];
			}
		
		$itrees{$scaf} = Set::IntervalTree->new()
			unless exists $itrees{$scaf};
	
		# loading itree #
		$itrees{$scaf}->insert(\@l, $l[$$pos_r[0]], $l[$$pos_r[1]] );
		}
	close IN;
	
	return \%itrees;
	}
	
sub coord_flip{
	my ($start, $end) = @_;
	
	die "ERROR: start value '$start' is not an integer!\n"
		unless $start =~ /^\d+$/;
	die "ERROR: end value '$end' is not an integer!\n"
		unless $end =~ /^\d+$/;
	
	if($start <= $end){ return $start, $end; }
	else{ return $end, $start; }
	}

sub set_default_pos{
	my @pos = @_;
	
	my @df = qw/1 2 3/;
	foreach my $i ( (0,1) ){
		$pos[$i] = $df[$i] unless defined $pos[$i];
		}
	map{ $_ -= 1 if defined $_ } @pos;
	
	return @pos[0..2];
	}



