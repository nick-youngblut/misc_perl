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

my ($verbose, $total_bool, $sort_opt, $prefix);
$prefix = "aln_var";
$sort_opt = "char";
GetOptions(
	   "total" => \$total_bool,
	   "prefix=s" => \$prefix,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults


### MAIN
# loading input #
print STDERR "...loading fasta files\n";
open my $ofh, ">$prefix\_by-pos.txt" or die $!;
my %stats;
my %var;
foreach my $infile (@ARGV){
	die " ERROR: $infile not found!" unless -e $infile;
	print STDERR "." if $verbose;
	
	my ($char_r, $gap_r, $total, $len) = load_fasta($infile);
	Nchar($infile, $char_r, $total, $total_bool, $len, \%stats, \%var, $ofh);
	Ngap($infile, $gap_r, $total, $total_bool, $len, \%stats, \%var, $ofh);
	
	write_ggplot_table($infile, $char_r, $gap_r, $len, $ofh);
	}

close $ofh;
print STDERR "\n" if $verbose;

# writing #
print STDERR "...writing out tables\n";
sort_N_write(\%var, \%stats, $sort_opt, $prefix);
write_stats(\%stats, $prefix);

### Subroutines
sub write_stats{
	my ($stats_r, $prefix) = @_;
	
	open OUT, ">$prefix\_sum.txt" or die $!;
	
	foreach my $infile (keys %{$stats_r->{"char"}}){
		print OUT join("\t", $infile, $stats_r->{"char"}{$infile},
					$stats_r->{"gap"}{$infile}), "\n";
		}
		
	close OUT;
	}

sub write_ggplot_table{
# writing posistion infor for plotting in ggplot #
	my ($infile, $char_r, $gap_r, $len, $ofh) = @_;
	
	foreach my $i (0..($len - 1)){
		# character #
		if (exists $char_r->{$i}){
			print $ofh join("\t", $infile, 
				sprintf("%.3f", $i / $len * 100), 
				"char", scalar keys %{$char_r->{$i}}), "\n";
			}
		else{
			print $ofh join("\t", $infile, 
				sprintf("%.3f", $i / $len * 100),
				"char", 0), "\n";
			}

		# gaps #
		$gap_r->{$i} = 0 unless exists $gap_r->{$i};	
		print $ofh join("\t", $infile, 
			sprintf("%.3f", $i / $len * 100), 
			"gap", $gap_r->{$i}), "\n";		# file, position, gap, gap_value
		}

	}

sub sort_N_write{
# sorting by gaps or positions #
	my ($var_r, $stats_r, $sort_opt, $prefix) = @_;
	
	# N-char #
	open OUT, ">$prefix\_Nchar.txt" or die $!;
	
	foreach my $infile (sort {$stats_r->{"char"}{$b} <=> $stats_r->{"char"}{$a}}
	 			keys %{$stats_r->{"char"}}){	
	 	print OUT join("\t", "$infile", sprintf("%.2f", $stats_r->{"char"}{$infile}),
	 		 join("|", @{$var_r->{$infile}{"char"}})), "\n";

	 	}
	close OUT;

	# N-gap #
	open OUT, ">$prefix\_Ngap.txt" or die $!;
	foreach my $infile (sort {$stats_r->{"gap"}{$b} <=> $stats_r->{"gap"}{$a}}
	 			keys %{$stats_r->{"char"}}){	
	 	#print join("\t", "char $infile", join("|", @{$var_r->{$infile}{"char"}})), "\n";
	 	print OUT join("\t", "$infile", $stats_r->{"gap"}{$infile}, 
	 		join("|", @{$var_r->{$infile}{"gap"}})), "\n";
	 	}
	close OUT; 
	}

sub Ngap{
# summing number of gaps at each position #
	my ($infile, $gap_r, $total, $total_bool, $len, $stats_r, $var_r) = @_;
	
	my @line;
	my $gap_total = 0;
	foreach my $i (0..($len - 1)){			# each position
		if( exists $gap_r->{$i} ){		# if at least 1 character at the position
			$gap_total += $gap_r->{$i};
			
			if($total_bool){
				push(@line, sprintf("%.1f", $gap_r->{$i} / $total * 100 ));
				}
			else{
				push(@line, sprintf("%.0f", $gap_r->{$i} ));
				}
			}
		else{
			push(@line, 0);
			}	
		}
	
	$var_r->{$infile}{"gap"} = \@line;
	$stats_r->{"gap"}{$infile} = $gap_total / scalar @line;
	}

sub Nchar{
# summing number of characters (except gaps) at each position #
	my ($infile, $char_r, $total, $total_bool, $len, $stats_r, $var_r) = @_;
	
	my @line;
	my $char_total = 0;
	foreach my $i (0..($len - 1)){			# each position
		if( exists $char_r->{$i} ){		# if at least 1 character at the position
			my $N = scalar keys %{$char_r->{$i}};
			$char_total += $N;
			
			if($total_bool){
				push(@line, sprintf("%.1f", $N / $total * 100 ));
				}
			else{
				push(@line, sprintf("%.0f", $N ));
				}
			}
		else{
			push(@line, 0);
			}
		}
		
	$var_r->{$infile}{"char"} = \@line;
	$stats_r->{"char"}{$infile} = $char_total / scalar @line;
	
		#print Dumper %$stats_r; 
		#print Dumper %$var_r; exit;
	#print join("\t", $infile, join("|", @line)), "\n";
	}

sub load_fasta{
	# version: 2.0
	# usage: load_fasta($fasta_file_name); returns a hash of sequences
	my ($fasta_in) = @_;
	open IN, $fasta_in or die $!;
	my (%char, %gap, $total, $len);
	my $fasta;
	while(<IN>){
		chomp;
 		$_ =~ s/#.+//;
 		next if  $_ =~ /^\s*$/;	
 		
 		if(/^>/ || eof IN){
 			next unless $fasta;
 			$fasta .= $_ if eof IN;
 			
 			# count positions #
 			my @line = split //, $fasta;
 			for my $i (0..$#line){
	 			if ($line[$i] eq "-"){			# if gap
 					$gap{$i}++;
 					}		# skipping gaps
 				else{
 					$char{$i}{$line[$i]}++;
 					}
 				}
 				
 			# length #
 			die " ERROR: $fasta_in has sequence that is not the same length as the others!\n"
 				if $len && $len != length $fasta;
 			$len = length $fasta;
 			
 			# purging #	
 			$fasta = "";
 			$total++;
 			}
 		else{
 			$fasta .= $_;
 			}
 		
		}
	close IN;

		#print Dumper %char; exit;
	return \%char, \%gap, $total, $len;
	} #end load_fasta


__END__

=pod

=head1 NAME

alignment_var.pl -- determine Character variability and Gap presence in >= 1 alignment

=head1 SYNOPSIS

alignment_var.pl [options] alignment(s).fasta

=head2 options

=over

=item -prefix

Output file prefix. ["aln_var"]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc alignment_var.pl

=head1 DESCRIPTION

Determine the variability in >= 1 alignment. Variability
in terms of number of different characters and
number of gaps at each position.

=head2 Output files:

=head3 *_by-pos.txt

For plotting the variability at each position for each
alignment with ggplot2. The columns are: file, position, char|gap, value

=head3 *_Nchar.txt

Shows number of characters for each position in each alignment.
The file is sorted in decreasing order by most character variability.
Positions are delimited by "|".
Columns are: file, value, position_values

value = average number of characters per position.

=head3 *_Ngap.txt

Shows number of gaps for each position in each alignment.
The file is sorted in decreasing order by most gaps.
Positions are delimited by "|".
Columns are: file, value, position_values

value = number of gaps per position.

=head3 *_sum.txt

Summary character variability and gap values for each alignment.
Columns are: file, ave_num_char_per_position, ave_num_gap_per_position

=head1 EXAMPLES

=head2 Usage: 

alignment_var.pl *fasta 

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

