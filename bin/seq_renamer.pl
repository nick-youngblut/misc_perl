#!/usr/bin/env perl

### packages/perl_flags
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use Spreadsheet::ParseExcel;

### global variables
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

### I/O
my ($namefile, $verbose, $outfile, $sort_bool);
GetOptions(
	   "name=s" => \$namefile,
	   "sort" => \$sort_bool,
	   "verbose" => \$verbose,
	   "outfile=s" => \$outfile,
	   "help|?" => \&pod2usage # Help
	   );

### Input error check
die " ERROR: Provide a renamer list; 2 columns (old_names, new_names)\n"
	unless $namefile;


### Main
my ($fasta_r, $seq_order_r) = load_fasta();
my $namelist_r = load_namefile($namefile);
change_names($fasta_r, $seq_order_r, $namelist_r);

#----------------------Subroutines----------------------#
sub load_fasta{
	# version: 2.0
	# usage: load_fasta($fasta_file_name); returns a hash of sequences
	# IN($fasta_name)
	# OUT($%)
	my(%fasta, $tmpkey, $tmpval);
	my @seq_order;
	while(<>){
		chomp;		
		next if /^\s*$/;
		
		if ($_ =~ />.+/){
			unless (! $tmpval){
				$tmpkey =~ s/^>//;
				$fasta{$tmpkey} = $tmpval; #loading first seq into hash
				push(@seq_order, $tmpkey);
				}
			$tmpkey = $_;
			$tmpval = ""
			}
		else{
			$tmpval = $tmpval . $_; #combining lines of sequence
			}
		}
	$tmpkey =~ s/^>//;
	$fasta{$tmpkey} = $tmpval; #loading last sequence
	#print Dumper(%fasta); exit;
	return \%fasta, \@seq_order;
	} 

sub load_namefile{
	# version: 0.1
	# usage: loading name file (.txt, .xls, .xlsx)
	# output: %
	my $namefile = shift;
	## checking file extension
	my (%namelist, %newnamechk);
	if($namefile =~ /.txt$/i){
		open(FILE, $namefile) or die $!;
		while(<FILE>){
			chomp;
			$_ =~ s/#.+//;
			if($_ eq ""){ next;}
			my @tmp = split(/\t/);
			die " ERROR: '$_' does not have >=2 columns!\n"
				unless scalar @tmp >= 2;
			if(exists($namelist{$tmp[0]})){ die "name file contains duplicate old names: $tmp[0]\n!";}
			else{ $namelist{$tmp[0]} = $tmp[1]; }
			$newnamechk{$tmp[1]} = $tmp[0];
			}
		}
	elsif($namefile =~ /.xls$/i){
		my $parser = Spreadsheet::ParseExcel->new();
		my $workbook = $parser->parse($namefile);
		if(!defined $workbook){ die $parser->error(), ".\n"; }
		for my $worksheet ($workbook->worksheets()){
			my ($row_min, $row_max) = $worksheet->row_range();
			my ($col_min, $col_max) = $worksheet->col_range();
			if($col_max==0 || $row_max==0){ die "ERROR: name list appears empty\n";}
			for my $row ($row_min .. $row_max){
				my $cell1 = $worksheet->get_cell( $row, 0 );
				my $cell2 = $worksheet->get_cell( $row, 1 );
				next unless $cell1;
				next unless $cell2;
				if(exists($namelist{$cell1->value()})){ die "name file contains duplicate old names!";}
				else{ 
					$namelist{$cell1->value()} = $cell2->value();
					}
				$newnamechk{$cell2->value()} = $cell1->value();
				}
			last;
			}
		}
	else{ #if no file extension, default: .txt
		open(FILE, $namefile) or die $!;
		while(<FILE>){
			chomp;
			$_ =~ s/#.+//;
			if($_ eq ""){ next;}
			my @tmp = split(/\t/);
			if(exists($namelist{$tmp[0]})){ die "name file contains duplicate old names!";}
			else{ $namelist{$tmp[0]} = $tmp[1]; }
			$newnamechk{$tmp[1]} = $tmp[0]; 
			}		
		}
		#print Dumper(%namelist); exit;
	return \%namelist;
	}

sub change_names{
	# usage: changing names
	my($fasta_r, $seq_order_r, $namelist_r) = @_;
	
	#foreach (sort keys %$fasta_r){
	foreach (@$seq_order_r){
		die " ERROR: \"$_\" not found in rename list!\n"
			unless exists $namelist_r->{$_};
		print join("\n", ">$namelist_r->{$_}", $fasta_r->{$_}), "\n"
			unless $namelist_r->{$_} =~ /^delete$/i;
		}
	}


__END__

=pod

=head1 NAME

seq_renamer.pl -- rename/delete sequence in a fasta file

=head1 SYNOPSIS

seq_renamer.pl [options] < file.fasta > file_renamed.fasta

=head2 options

=over

=item -n 	Name file (2 column: old_name, new_name).

=item -h	This help message

=back

=head2 For more information:

perldoc seq_renamer.pl

=head1 DESCRIPTION

Batch rename/delete sequences in a fasta file.

The name file can be .txt (tab-delimited) or .xls (excel).

.xls format: 1st worksheet, header required.


=head1 EXAMPLES

=head2 Basic usage:

seq_renamer.pl -n rename_list.txt < file.fasta > renamed.fasta

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

