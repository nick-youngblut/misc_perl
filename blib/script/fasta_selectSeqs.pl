#!/usr/bin/env perl

=pod

=head1 NAME

fasta_selectSeqs.pl -- selecting sequences from a fasta

=head1 SYNOPSIS

fasta_selectSeqs.pl [flags] > filtered.fasta

=head2 Required flags

=over

=item -fasta  <char>

Fasta file name. ('-' = STDIN)

=item list  <char>

List of sequence names (no '>'; must be extact match). ('-' = STDIN)

=back

=head2 Optional flags

=over

=item -verbose  <bool>

Verbose output. [TRUE]

=item -h  <bool>

Print this help message & exit. [FALSE]

=back

=head2 For more information:

perldoc fasta_selectSeqs.pl

=head1 DESCRIPTION

Simple script for selecting sequences
of interest from a fasta.

Only the sequences listed in '-list'
will be written to STDOUT.

=head1 EXAMPLES

=head2 Basic usage:

fasta_selectSeqs.pl -f in.fasta -l target_seqs.txt > in_filtered.fasta

=head1 AUTHOR

Nick Youngblut <ndy2@cornell.edu>

=head1 AVAILABILITY

email me

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

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b);
my ($fasta_in, $list_in);
GetOptions(
	   "fasta=s" => \$fasta_in,
	   "list=s" => \$list_in,
	   "verbose" => \$verbose_b,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error ---#
die "ERROR: provide a fasta file\n"
  unless defined $fasta_in;
die "ERROR: provide a list of sequences\n"
  unless defined $list_in;
die "ERROR: -fasta == -list!\n"
  if $fasta_in eq $list_in;

#--- setting defaults ---#

#--- MAIN ---#
my $list_r = load_list($list_in);
select_from_fasta($fasta_in, $list_r);

#--- Subroutines ---#
sub select_from_fasta{
# selecting sequences from a fasta
  my ($fasta_in, $list_r) = @_;

  my $fh;
  if($fasta_in eq "-"){
    $fh = *STDIN;
  }
  else{
    open $fh, $fasta_in or die $!;
  }

  my %fasta;
  my $tmpkey;
  while(<$fh>){
    chomp;
    next if /^\s*$/;

    if(eof($fh)){
      $fasta{$tmpkey} .= "";
      print join("\n", ">$tmpkey", $fasta{$tmpkey}), "\n"
	if exists $list_r->{$tmpkey};
    }
    elsif(/^>/){
      if(defined $tmpkey){
	print join("\n", ">$tmpkey", $fasta{$tmpkey}), "\n"
	    if exists $list_r->{$tmpkey};
	%fasta = ();
      }

      ($tmpkey = $_) =~ s/^>//;
      $fasta{$tmpkey} = "";
    }
    else{
      $fasta{$tmpkey} .= $_;
    }
  }
  close $fh or die $!;

}

sub load_list{
# loading list of sequences to keep
  my ($list_in) = @_;

  my $fh;
  if($list_in eq "-"){
    $fh = *STDIN;
  }
  else{
    open $fh, $list_in or die $!;
  }

  my %list;
  while(<$fh>){
    chomp;
    next if /^\s*$/;
    print "WARNING: '$_' is not unique in sequence list\n"
      if !$verbose_b && exists $list{$_};
    $list{$_} = 1;
  }
  close $fh or die $!;

  return \%list;
}
