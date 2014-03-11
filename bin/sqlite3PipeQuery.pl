#!/usr/bin/env perl

=pod

=head1 NAME

sqlite3PipeQuery.pl -- script to provide piping sql queries into an sqlite3 DB

=head1 SYNOPSIS

echo "select * from my_table" | sqlite3PipeQuery.pl > results.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -h	This help message

=back

=head2 For more information:

perldoc sqlite3PipeQuery.pl

=head1 DESCRIPTION

Simple script for piping sql queries 
and getting a table of results.

=head2 Output

tab-delimited table

=head1 EXAMPLES

=head2 Write all spacers to a fasta:

echo "select * from my_table" | sqlite3PipeQuery.pl -d my_db.sqlite | less

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

#--- modules ---#
# core #
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;


#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
GetOptions(
	   "database=s" => \$database_file,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
die "ERROR: provide a database file!\n"
  unless defined $database_file;
die "ERROR: cannot find $database_file\n"
  unless -e $database_file;

#--- MAIN ---#
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr)
  or die " Can't connect to $database_file!\n";

my @query = <>;
my $query = join(" ", @query);

my $ret = $dbh->selectall_arrayref($query);

foreach my $row (@$ret){
  map{ $_ = 'NULL' unless defined $_ } @$row;
  print join("\t", @$row), "\n";
}


#--- disconnect ---#
$dbh->disconnect();
exit;





