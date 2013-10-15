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

my $delim = "\t";
my $regex;
my $replace = "";
GetOptions(
	   "regex=s" => \$regex,
	   "replace=s" => \$replace,
	   "delimiter=s" => \$delim, 
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a file name!\n" unless $ARGV[0];
my $filename = $ARGV[0];
$regex = qr/$regex/ if $regex;
$filename =~ s/$regex/$replace/g if $regex;

### MAIN
while(<>){
	if(/^\s*$/){
		print;
		next;
		}
	print join($delim, $filename, $_);
	}


__END__

=pod

=head1 NAME

addFileName.pl -- add file name to each line of a table (good w/ xargs)

=head1 SYNOPSIS

addFileName.pl [options] FILE_NAME < FILE > FILE_name-added

=head2 options

=over

=item -delim  <char>

Column delimiter in table. ['\t']

=item -regex  <char>

Regex for editing file name. ['']

=item -replace  <char>

Replacement for regex matches. ['']

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc addFileName.pl

=head1 DESCRIPTION

Quick script for adding file name (or part of file name)
as the 1st colum in the provided table file.
Useful w/ xargs. 

Blank lines will be skipped.

=head1 EXAMPLES

=head2 xargs

find ./blastp_tables/ -name "*blast.txt" | xargs -n 1 -I % bash -c "addFileName.pl < % % -reg '.+maxbit_|\.blast.+' " | less

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

