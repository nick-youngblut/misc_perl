#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use POSIX qw/strftime/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $ntime, $pid, $email);
GetOptions(
	   "step=i" => \$ntime,
	   "proc=s" => \$pid,			# process to run
	   "mail=s" => \$email,			
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a pid\n" if ! $pid;
$ntime = 60 if ! $ntime;

my $out = `ps aux | egrep "^[a-zA-Z0-9]+ +$pid +"`;
die " ERROR: PID '$pid' not found!\n" if ! $out;

print STDERR "PID: $pid\n";
print STDERR "Log every: $ntime seconds\n";

my $start = time();

### MAIN
while(1){
	my $out = `ps aux | egrep "^[a-zA-Z0-9]+ +$pid +"`;
	$out =~ s/ +/\t/g;
	my @out = split /\n|\n/, $out;
	print join("\n", @out), "\n";

	if(! @out){			# pid gone
		# time to completion #
		my $end = time();
		my $duration = $end - $start;
		#my $durm = sprintf("%.2f", $duration / 60);
		#my $durh =  sprintf("%.2f", $durm / 60);
		
		# notification #	
		`echo $0 done | mailx -s "PID: $pid done\nFollowed the job for $duration seconds." $email` if $email;
		print STDERR "$pid job done\nFollowed the job for $duration seconds\n";
		exit; 
		}
	sleep $ntime;
	}

__END__

=pod

=head1 NAME

proc_track.pl -- Track 'top' output for a process (pid)

=head1 SYNOPSIS

proc_track.pl [-s] [-m] -p > proc.log

=head2 options

=over

=item -p 	PID (process ID). Use 'top' to find it.

=item -s 	Seconds between tracking PID. [60]

=item -m 	Email for a notification when the PID is finished.

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc proc_track.pl

=head1 DESCRIPTION

Basically log continual 'ps aux' 

=head1 EXAMPLES

=head2 with email reminder

proc_track.pl -p 1214 -m nyoungb2@illinois.edu > pid1214.log

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

