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

my ($verbose, $ntime, $pid, $email, $cmd, $outfile);
GetOptions(
	   "step=i" => \$ntime,
	   "proc=s" => \$pid,			# pid to track
	   "mail=s" => \$email,			
	   "command=s" => \$cmd, 		# command
	   "output=s" => \$outfile,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a pid or a command ('-p' or '-c')\n" if (! $pid && ! $cmd) || ($pid && $cmd);
$ntime = 60 if ! $ntime;

if($pid){
	my $out = `ps aux | egrep "^[a-zA-Z0-9]+ +$pid +"`;
	die " ERROR: PID '$pid' not found!\n" if ! $out;

	print STDERR "PID: $pid\n";
	print STDERR "Log every: $ntime seconds\n";
	}

my $start = time();

### MAIN
if(! $pid){		# if running command
	my @pids;	
	my $pid = fork();
	if($pid){		# if parent
		push(@pids, $pid);
		follow_pid($pid, $outfile);
		}
	else{			# if child
		`$cmd`;
		exit;
		}
	foreach my $child (@pids){
		my $pid = waitpid($child, 0);
		#print "$pid finished\n";
		}
	}
else{		# if following a pid provided pid
	follow_pid($pid, $outfile);
	}

# follow pid #
sub follow_pid{
	my ($pid, $outfile) = @_;
	
	$outfile = "PID$pid.log" if ! $outfile;
	open OUT, ">>$outfile" or die $!;
	
	while(1){
		my $out = `ps aux | egrep "^[a-zA-Z0-9]+ +$pid +"`;
		$out =~ s/ +/\t/g;
		my @out = split /\n|\n/, $out;
	
		if(! @out || $out[0] =~ /<defunct>$/){			# pid gone
			# time to completion #
			my $end = time();
			my $duration = $end - $start;

			# notification #	
			`echo "PID: $pid done\nFollowed the job for $duration seconds." | mailx -s "$0 job done" $email` if $email;
			print STDERR "PID: $pid job done\nFollowed the job for $duration seconds\n";
			last;
			}

		print OUT join("\n", @out), "\n";	
		
		sleep $ntime;
		}
	
	close OUT;
	}
	

__END__

=pod

=head1 NAME

proc_track.pl -- Tracking of 'ps aux' for a process.

=head1 SYNOPSIS

=head2 Provide a command

proc_track.pl [-s] [-m] [-o] -c

=head2 Provide a PID

proc_track.pl [-s] [-m] [-o] -p

=head2 options

=over

=item -c 	Command to execute and follow.

=item -p 	PID (process ID). Use 'top' to find it.

=item -s 	Seconds between tracking PID. [60]

=item -m 	Email for a notification when the PID is finished.

=item -o 	Output file name. [PID#.log]

=item -h	This help message

=back

=head2 For more information:

perldoc proc_track.pl

=head1 DESCRIPTION

Continually log 'ps aux' for a process. If a command is executed,
this script will fork a child process and execute the command; the
parent will then follwo the child process.

Output is appended to old output if file already exists.

=head1 EXAMPLES

=head2 Providing a command

proc_track.pl -m nyoungb2@illinois.edu -s 10 -c "sleep 60"

=head2 Providing a PID

proc_track.pl -m nyoungb2@illinois.edu -p 1214

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

