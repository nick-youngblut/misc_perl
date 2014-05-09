#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Path;
use File::Spec;
use Parallel::ForkManager;


### I/O
my ($uncomp, $verbose, $keep_bool);
my $fork = 1;
GetOptions(
	   "uncomp" => \$uncomp,
	   "fork=i" => \$fork,
	   "keep" => \$keep_bool,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

# loading file names #
my @infiles;
if(! @ARGV){
	while(<>){
		chomp;
		push(@infiles, $_);
		}
	#die " ERROR: no files designated\n" if ! @ARGV;
	}
else{ @infiles = @ARGV; }
map{ $_ = File::Spec->rel2abs($_) } @infiles;

# checking for file existence #
map{ die " ERROR: cannot find $_\n" unless -e $_ or -d $_ } @infiles;


# splitting files into directory level #
my %file_bylevel;
foreach(@infiles){
	my $lev = 0;
	$lev++ while $_ =~ /\//g;
	 
	if(-d $_){
		push(@{$file_bylevel{$lev}{"dir"}}, $_);
		}
	else{ 
		if(-B $_ && ! $uncomp){
			print STDERR "$_ appears already compressed. Skipping\n";
			next;
			}
		elsif(! -B $_ && $uncomp){
			print STDERR "$_ appears already uncompressed. Skipping\n";
			next;		
			}
		push(@{$file_bylevel{$lev}{"file"}}, $_); 
		}
	}

# forking and compressing 
my $pm;
if($fork){ $pm = new Parallel::ForkManager($fork); }
else{ $pm = new Parallel::ForkManager(scalar(@infiles)); }

foreach my $lev (sort {$b <=> $a} keys %file_bylevel){					# deepest level 1st
	foreach my $type (sort {$b cmp $a} keys %{$file_bylevel{$lev}}){	# files before directories
		foreach ( @{$file_bylevel{$lev}{$type}} ){
			$pm->start and next;
			my $infile;
		
			if(! $uncomp){	# if compressing #
				if(-e $_ or -d $_){
					my @p = File::Spec->splitpath($_);
					chdir($p[0] . $p[1]) or die $!;
					$infile = $p[2];
					}
				else{ die $!; }
		
				(my $out = $infile) =~ s/\/*$/.tar.gz/;
				if($out eq $infile){ die $!; }
				if(-e $out or -d $out){ die "  ERROR: compressed file already exists!\n";}
				else{ `tar -czpf $out $infile`; }
				unlink($infile) or rmtree($infile, 0, 1) unless $keep_bool;
				if(scalar(@infiles) > 1){ print STDERR "  $infile compressed\n"; }
				}
			else{ 	# if uncompressing #
				next if ! -B $_;						# if not compressed
				my @p = File::Spec->splitpath($_);		
				chdir($p[0] . $p[1]) or die $!; 
				my $out = `tar -zxvf $_ 2>&1`;
				die " ERROR:\n$out\n" if ($out =~ /error/i);
				unlink($_) or rmtree($_, 0, 1) unless $keep_bool;
				if(scalar(@infiles) > 1){ print STDERR "  $_ uncompressed\n"; }		
				}
			$pm->finish;
			}
		$pm->wait_all_children;		# waiting before next file type or level	
		}
	}

__END__

=pod

=head1 NAME

comp.pl -- convenience compression script

=head1 SYNOPSIS

comp.pl [flags] file(s)

=back

=head2 Flags

=over

=item -f 	Max number of parallel file compressions. [1]

=item -u 	Uncompress files instead of compress. [FALSE]

=item -k 	Keep original file? [FALSE]

=item -v	Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc comp.pl

=head1 DESCRIPTION

Perl wrapper for parallel calling  of tar and gzip 
(makes and uncompresses *tar.gz or *tgz files).
By default the original file is deleted after
(un)compression has completed. 

=head2 WARNING!

No spaces in names allowed!

=head1 EXAMPLES

=head2 Compression:

comp.pl file.txt

=head2 Uncompression:

comp.pl -u file.tar.gz

=head2 Compression (keep original file):

comp.pl -k file.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut



