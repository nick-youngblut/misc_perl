#!/usr/bin/env perl

=pod

=head1 NAME

template.pl -- script template

=head1 SYNOPSIS

tbl2asn_batch.pl [flags] < input > output

=head2 Required flags

=over

=item -p

Path to all of the necessary file (eg. *fsa & *cmt)

=item -c

Configuration file. Similar to Circos cfg files. 
Use '-x' to write an example to STDOUT.

=back

=head2 Optional flags

Organism name

=over -s

regular expression to pull the strain name from the
organism name. 2 arguments required [-s '.+ ' '']

=over

=item -x 

Write an example config file to STDOUT.

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc tbl2asn_batch.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2

=head1 EXAMPLES

=head2 Basic usage:

tbl2asn_batch.pl < input > output

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
use Config::General 2.50;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b, $config_file);
my ($organism, $path);
my @strain_regex = ('.+ ', '');
GetOptions(
	   "organism=s" => \$organism,
	   "strain=s{2,2}" => \@strain_regex,
	   "path=s" => \$path,
	   "config=s" => \$config_file,
	   "x" => \&write_config,
	   "verbose" => \$verbose_b,
	   "help|?" => \&pod2usage # Help
	  );

#--- I/O error ---#

die "ERROR: provide a config file\n"
  unless defined $config_file;
die "ERROR: provide a path to the necesary files\n"
  unless defined $path;

$path = File::Spec->rel2abs($path);
$strain_regex[0] = qr/$strain_regex[0]/;
(my $strain = $organism) =~ s/$strain_regex[0]/$strain_regex[1]/
  if $organism;

#--- MAIN ---#

# getting configuration options
my $conf = Config::General->new(-SplitPolicy => 'equalsign',
			     -ConfigFile => $config_file,
			     -LowerCaseNames => 0);
my %config = $conf->getall;


# making tbl2asn command
my $cmd = make_tbl2asn_cmd(\%config, $path, $organism, $strain);

# calling tbl2asn
call_tbl2asn($cmd);

#--- Subroutines ---#
sub write_config{
  while(<DATA>){ print; }
  exit;
}

sub call_tbl2asn{
# calling tbl2asn
  my ($cmd) = @_;

  $cmd = "echo | $cmd";
  print STDERR $cmd, "\n";
  system($cmd);
  if ($? == -1) {
    print "ERROR: failed to execute: $!\n";
  }
  elsif ($? & 127) {
    printf "ERROR: child died with signal %d, %s coredump\n",
      ($? & 127),  ($? & 128) ? 'with' : 'without';
  }
  else {
    printf "ERROR: child exited with value %d\n", $? >> 8;
  }
}

sub make_tbl2asn_cmd{
# making tbl2asn commadn 
  my ($config_r, $path, $organism, $strain) = @_;

  # path
  my $cmd = "tbl2asn -p $path";
  # template file
  $cmd .= join(" ",' -t', $config_r->{'-t'})
    if defined $config_r->{'-t'};
  # -y flag
  $cmd .= join(" ", ' -y', $config_r->{'-y'})
    if defined $config_r->{'-y'};
  # gap mode
  if(defined $config_r->{gap_mode}){
    foreach my $key (keys %{$config_r->{gap_mode}}){
      $cmd .= join(" ", " $key", $config_r->{gap_mode}{$key});
    }
  } 
  # other options
  if(defined $config_r->{'other'}){
    foreach my $key (keys %{$config_r->{'other'}}){
      $cmd .= join(" ", " $key", $config_r->{other}{$key});
    }
  }
  # -j options
  if(defined $config_r->{'-j'}){
    my $jflag = " -j \"";
    # fixed options
    if(defined $config_r->{'-j'}{fixed}){
      foreach my $key (keys %{$config_r->{'-j'}{fixed}}){
	$jflag .= join("=", " [$key", "$config_r->{'-j'}{fixed}{$key}]");
       }
    }
    # user input options
    if(defined $config_r->{'-j'}{var}){
      if(defined $config_r->{'-j'}{var}{organism}){
	die "ERROR: organism provided in the config file, but not as a flag for this script\n"
	  unless defined $organism;
	$jflag .=join("=", " [organism=$organism]");
      }
      if(defined $config_r->{'-j'}{var}{strain}){
	die "ERROR: strain provided in the config file, but not as a flag for this script\n"
	  unless defined $strain;
	$jflag .=join("=", " [strain=$strain]");
      }
    }
    $cmd .= " $jflag\"";     
  }

  #print Dumper $cmd; exit;
  return $cmd;
}


__DATA__
# tbl2asn_batch.pl configuration file
## WARNING: use single quotes around strings if quotes needed!

## template file
-t = template.sbt

## -y flag
-y = 'Source DNA available from Rachel J. Whitaker (rwhitaker@life.illinois.edu). Culture available from Rachel J. Whitaker (rwhitaker@life.illinois.edu)'

## gap mode
<gap_mode>
-a = r2k
-l = paired-ends
</gap_mode>

## other options
<other>
-M = n
-Z = discrep
-X = C
</other>

## -j options
### use double quotes in this section
### '?' will be replaced by script arguments
<-j>
 <var>
 organism = ?
 strain = ?
 </var>
 <fixed>
 country = "USA: Oregon"
 collection-date = 22-Aug-2011
 isolation-source = sediment
 </fixed>
</-j>
