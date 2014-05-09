#!/usr/bin/env perl

=pod

=head1 NAME

hhblits-search.pl -- run hhblits then hhsearch on an AA sequence

=head1 SYNOPSIS

echo "MNKIEKDIFKNFVVFLMTPENIKLKTLEMTFEGAG" | hhblits-search.pl [flags] 

=head2 Required flags

=over

=item -prefix  <char>

Prefix for output file. A new directory will be made if file path not found.

=back

=head2 Optional flags

=over

=item -database  <char>

hhblits database ('nr' or 'uniprot'). [nr]

=item -fasta  <char>

Multi-fasta file of AA sequences (instead of piping 1 sequence from STDIN).
Each sequence will produce separate output files. 

=item -cpu  <int>

Number of CPUs to use for hhblits & hhsearch. [1]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc hhblits-search.pl

=head1 DESCRIPTION

Help identify homology to known proteins.
Provide an AA sequence in fasta format
via STDIN (see examples).
A sequence name doesn't have to be provided.

The pipeline is as follows:

=over

=item 1) call hhblits on sequence to produce a a3m file.

=item 2) call hhsearch using a3m file on COG, pfam, and CDD databases.

=back 

=head1 EXAMPLES

=head2 Basic usage:

printf ">test\nMNKIEKDIFKNFVVFLMTPENIKLKTLEMTF" | hhblits-search.pl -pre output

=head2 No name in sequence

printf "MNKIEKDIFKNFVVFLMTPENIKLKTLEMTF" | hhblits-search.pl -pre output

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

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

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b, $outfile, $fasta_in);
my $db = "nr";
my $cpu = 1;
GetOptions(
	   "database=s" => \$db,
	   "prefix=s" => \$outfile,
	   "fasta=s" => \$fasta_in,
	   "cpu=i" => \$cpu,
	   "verbose" => \$verbose_b,
	   "help|?" => \&pod2usage # Help
	  );

#--- I/O error ---#
die "ERROR: provide an output file\n"
	unless defined $outfile;
$outfile = File::Spec->rel2abs($outfile);

#--- MAIN ---#
# making an output directory if not present
my @parts = File::Spec->splitdir($outfile);
mkdir $parts[1] unless -d $parts[1];		# making output directory if not found

# loading stdin
my $fasta_r;
if(defined $fasta_in){
  $fasta_r = load_fasta($fasta_in);
}
else{
  $fasta_r = load_seq();
}

foreach my $seq (@$fasta_r){
  
  # outfile
  my $outfile2;
  if( scalar @$fasta_r > 1){
    $outfile2 = multi_outfile($outfile, $seq);
  }
  else{
  $outfile2 = $outfile;
  }

  # calling hhblits
  my $outa3m = call_hhblits($seq, $outfile2, $db, $cpu);

  # calling hhsearch
  call_hhsearch($outa3m, $outfile2, $cpu);
}


#--- Subroutines ---#
sub multi_outfile{
  my ($outfile, $seq) = @_;

  my @l = split /\\n/, $seq;
  die "ERROR: '$seq' not formatted correctly $!\n"
    unless scalar @l == 2;

  $l[0] =~ s/^>//;
  $l[0] =~ s/[| \/]/_/g;
  $outfile = $outfile . "__$l[0]";

  #print Dumper $outfile; exit;
  return $outfile;
}


sub call_hhsearch{
	my ($outa3m, $outfile, $cpu) = @_;
	
	my %dbs = ("COG" => "/home/gtl-shared-2/hh-suite_db/COG_18Feb11.hhm",
			   "pfam" => "/home/gtl-shared-2/hh-suite_db/pfamA_22Nov13.hhm",
			   "cdd" => "/home/gtl-shared-2/hh-suite_db/cdd_18Feb11.hhm");
			   
	foreach my $db (sort keys %dbs){
		(my $out_db = $outfile) =~ s/$/_$db.hrr/;
		my $cmd = "hhsearch -i $outa3m -d $dbs{$db} -o $out_db -cpu $cpu";
		print STDERR "$cmd\n";
		system($cmd);
		}
	}

sub call_hhblits{
# calling hhblits
	my ($seq, $outfile, $db, $cpu) = @_;
	
	# selecting db
	my $db_cmd;
	if($db =~ /nr/i){
		$db_cmd = "/home/gtl-shared-2/hh-suite_db/nr20/nr20_12Aug11";
		}
	elsif($db =~ /uniprot/i){
		$db_cmd = "/home/gtl-shared-2/hh-suite_db/uniprot20_2013_03/uniprot20_2013_03";
		}
	else{
		die "ERROR: do not recognize database \"$db\"\n";
		}
	
	# outfile (a3m)
	(my $outa3m = $outfile) =~ s/$/.a3m/;
	
	# calling 
	my $cmd = "printf \"$seq\" | hhblits -i stdin -cpu $cpu -d $db_cmd -oa3m $outa3m";
	print STDERR "$cmd\n";
	system($cmd);
	
	return $outa3m; 
	}

sub load_fasta{
  my ($fasta_in) = @_;
  
  open IN, $fasta_in or die $!;
  my (%fasta, $tmpkey);
  while(<IN>){
    chomp;
    next if /^\s*$/;

    if(/^>/){
      $tmpkey = $_;
      $fasta{$tmpkey} = "";
    }
    else{
      $fasta{$tmpkey} .= $_;
    }
  }
    
  # converting to an array
  my @fasta;
  foreach my $name (keys %fasta){
    push @fasta, join("\\n", $name, $fasta{$name});
  }

  #print Dumper @fasta; exit;
  return \@fasta;
}

sub load_seq{
  my $seq;
  while(<>){
    chomp;
    next if /^\s*$/;
    if($.==1 ){
      if(/^>/){
	$seq = "$_\\n";
      }
      else{		# if no sequence name provided, adding name
	$seq = ">sequence\\n$_";
      }
    }
    else{
      $seq .= $_;
    }
  }
  #print Dumper $seq; exit;
  return [$seq];
}

