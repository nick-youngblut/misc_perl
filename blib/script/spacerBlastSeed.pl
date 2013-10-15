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

my ($verbose, $fasta_in, $genome_in);
my $trunc = 8;
my $blast_params = "-evalue 100";
my $pam = "CCN";
GetOptions(
	   "spacers=s" => \$fasta_in,
	   "truncate=i" => \$trunc,			# truncate this
	   "genome=s" => \$genome_in,
	   "blast_params=s" => \$blast_params, 
	   "pam=s" => \$pam,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a fasta of spacers!\n" unless $fasta_in;
die " ERROR: $fasta_in doesn't exist!\n" unless -e $fasta_in;
die " ERROR: provide a fasta of a genome!\n" unless $genome_in;
die " ERROR: $genome_in doesn't exist!\n" unless -e $genome_in;
print STDERR "...PAM used: '$pam'\n";

### MAIN
# load fasta & truncate #
my $fasta_r = load_fasta($fasta_in);
truncate_spacers($fasta_r, $trunc);

# blastn #
my $blast_hits_r = call_blastn($blast_params, $genome_in);

# load genome; screen pams #
my $genome_r = load_fasta($genome_in);
screen_blast_by_mismatch($blast_hits_r);
screen_blast_by_pam($blast_hits_r, $genome_r, $pam);

# load spacer sequence; filter by mismatches #
screen_blast_by_full_length_mismatch($blast_hits_r, $genome_r, $fasta_r);

# writing spacer blasts #
foreach my $hit (keys %$blast_hits_r){
	print join("\t", @{$blast_hits_r->{$hit}}), "\n";
	}


### Subroutines
sub screen_blast_by_full_length_mismatch{
	my ($blast_hits_r, $genome_r, $fasta_r) = @_;
	
	print STDERR "# filtering seed sequence blast by mismatch #\n";
	print STDERR " Number of starting blast hits:\t", scalar keys %$blast_hits_r, "\n";

	foreach my $hit (keys %$blast_hits_r){
		my $sstart = ${$blast_hits_r->{$hit}}[8];
		my $send = ${$blast_hits_r->{$hit}}[9];
		
		my $queryID = ${$blast_hits_r->{$hit}}[0];
		my $subjectID = ${$blast_hits_r->{$hit}}[1];
		
		die " ERROR: $queryID not in full-length fasta!\n"
			unless exists $fasta_r->{$queryID};
		my $spacer_len = length $fasta_r->{$queryID};

		# strand #
		my $strand = "+";
		$strand = "-" if $sstart > $send;

		# get protospacer #
		my $proto;
		if($strand eq "+"){
			$proto = substr($genome_r->{$subjectID}, $sstart -1, $spacer_len);
			}
		elsif($strand eq "-"){
			$proto = substr($genome_r->{$subjectID}, $sstart - $spacer_len, $spacer_len);
			$proto = revcomp($proto) if $strand eq "-";		
			}
		else{ die " LOGIC ERROR: $!\n"; }

		
		# write spacer-protospacer #
		open OUT, ">spacer-proto_pre-align.fasta" or die $!;
		print OUT ">$queryID\_spacer\n$fasta_r->{$queryID}\n"; 		# spacer
		print OUT ">$queryID\_proto\n$proto\n";
		close OUT;
		
		# call mafft #
		my $cmd = "mafft --quiet spacer-proto_pre-align.fasta |";
		open PIPE, $cmd or die $!;
		
		my @seqs;
		my $seq_cnt = -1; 
		while(<PIPE>){
			chomp;
			if(/^>/){
				$seq_cnt++;
				next;
				}
			else{
				push(@{$seqs[$seq_cnt]}, split //); 
				}
			}
		close PIPE;
		die " ERROR: too many sequences!\n" if scalar keys @seqs > 2;
		
		# screen by mismatches #
			#print "strand: $strand\n";
		my $mismatch = count_mismatches(\@seqs, $queryID);
		delete $blast_hits_r->{$hit} if $mismatch > 9;
		}
	
	print STDERR " Number of protospacers (beyond seed) w/ <9 mismatches:\t\t", scalar keys %$blast_hits_r, "\n";
	}

sub count_mismatches{
	my ($seqs_r, $queryID) = @_;
	
	unless ($$seqs_r[0] && $$seqs_r[1] && (scalar @{$$seqs_r[0]} == scalar @{$$seqs_r[1]})){
		print " ERROR $!\n";
		print join("",  @{$$seqs_r[0]}), "\n";
		print join("",  @{$$seqs_r[1]}), "\n";		
		exit;
		}
	
	my $mis_cnt = 0;
	for my $i (($trunc - 1)..$#{$$seqs_r[0]}){		# gaps count as mismatches
		$mis_cnt++ if $$seqs_r[0][$i] ne $$seqs_r[1][$i];
		}
			
	if($verbose){
		print STDERR "spacer: $queryID, N-mismatches: $mis_cnt\n";
		print STDERR join("",  @{$$seqs_r[0]}), "\n";
		print STDERR join("",  @{$$seqs_r[1]}), "\n";	
		}
		
	return $mis_cnt;
	}

sub screen_blast_by_mismatch{
	my ($blast_hits_r) = @_;

	print STDERR "# filtering seed sequence blast by mismatch #\n";
	print STDERR " Number of starting blast hits:\t", scalar keys %$blast_hits_r, "\n";


	foreach my $hit (keys %$blast_hits_r){		
		my $mismatch = ${$blast_hits_r->{$hit}}[4];
		delete $blast_hits_r->{$hit} unless $mismatch <= 1;	
		}

	my $end_hit_cnt = scalar keys %$blast_hits_r;
	
	print STDERR " Number of blast hits with <= mismatch:\t\t", scalar keys %$blast_hits_r, "\n";	

	}

sub screen_blast_by_pam{
	my ($blast_hits_r, $genome_r, $pam) = @_;
	
	print STDERR "# filtering by PAM #\n";
	print STDERR " Number of starting blast hits:\t", scalar keys %$blast_hits_r, "\n";
	
	foreach my $hit (keys %$blast_hits_r){		
		my $scaf = ${$blast_hits_r->{$hit}}[1];
		my $sstart = ${$blast_hits_r->{$hit}}[8];
		my $send = ${$blast_hits_r->{$hit}}[9];
		
		# scaf #
		die " ERROR: $scaf not found in genome!\n"
			unless exists $genome_r->{$scaf};
		
		# strand #
		my $strand = "+";
		if($sstart > $send){ # start > end
			$strand = "-";
			($sstart, $send) = ($send, $sstart);
			}
		
		# getting pam #
		my $ext_len = length $pam;
		my $pam_gen = substr($genome_r->{$scaf}, $sstart - $ext_len, $ext_len);
		$pam_gen = revcomp($pam_gen) if $strand eq "-";
		
		# checking pam #
		(my $pam_regex = $pam) =~ s/N/[ATGC]/g;
		$pam_regex = qr/$pam_regex/;
		delete $blast_hits_r->{$hit} unless $pam_gen =~ /$pam_regex/;		
		}

	my $end_hit_cnt = scalar keys %$blast_hits_r;
	
	print STDERR " Number of blast hits with PAM:\t\t", scalar keys %$blast_hits_r, "\n";
	}

sub revcomp{
	# reverse complements DNA #
	my $seq = shift;
	$seq = reverse($seq);
	$seq =~ tr/[a-z]/[A-Z]/;
	$seq =~ tr/ACGTNBVDHKMRYSW\.-/TGCANVBHDMKYRSW\.-/;
	return $seq;
	}

sub call_blastn{
	my ($blast_params, $genome_in) = @_;
	
	#my $outfmt = "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen'";
	
	my $cmd = "blastn -task 'blastn-short' -query spacer_trunc.fasta -subject $genome_in -outfmt 6 $blast_params |";
	
	open PIPE, $cmd or die $!;
	my %blast_hits;
	while(<PIPE>){
		chomp;
		
		my @line = split /\t/;
		$blast_hits{$.} = \@line;
		}
	close PIPE;
	
		#print Dumper %blast_hits; exit;
	return \%blast_hits;
	}

sub write_pam_spacer{
	my ($fasta_r, $fasta_in) = @_;
	
	$fasta_in =~ s/\.[^.]+$|$//;
	$fasta_in .= "_pam.fasta";
	
	open OUT, ">$fasta_in" or die $!;
	foreach my $name(keys %$fasta_r){
		print OUT "$name\n$fasta_r->{$name}\n";
		}
	close OUT;
		#print Dumper $fasta_in; exit;
	}

sub truncate_spacers{
	my ($fasta_r, $trunc) = @_;
	
	open OUT, ">spacer_trunc.fasta" or die $!;
	
	foreach my $name (keys %$fasta_r){
		#print "$name\n";
		my $trunc_seq = substr($fasta_r->{$name},0,$trunc);

		print OUT ">$name\n$trunc_seq\n";
		}
	close OUT;
	
	}

sub load_fasta{
	my ($fasta_in) = @_;
	
	my $cmd = "perl -pi -e 's/\r/\n/g' $fasta_in";
	`$cmd`;
	
	open IN, $fasta_in or die $!;
	my (%fasta, $tmpkey);
	while(<IN>){
		chomp;
 		s/#.+//;
 		next if /^\s*$/;	
 		
 		if(/>.+/){
 			$_ =~ s/^>//;
 			$fasta{$_} = "";
 			$tmpkey = $_;	# changing key
 			}
 		else{ $fasta{$tmpkey} .= $_; }
		}
	close IN;
	
		#print Dumper %fasta;
	return \%fasta;
	}
	


__END__

=pod

=head1 NAME

spacerBlastSeed.pl -- use blastn to find spacer seed sequences

=head1 SYNOPSIS

spacerBlastSeed.pl [flags] > blast_filtered.fasta

=head2 required flags

=over

=item -spacer  <char>

Fasta of spacer sequences. 

=item -genome  <char>

Fasta of a genome to query.

=back 

=head2 Optional flags

=over

=item -truncate  <int>

Trancate spacer to just SEED for blasting (bp from 5' end). [8]

=item -blast  <char>

Blast parameters besides -outfmt. ['-evalue 100']

=item -pam  <char>

PAM sequence to filter SEED sequence blast hits. [CCN]

=item -verbose  <bool>

Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc spacerBlastSeed.pl

=head1 DESCRIPTION

Workflow:

=over

=item 1) Truncate spacer to seed sequence

=item 2) Blast seed sequence (blastn-short)

=item 3) Filter seed blast by mismatches (<= 1)

=item 4) Filter seed blast by adjacent PAM (must match '-pam')

=item 5) Filter seed blast by mismatches to rest of protospacer (<= 9)

=back


=head1 EXAMPLES

=head2 Basic usage:

spacerBlastSeed.pl -s spacers.fasta -g genome.fasta > blast_filtered.fasta

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

