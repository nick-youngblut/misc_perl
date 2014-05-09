#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::TreeIO;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $species_in);
GetOptions(
	   "species=s" => \$species_in,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );
	   
### I/O error & defaults
map{ die " ERROR: $_ not found!\n" unless -e $_ } @ARGV;

### MAIN
my $lca_r = load_species_tree($species_in);
foreach my $parse_in (@ARGV){
	parse_notung_reconcile($parse_in, $lca_r);
	}

### Subroutines
sub parse_notung_reconcile{
# parsing notung reconcile #
	my ($parse_in, $lca_r) = @_;

	open my $in_fh, $parse_in or die $!;
	while(<$in_fh>){
		chomp;
		next if /^\s*$/;
		
		if(/^#D/){
			parse_duplication_block($parse_in, $in_fh, $lca_r);
			}
		elsif(/^#CD/){
			parse_codivergence_block($parse_in, $in_fh, $lca_r);
			}
		elsif(/^#T/){
			parse_transfer_block($parse_in, $in_fh, $lca_r);
			}
		elsif(/^#L/){
			parse_loss_block($parse_in, $in_fh, $lca_r);
			}
		}
	close $in_fh or die $!;
	
	}

sub parse_loss_block{
	my ($parse_in, $in_fh, $lca_r) = @_;
	
	while(<$in_fh>){
		chomp;
		last if /^\s*$/;
		
		my @line = split /\t/;
		next unless $line[2]; 		# if no losses
		
		$line[1] = $lca_r->{$line[1]} 
			if exists $lca_r->{$line[1]};
	
		print join("\t", $parse_in, "NA", $line[1],
				"Loss", "NA"), "\n";
		}
	}

sub parse_transfer_block{
	my ($parse_in, $in_fh, $lca_r) = @_;
	
	while(<$in_fh>){
		chomp;
		last if /^\s*$/;
		
		my @line = split /\t/;
		$line[3] = $lca_r->{$line[3]} 
			if exists $lca_r->{$line[3]};
		$line[4] = $lca_r->{$line[4]} 
			if exists $lca_r->{$line[4]};

		print join("\t", $parse_in, $line[1], $line[3],
				"Transfer", $line[4]), "\n";
		}
	}

sub parse_codivergence_block{
# parsing codivergences #
	my ($parse_in, $in_fh, $lca_r) = @_;
	
	my $last_blank = 1; 		
	while(<$in_fh>){
		chomp;

		# end of block = 2 blank lines in a row #
		if(/^\s*$/){
			last if $last_blank;
			$last_blank = 1;
			}
		else{ $last_blank = 0; }
		next if /^\s*$/;
		
		# writing out line #
		my @line = split /\t/;
		$line[1] = $lca_r->{$line[1]} 
			if exists $lca_r->{$line[1]};
		print join("\t", $parse_in, "NA", $line[1],
				"Co-Divergence", "NA"), "\n";
		}
	}

sub parse_duplication_block{
# parsing duplications #
	my ($parse_in, $in_fh, $lca_r) = @_;
	
	my $last_blank = 1; 		
	while(<$in_fh>){
		chomp;

		# end of block = 2 blank lines in a row #
		if(/^\s*$/){
			last if $last_blank;
			$last_blank = 1;
			}
		else{ $last_blank = 0; }
		next if /^\s*$/;
		
		# writing out line #
		my @line = split /\t/;
		$line[1] = $lca_r->{$line[1]} 
			if exists $lca_r->{$line[1]};		
		print join("\t", $parse_in, "NA", $line[1],
				"Duplication", "NA"), "\n";
		}
	}


sub load_species_tree{
# loading species tree; getting node LCA for each node ID #
	my ($species_in) = @_;
	my $treeio = Bio::TreeIO->new( -file => $species_in, 
					-format => "newick");
	
	# getting LCA #
	my %lca;
	my $node_cnt;
	while (my $tree = $treeio->next_tree){
		for my $node ($tree->get_nodes){
			next if $node->is_Leaf;
			my $node_id = $node->id;
			$node_cnt++;

			# getting all leaves of node #
			my @children;
			for my $child ($node->get_all_Descendents){
				push @children, $child if $child->is_Leaf;
				}
			
			# finding a pair of leaves where LCA is node #
			for (my $i=0; $i<=$#children; $i++){
				last if exists $lca{$node_id};
				for (my $ii=$#children; $i>=0; $i--){			# working backward
					last if exists $lca{$node_id};
					
					my $lca_id = $tree->get_lca( @children[($i, $ii)] )->id;
					if($lca_id =~ /[A-Za-z]/ && $lca_id eq $node_id){
						$lca{$node_id} = join("|", $children[$i]->id, $children[$ii]->id);
						}
					elsif( $lca_id == $node_id){
						$lca{$node_id} = join("|", $children[$i]->id, $children[$ii]->id);
						}
					}
				}
			
			$lca{$node_id} = "NA" unless exists $lca{$node_id};
			}
		}
		
		#print Dumper %lca; exit;
	return \%lca;
	}



__END__

=pod

=head1 NAME

Notung_reconcile_parse.pl -- parse Notung '--reconcile' output into a *txt table

=head1 SYNOPSIS

Notung_reconcile_parse.pl [flags] *parsable.txt > reconcile.txt

=head2 Required flags

=over

=item -s

Species tree produced by Notung '-stpruned --treeoutput newick'

=back

=head2 Optional flags

=over

=item -h	This help message

=back

=head2 For more information:

perldoc Notung_reconcile_parse.pl

=head1 DESCRIPTION

The output is a tab-delimited table. Blank value = not applicable.

=head2 Columns of output

=over

=item * 	Tree ID (file name)

=item * 	Gene tree node name

=item * 	Species tree node name

=item * 	Category (Duplication, Transfer, Speciation, or Leaf)

=item * 	Recipient (species tree node name)

=back

=head1 EXAMPLES

=head2 Usage: 

Notung_reconcile_parse.pl -s species.nwk genetree1.parsable.txt > reconcile.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

