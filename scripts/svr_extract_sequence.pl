#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use SAPserver;
use ScriptThing;
use SeedUtils;
use gjoseqlib;

#
#	This is a SAS Component.
#

=head1 svr_extract_sequence

    svr_extract_sequence input-fasta < locations.tbl > fasta

Extract substrings of sequence data from an input fasta file.

This script takes as input a tab-delimited file of locations, one per line.
These locations may either be location strings in either of the SEED standard
forms, or they may be tab-delimited lines containing contig id, begin, and end.

For each location, we will extract the given substring of sequence 
data from the fasta file provided as a command line argument.

=head2 Command-Line Options

=back

=cut

# Parse the command-line options.

if (@ARGV != 1) {
    print "usage: svr_extract_sequence fasta-file < locations.tbl > fasta-file\n";
    exit;
}

my $fasta = shift;

if (!open(FASTA, "<", $fasta))
{
    die "Cannot open fasta file $fasta: $!";
}

#
# Scan the input and build the set of contig/begin/end sets we are looking for.
#

my %wanted;

while (<STDIN>)
{
    chomp;
    my ($id, @list) = split(/\t/);
    my($contig, $begin, $end);
    if (@list == 1)
    {
	my $loc = $list[0];
	my $strand;
	($contig, $begin, $end, $strand) = &SeedUtils::parse_location($loc);
	if ($strand eq '-')
	{
	    ($begin, $end) = ($end, $begin);
	}
    }
    elsif (@list == 3)
    {
	($contig, $begin, $end) = @list;
    }
    else
    {
	warn "Invalid input at line $.\n";
	next;
    }
    push(@{$wanted{$contig}}, [$id, $begin, $end]);
}

while (my($id, $def, $seq) = gjoseqlib::read_next_fasta_seq(\*FASTA))
{
    for my $wanted (@{$wanted{$id}})
    {
	my($str_id, $beg, $end) = @$wanted;
	my $subseq = gjoseqlib::DNA_subseq($seq, $beg, $end);
	gjoseqlib::print_alignment_as_fasta([$str_id, undef, $subseq]);
    }
}
