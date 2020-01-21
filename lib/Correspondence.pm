
package Correspondence;
use Data::Dumper;
use CorrTableEntry;

use strict;
use SeedUtils;

my %map_cache;

sub fill_cache
{
    my($sap, $genomes) = @_;

    for my $i (0..$#$genomes)
    {
	my $genome1 = $genomes->[$i];
	for my $j ($i+1 .. $#$genomes)
	{
	    my $genome2 = $genomes->[$j];
	    my $corr = $sap->gene_correspondence_map(-genome1 => $genome1,
						     -genome2 => $genome2,
						     -fullOutput => 1);
	    warn "Loaded $genome1 $genome2\n";
	    $map_cache{$genome1, $genome2} = $corr;
	    $map_cache{$genome2, $genome1} = $corr;
	}
    }	
}

sub get_correspondence_entry
{
    my($sap, $peg, $genome2) = @_;

    my $genome1 = SeedUtils::genome_of($peg);

    my $corr = $map_cache{$genome1, $genome2};

    if (!defined($corr))
    {
	$corr = $sap->gene_correspondence_map(-genome1 => $genome1,
					      -genome2 => $genome2,
					      -fullOutput => 1);
	$map_cache{$genome1, $genome2} = $corr;
    }

    my($tuple) = grep { $_->[0] eq $peg } @$corr;

    if (!defined($tuple))
    {
	warn "No entry for $peg\n";
	return undef;
    }
    bless $tuple, 'CorrTableEntry';
    return $tuple;
}

sub is_connected_to
{
    my($sap, $peg, $genome2, %args) = @_;

    # args are:
    # $bbh_required, $context_size, $num_matching_functions_in_context, $max_psc, $min_iden, $min_coverage) = @_;

    my $genome1 = SeedUtils::genome_of($peg);

    my $tuple = get_correspondence_entry($sap, $peg, $genome2);
    if (!defined($tuple))
    {
	return undef;
    }

    if (entry_meets_criteria($sap, $tuple, %args))
    {
	return $tuple->id2;
    }
}

sub entry_meets_criteria
{
    my($sap, $tuple, %args) = @_;

    my $peg = $tuple->id1;

    if ($args{-bbhRequired} && ($tuple->hitinfo ne '<=>'))
    {
	warn "Reject $peg to ", $tuple->id2, " due to missing BBH\n";
	return undef;
    }

    if (defined($args{-contextSize}) && ($tuple->npairs < $args{-contextSize}))
    {
	warn "Reject $peg to ", $tuple->id2, " due to short context\n";
	return undef;
    }

    if (defined($args{-maxPsc}) && ($tuple->psc > $args{-maxPsc}))
    {
	warn "Reject $peg to ", $tuple->id2, " due to psc\n";
	return undef;
    }

    if (defined($args{-minIden}) && ($tuple->iden < $args{-minIden}))
    {
	warn "Reject $peg to ", $tuple->id2, " due to iden\n";
	return undef;
    }

    if (defined(my $min_cover1 = $args{-minCoverage1}))
    {
	my $coverage = ($tuple->end1 - $tuple->beg1 + 1) / $tuple->ln1;
	if ($coverage < $min_cover1)
	{
	    warn "Reject $peg to ", $tuple->id2, " due to coverage1\n";
	    return undef;
	}
    }

    if (defined(my $min_cover2 = $args{-minCoverage2}))
    {
	my $coverage = ($tuple->end2 - $tuple->beg2 + 1) / $tuple->ln2;
	if ($coverage < $min_cover2)
	{
	    warn "Reject $peg to ", $tuple->id2, " due to coverage2\n";
	    return undef;
	}
    }

    if (defined(my $min_matching = $args{-numMatchingFunctionsInContext}))
    {
	my $count = $tuple->num_matching_functions;
	if (!defined($count))
	{
	    if (!defined($sap))
	    {
		die "-numMatchingFunctionsInContext requires either a SAP object or a 19-column correspondence entry";
	    }
	    
	    my @pairs = $tuple->pairs;
	    my @pegs = map { @$_ } @pairs;
	    my $fns = $sap->ids_to_functions(-ids => \@pegs);
	    
	    $count = 0;
	    for my $pair (@pairs)
	    {
		my $fn1 = $fns->{$pair->[0]};
		my $fn2 = $fns->{$pair->[1]};
		$count++ if ($fn1) && ($fn1 eq $fn2);
	    }
	}
	    
	if ($count < $min_matching)
	{
	    warn "Reject $peg to ", $tuple->id2, " due to min matching functions\n";
	    return undef;
	}
    }

    return 1;
}

1;
