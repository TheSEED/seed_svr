#
# Copyright (c) 2003-2014 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
#
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License.
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#

package AlignsAndTreesServer;
#===============================================================================
#  perl functions for loading and accessing Alignments and Trees
#
#  Usage:  use AlignsAndTreesServer;
#
#-------------------------------------------------------------------------------
#  Alignments
#-------------------------------------------------------------------------------
#
#    @alignIDs = all_alignIDs( [$SAPserverO] );
#   \@alignIDs = all_alignIDs( [$SAPserverO] );
#
#    @alignIDs = aligns_with_md5ID( [$SAPserverO,] $md5 );
#   \@alignIDs = aligns_with_md5ID( [$SAPserverO,] $md5 );
#
#    @md5IDs   = md5IDs_in_align( [$SAPserverO,] $alignID );
#   \@md5IDs   = md5IDs_in_align( [$SAPserverO,] $alignID );
#
#   \%aligns_of_md5ID = md5IDs_to_aligns( [$ALITREserverO,] @md5IDs );
#
#        $aligns_of_md5ID{ $md5ID } = \@alignIDs
#
#   \@seqs                       = md5_alignment_by_ID( [ALITREserverO,] $alignID );
# ( \@seqs, \%md5_row_metadata ) = md5_alignment_by_ID( [ALITREserverO,] $alignID );
#           \%md5_row_metadata   = md5_alignment_metadata( [ALITREserverO,] $alignID );
#
#        $md5_row_metadata{ $seqID } = [ $md5ID, $peg_length, $trim_beg, $trim_end, $location_string ]
#
#   \%md5_row_metadata = alignment_metadata_by_md5( [$ALITREserverO,] $alignID, \@md5IDs );
#
#        $md5_row_metadata{ $seqID } = [ $md5ID, $peg_length, $trim_beg, $trim_end, $location_string ]
#
#    @metadata = alignments_metadata_by_md5( [$SAPserverO,] \@md5IDs );
#   \@metadata = alignments_metadata_by_md5( [$SAPserverO,] \@md5IDs );
#
#        @metadata = ( [ $alignID, $seqID, $md5, $peg_length, $trim_beg, $trim_end, $location_string ], ... )
#
#    @alignID_coverages = alignment_coverages_of_md5( [$ALITREserverO,] $md5 );
#   \@alignID_coverages = alignment_coverages_of_md5( [$ALITREserverO,] $md5 );
#    @alignID_coverages = alignment_coverages_of_md5( [$ALITREserverO,] $md5, @alignIDs );
#   \@alignID_coverages = alignment_coverages_of_md5( [$ALITREserverO,] $md5, @alignIDs );
#
#        Return values are:  [ $alignID, $md5, $len, \@coverages ]
#        A coverage is:      [ $beg, $end, $loc ]   
#
#    @alignIDs = aligns_with_pegID( [$SAPserverO,] $pegID );
#   \@alignIDs = aligns_with_pegID( [$SAPserverO,] $pegID );
#
#   \%aligns_of_pegID = pegIDs_to_aligns( [@servers,] @pegIDs );
#
#        $aligns_of_pegID{ $pegID } = \@alignIDs
#
#    @pegIDs = pegIDs_in_align( [$SAPserverO,] $alignID );
#   \@pegIDs = pegIDs_in_align( [$SAPserverO,] $alignID );
#
#   \@seqs                   = peg_alignment_by_ID( [ALITREserverO,] $alignID );
# ( \@seqs, \%peg_metadata ) = peg_alignment_by_ID( [ALITREserverO,] $alignID );
#
#   \%peg_metadata = peg_alignment_metadata( [SAPserverO,] $alignID );
#
#        $peg_metadata{ $seqID } = [ $pegID, $peg_length, $trim_beg, $trim_end, $location_string ]
#
#  seqIDs are based on the md5ID or pegID, appending a hyphen and positive
#  integer when necessary to disambiguate multiple sequences in an alignment
#  or tree that are derived from the same underlying protein.
#
#-------------------------------------------------------------------------------
#  Trees
#-------------------------------------------------------------------------------
#
#    @treeIDs = all_treeIDs( [$SAPserverO] );
#   \@treeIDs = all_treeIDs( [$SAPserverO] );
#
#    @treeIDs = trees_with_md5ID( [$SAPserverO,] $md5 );
#   \@treeIDs = trees_with_md5ID( [$SAPserverO,] $md5 );
#
#    @md5IDs = md5IDs_in_tree( [$SAPserverO,] $treeID );
#   \@md5IDs = md5IDs_in_tree( [$SAPserverO,] $treeID );
#
#   \%trees_of_md5ID = md5IDs_to_trees( [$ALITREserverO,] @md5IDs );
#
#   \%trees_of_pegID = pegIDs_to_trees( [@servers,] @pegIDs );
#
#    $tree                   = md5_tree_by_ID( $treeID );
#  ( $tree, \%md5_metadata ) = md5_tree_by_ID( $treeID );
#
#    @treeIDs = trees_with_pegID( [$SAPserverO,] $pegID );
#   \@treeIDs = trees_with_pegID( [$SAPserverO,] $pegID );
#
#    @pegIDs = pegIDs_in_tree( [$SAPserverO,] $treeID );
#   \@pegIDs = pegIDs_in_tree( [$SAPserverO,] $treeID );
#
#    $tree                   = peg_tree_by_ID( [ALITREserverO,] $treeID );
#  ( $tree, \%peg_metadata ) = peg_tree_by_ID( [ALITREserverO,] $treeID );
#
#===============================================================================
#  Support for counting distinct roles in alignments and trees:
#
#    @role_count_pairs = roles_in_align( [$SAPserverO,] $alignID );
#   \@role_count_pairs = roles_in_align( [$SAPserverO,] $alignID );
#
#    $role  = majority_role_in_align( [$SAPserverO,] $alignID );
#
#    @role_count_pairs = roles_in_tree( [$SAPserverO,] $treeID );
#   \@role_count_pairs = roles_in_tree( [$SAPserverO,] $treeID );
#
#    $role  = majority_role_in_tree( [$SAPserverO,] $treeID );
#
#  Support for function function lookup for md5 IDs:
#
#    \%md5_function = md5s_to_functions(              @md5s );
#    \%md5_function = md5s_to_functions( $SAPserverO, @md5s );
#
#  When an md5ID maps to multiple pegIDs, only one pegID's function is counted.
#  
#-------------------------------------------------------------------------------
#  Support for getting projections between MD5 IDs:
#
#   \%md5_projections = get_md5_projections( [$SAPserverO,] @md5s [, \%opts] );
#
#-------------------------------------------------------------------------------
#  Get data on alignments with a protein
#
#   @align_data = data_on_aligns_with_prot( [$ALITREserverO,] [$SAPserverO,] $fid [, \%opts] );
#   @align_data = data_on_aligns_with_prot( [$ALITREserverO,] [$SAPserverO,] $md5 [, \%opts] );
#
#   or with $opts->{fid} or $opts->{md5}
#
#   @align_data = data_on_aligns_with_prot( [$ALITREserverO,] [$SAPserverO,]      [  \%opts] );
#
#   Alignment data:
#
#   [ $alignID, $coverage, \@roles ]
#
#        $coverage = [ $beg, $end, $len ]
#       \@roles    = [ [ $role1, $cnt1 ], [ $role2, $cnt2 ], ... ]
#
#   Options:
#
#       coverage => $frac  # minimum fraction of protein residues included
#       fid      => $fid   # feature id, if not an argument
#       md5      => $md5   # protein sequence md5, if not an argument
#       noroles  => $bool  # omit roles and counts (faster)
#       roles    => $nrole # limit on number of roles returned
#       sap      => $sapO  # Sapling server object
#
#-------------------------------------------------------------------------------
#  Support for md5 <-> pegID interconversion:
#
#  Each function supports an optional first argument that defines a SAPserver.
#  If it is supplied, it becomes the default for subsequence calls.
#
#    @pegIDs         = md5_to_pegs(               $md5 );  # One md5
#   \@pegIDs         = md5_to_pegs(               $md5 );
#    @pegIDs         = md5_to_pegs(  $SAPserverO, $md5 );
#   \@pegIDs         = md5_to_pegs(  $SAPserverO, $md5 );
#
#    %pegIDs_of_md5  = md5s_to_pegs(              @md5s );  # One or more md5s
#   \%pegIDs_of_md5  = md5s_to_pegs(              @md5s );
#    %pegIDs_of_md5  = md5s_to_pegs( $SAPserverO, @md5s );
#   \%pegIDs_of_md5  = md5s_to_pegs( $SAPserverO, @md5s );
#
#    $md5            = peg_to_md5(               $pegID );  # One pegID
#    $md5            = peg_to_md5(  $SAPserverO, $pegID );
#
#   \%md5s_of_pegIDs = pegs_to_md5(              @pegIDs );  # One or more pegIDs
#   \%md5s_of_pegIDs = pegs_to_md5( $SAPserverO, @pegIDs );
#
#-------------------------------------------------------------------------------

use strict;
use Data::Dumper;
use SeedUtils;
use SAPserver;
use ALITREserver;
use gjonewicklib;
#   BlastInterface is invoked with require, for the few functions that use it.

my $SAPserverO;
my $ALITREserverO;

our $ReferenceFIG;

my %peg_to_md5;
my %md5_to_pegs;

sub md5_to_pegs
{
    $SAPserverO = shift if eval { $_[0]->isa( 'SAPserver' ) };

    my ( $md5 ) = @_;
    return wantarray ? () : [] if ! $md5;

    my @pegIDs;
    if ( exists $md5_to_pegs{ $md5 } )
    {
        @pegIDs = @{ $md5_to_pegs{ $md5 } };
    }
    else
    {
	if ($ReferenceFIG)
	{
	    @pegIDs = $ReferenceFIG->pegs_with_md5($md5);
	}
	else
	{
	    $SAPserverO ||= SAPserver->new() or die "Could not get a new SAPserver\n";
	    my $pegIDsH = $SAPserverO->proteins_to_fids( -prots => [ $md5 ] );
	    @pegIDs = $pegIDsH ? @{ $pegIDsH->{ $md5 } } : ( );
	}
	$md5_to_pegs{ $md5 } = [ @pegIDs ];
    }

    wantarray ? @pegIDs : \@pegIDs;
}


sub md5s_to_pegs
{
    $SAPserverO = shift if eval { $_[0]->isa( 'SAPserver' ) };

    my @md5 = @_;
    return wantarray ? () : {} if ! @md5;

    #  Remove the md5 ids already in the cache

    @md5 = grep { ! exists $md5_to_pegs{ $_ } } @md5;

    #  Get the remaining md5 ids from the Sapling:

    if ( @md5 )
    {
	if ($ReferenceFIG)
	{
	    for my $md5 (@md5)
	    {
		my @pegs = $ReferenceFIG->pegs_with_md5($md5);
		push(@{$md5_to_pegs{$md5}}, $_) foreach @pegs;
	    }
	}
	else
	{
	    $SAPserverO ||= SAPserver->new() or die "Could not get a new SAPserver\n";
	    my $pegIDsH = $SAPserverO->proteins_to_fids( -prots => \@md5 ) || {};
	    foreach ( @md5 ) { $md5_to_pegs{ $_ } = $pegIDsH->{ $_ } }
	}
    }

    #  Return the whole cache

    wantarray ? %md5_to_pegs : \%md5_to_pegs;
}


sub peg_to_md5
{
    $SAPserverO = shift if eval { $_[0]->isa( 'SAPserver' ) };

    my ( $pegID ) = @_;
    return undef if ! $pegID;

    return $peg_to_md5{ $pegID } if exists $peg_to_md5{ $pegID };

    my $md5;
    if ($ReferenceFIG)
    {
	$md5 = $ReferenceFIG->md5_of_peg($pegID);
    }
    else
    {
	$SAPserverO ||= SAPserver->new() or die "Could not get a new SAPserver\n";
	my $md5H = $SAPserverO->fids_to_proteins( -ids => [ $pegID ] );
	$md5  = $md5H ? $md5H->{ $pegID } : undef;
    }

    $peg_to_md5{ $pegID } = $md5;
}


sub pegs_to_md5
{
    $SAPserverO = shift if eval { $_[0]->isa( 'SAPserver' ) };

    my @pegID = @_;
    return wantarray ? () : {} if ! @pegID;

    #  Remove the pegIDs already in the cache

    @pegID = grep { ! exists $peg_to_md5{ $_ } } @pegID;

    #  Get the remaining pegIDs from the Sapling:

    if ( @pegID )
    {
	if ($ReferenceFIG)
	{
	    my $ret = $ReferenceFIG->md5_of_peg_bulk(\@pegID);
	    for my $peg (@pegID)
	    {
		$peg_to_md5{$peg} = $ret->{$peg};
	    }
	}
	else
	{
	    $SAPserverO ||= SAPserver->new() or die "Could not get a new SAPserver\n";
	    my $md5H = $SAPserverO->fids_to_proteins( -ids => \@pegID ) || {};
	    foreach ( @pegID ) { $peg_to_md5{ $_ } = $md5H->{ $_ } }
	}
    }

    #  Return the whole cache

    wantarray ? %peg_to_md5 : \%peg_to_md5;
}


#===============================================================================
#  Alignments
#
#  Each function supports an optional first argument that defines a ALITREserver.
#  If it is supplied, it becomes the default for subsequence calls.
#===============================================================================
#
#    @alignIDs = all_alignIDs( [$ALITREserverO] );
#   \@alignIDs = all_alignIDs( [$ALITREserverO] );
#
#-------------------------------------------------------------------------------
sub all_alignIDs
{
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new()
        or die "Could not get a new ALITREserver\n";

    my $alignIDs = $ALITREserverO->all_alignIDs() || [];

    wantarray ? @$alignIDs : $alignIDs;
}

#-------------------------------------------------------------------------------
#  md5 based alignments:
#-------------------------------------------------------------------------------
#
#    @alignIDs = aligns_with_md5ID( [$ALITREserverO,] $md5 );
#   \@alignIDs = aligns_with_md5ID( [$ALITREserverO,] $md5 );
#
#-------------------------------------------------------------------------------
sub aligns_with_md5ID
{
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new()
        or die "Could not get a new ALITREserver\n";

    my $md5      = $_[0] || '';
    my $alignH   = $ALITREserverO->aligns_with_md5ID( -ids => [ $md5 ] ) || {};
    my $alignIDs = $alignH->{ $md5 } || [];

    wantarray ? @$alignIDs : $alignIDs;
}


#-------------------------------------------------------------------------------
#  Run through a list of MD5 IDs, and produce a hash table that maps
#  each MD5 ID to a list of alignments that contain it.
#
#   \%aligns_of_md5ID = md5IDs_to_aligns( [$ALITREserverO,] @md5IDs );
#
#-------------------------------------------------------------------------------
sub md5IDs_to_aligns
{
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new()
        or die "Could not get a new ALITREserver\n";

    my $aligns_of_md5ID = $ALITREserverO->aligns_with_md5ID( -ids => \@_ ) || {};

    $aligns_of_md5ID;
}


#-------------------------------------------------------------------------------
#
#    @md5IDs = md5IDs_in_align( [$ALITREserverO,] $alignID );
#   \@md5IDs = md5IDs_in_align( [$ALITREserverO,] $alignID );
#
#-------------------------------------------------------------------------------
my %atID_to_md5s;
sub md5IDs_in_align
{
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };

    ( my $alignID = $_[0] || '' ) =~ s/^aln//;
    my $md5IDs;

    if ( exists $atID_to_md5s{ $alignID } )
    {
        $md5IDs = $atID_to_md5s{ $alignID };
    }
    else
    {
        $ALITREserverO ||= ALITREserver->new()
            or die "Could not get a new ALITREserver\n";

        my $md5H = $ALITREserverO->md5IDs_in_align( -ids => [ $alignID ] ) || {};
        $md5IDs  = $md5H->{ $alignID } || [];

        $atID_to_md5s{ $alignID } = $md5IDs;
    }

    wantarray ? @$md5IDs : $md5IDs;
}


#-------------------------------------------------------------------------------
#
#    %alignID_to_md5s = md5IDs_in_aligns( [$ALITREserverO,] @alignIDs );
#   \%alignID_to_md5s = md5IDs_in_aligns( [$ALITREserverO,] @alignIDs );
#
#-------------------------------------------------------------------------------
sub md5IDs_in_aligns
{
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };

    s/^aln// for ( my @alignIDs = @_ );
    my @newIDs = grep { ! exists $atID_to_md5s{ $_ } } @alignIDs;
    
    if ( @newIDs )
    {
        $ALITREserverO ||= ALITREserver->new()
            or die "Could not get a new ALITREserver\n";
        
        my $md5H = $ALITREserverO->md5IDs_in_align( -ids => \@newIDs ) || {};
        foreach ( @newIDs)
        {
            my $md5s = $md5H->{ $_ };
            $atID_to_md5s{ $_ } = defined $md5s ? $md5s : [];
        }
    }

    wantarray ? map { $_ => $atID_to_md5s{ $_ } } @alignIDs : \%atID_to_md5s;
}


#-------------------------------------------------------------------------------
#
#    @unique_md5IDs = unique_md5IDs_in_aligns( [$ALITREserverO,] @alignIDs );
#   \@unique_md5IDs = unique_md5IDs_in_aligns( [$ALITREserverO,] @alignIDs );
#
#-------------------------------------------------------------------------------
sub unique_md5IDs_in_aligns
{
    $ALITREserverO   = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";

    s/^aln// for ( my @alignIDs = @_ );
    my $id_to_md5s = md5IDs_in_aligns( $ALITREserverO, @alignIDs );
    my %seen;
    my @unique_md5IDs = grep { ! $seen{ $_ }++ } map { @{ $id_to_md5s->{ $_ } } } @_;

    wantarray ? @unique_md5IDs : \@unique_md5IDs;
}


#-------------------------------------------------------------------------------
#  Alignment in md5IDs
#
#   \@seqs               = md5_alignment_by_ID( [$ALITREserverO,] $alignID );
# ( \@seqs, \%metadata ) = md5_alignment_by_ID( [$ALITREserverO,] $alignID );
#
#      $metadata{ $seqID } = [ $md5ID, $peg_length, $trim_beg, $trim_end, $location_string ]
#
#-------------------------------------------------------------------------------
sub md5_alignment_by_ID
{
    $ALITREserverO   = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";

    ( my $alignID = $_[0] || '' ) =~ s/^aln//;
    #***************************************************************************
    #  This function returns both the alignment and metadata for each alignID.
    #***************************************************************************
    my $alignH = $ALITREserverO->md5_alignment_by_ID( -ids => [ $alignID ] );
    return wantarray ? () : [] unless $alignH && $alignH->{ $alignID };

    my ( $align, $meta ) = @{ $alignH->{ $alignID } };
    my %indexed_seq = map { $_->[0] => $_ } @$align;
    my @align = map { $indexed_seq{ $_ } } md5_tree_order( $ALITREserverO, $alignID );
    wantarray ? ( \@align, $meta ) : \@align;
}


sub md5_tree_order
{
    my $md5_tree = md5_tree_by_ID( @_ );
    my @tips = $md5_tree ? gjonewicklib::newick_tip_list( gjonewicklib::aesthetic_newick_tree( gjonewicklib::reroot_newick_to_midpoint_w( $md5_tree ) ) ) : ();
    wantarray ? @tips : \@tips;
}


#-------------------------------------------------------------------------------
#  Alignment metadata in md5IDs
#
#   \%md5_row_metadata = md5_alignment_metadata( [ALITREserverO,] $alignID );
#
#        $md5_row_metadata{ $seqID } = [ $md5ID, $peg_length, $trim_beg, $trim_end, $location_string ]
#
#   where:
#
#       $ALITREserverO is an alignment and tree server object.
#       $alignID is an alignment whose MD5 relationship data is desired.
#
#-------------------------------------------------------------------------------
sub md5_alignment_metadata
{
    $ALITREserverO   = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";

    ( my $alignID = $_[0] || '' ) =~ s/^aln//;

    my $metaH = $ALITREserverO->md5_alignment_metadata( -ids => [ $alignID ] ) || {};

    $metaH->{ $alignID } || {};
}


#-------------------------------------------------------------------------------
#  Get row metadata for one or more specific md5s in a specific alignment:
#
#    %metadata = alignment_metadata_by_md5( [$ALITREserverO,] $alignID, \@md5IDs );
#   \%metadata = alignment_metadata_by_md5( [$ALITREserverO,] $alignID, \@md5IDs );
#
#   where:
#
#       $ALITREserverO is an alignment and tree server object.
#       $alignID is an alignment whose MD5 relationship data is desired.
#      \@md5IDs is a reference to a list of the md5IDs for which metadata are desired.
#
#       $metadata{ $seqID } = [ $md5, $peg_length, $trim_beg, $trim_end, $location_string ]
#-------------------------------------------------------------------------------
sub alignment_metadata_by_md5
{
    $ALITREserverO   = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";

    my ( $alignID, $md5IDs ) = @_;
    my $metadata;

    if ( $alignID && $md5IDs && ref($md5IDs) eq 'ARRAY' && @$md5IDs )
    {
        $metadata = $ALITREserverO->alignment_metadata_by_md5( -ids => [ $alignID, @$md5IDs ] );
    }
    
    $metadata;
}


#-------------------------------------------------------------------------------
#
#    @metadata = alignments_metadata_by_md5( [$ALITREserverO,] \@md5IDs );
#   \@metadata = alignments_metadata_by_md5( [$ALITREserverO,] \@md5IDs );
#
#   where:
#
#       $ALITREserverO is an alignment and tree server object.
#      \@md5IDs is a list of the md5IDs for which the data are desired.
#
#       @metadata = ( [ $alignID, $seqID, $md5, $peg_length, $trim_beg, $trim_end, $loc_string ], ... )
#-------------------------------------------------------------------------------
sub alignments_metadata_by_md5
{
    $ALITREserverO   = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";

    my ( $md5IDs ) = @_;
    my $metadata = [];

    if ( $md5IDs && ref($md5IDs) eq 'ARRAY' && @$md5IDs )
    {
        $metadata = $ALITREserverO->alignments_metadata_by_md5( -ids => $md5IDs ) || [];
    }

    wantarray ? @$metadata : $metadata;
}


#-------------------------------------------------------------------------------
#  Coverage of an md5 sequence by one or more alignments. If the alignIDs are
#  not supplied, they are retrieved from the server. An md5 can occur more
#  than once in an alignment, so the value for each alignment is a list of
#  coverages.
#
#    @alignID_coverages = alignment_coverages_of_md5( [$ALITREserverO,] $md5 );
#   \@alignID_coverages = alignment_coverages_of_md5( [$ALITREserverO,] $md5 );
#    @alignID_coverages = alignment_coverages_of_md5( [$ALITREserverO,] $md5, @alignIDs );
#   \@alignID_coverages = alignment_coverages_of_md5( [$ALITREserverO,] $md5, @alignIDs );
#
#  Return values are:  [ $alignID, $md5, $len, \@coverages ]
#  A coverage is:      [ $beg, $end, $loc ]   
#
#  Only alignments with a region covered by the md5 are returned.
#-------------------------------------------------------------------------------
sub alignment_coverages_of_md5
{
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new()
        or die "Could not get a new ALITREserver\n";

    my $md5 = shift or return wantarray ? () : [];
    my %keep = map { $_ => 1 } @_;

    # [ $alignID, $seqID, $md5, $peg_length, $trim_beg, $trim_end, $location_string ]
    my %metadata;
    foreach ( @{ alignments_metadata_by_md5( [$md5] ) || [] } )
    {
        push @{ $metadata{ $_->[0] } }, $_;
    }

    my @alignID_coverages;
    foreach my $alignID ( sort keys %metadata )
    {
        next if %keep && ! $keep{ $alignID };

        my @rows = @{ $metadata{$alignID} };
        if ( @rows )
        {
            my $len  = $rows[0]->[3];
            push @alignID_coverages, [ $alignID,
                                       $md5,
                                       $len,
                                       [ map { [@$_[4..6]] } @rows ]
                                     ];
        }
    }

    wantarray ? @alignID_coverages : \@alignID_coverages;
}


#-------------------------------------------------------------------------------
#  peg based alignments:
#-------------------------------------------------------------------------------
#  Alignments with a given pegID
#
#    @alignIDs = aligns_with_pegID( [@servers,] $pegID );
#   \@alignIDs = aligns_with_pegID( [@servers,] $pegID );
#
#    @servers can be $ALITREserverO and/or $SAPserverO.
#-------------------------------------------------------------------------------
sub aligns_with_pegID
{
    #  I will take neither, either or both objects in any order.
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $SAPserverO    = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    #  And then I will supply the missing ones.
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";
    $SAPserverO    ||= SAPserver->new()    or die "Could not get a new SAPserver\n";

    aligns_with_md5ID( $ALITREserverO, peg_to_md5( $SAPserverO, $_[0] || '' ) );
}


#-------------------------------------------------------------------------------
#  Run through a list of peg IDs, and produce a hash table that maps
#  each pegID to a list of alignments that contain it.
#
#   \%aligns_of_pegID = pegIDs_to_aligns( [@servers,] @pegIDs );
#
#    @servers can be $ALITREserverO and/or $SAPserverO.
#-------------------------------------------------------------------------------
sub pegIDs_to_aligns
{
    #  I will take neither, either or both objects in any order.
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $SAPserverO    = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    #  And then I will supply the missing ones.
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";
    $SAPserverO    ||= SAPserver->new()    or die "Could not get a new SAPserver\n";

    my @pegIDs = @_;
    my $md5s_of_pegIDs = pegs_to_md5( $SAPserverO, @pegIDs );
    my @md5IDs = values %$md5s_of_pegIDs;
    my $aligns_of_md5ID = md5IDs_to_aligns( $ALITREserverO, @md5IDs );
    my %aligns_of_pegID = map { $_ => $aligns_of_md5ID->{ $md5s_of_pegIDs->{ $_ } } } @pegIDs;
    
    \%aligns_of_pegID;
}


#-------------------------------------------------------------------------------
#  pegIDs in a given alignment
#
#    @pegIDs = pegIDs_in_align( [@servers,] $alignID );
#   \@pegIDs = pegIDs_in_align( [@servers,] $alignID );
#
#       @servers can be $ALITREserverO and/or $SAPserverO.
#-------------------------------------------------------------------------------
sub pegIDs_in_align
{
    #  I will take neither, either or both objects in any order.
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $SAPserverO    = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    #  And then I will supply the missing ones.
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";
    $SAPserverO    ||= SAPserver->new()    or die "Could not get a new SAPserver\n";

    my @md5s = md5IDs_in_align( $ALITREserverO, $_[0] || '' );
    #  Do a batch translation
    my $md5_to_pegs = md5s_to_pegs( $SAPserverO, @md5s );
    my %seen_peg;
    my @pegIDs = grep { ! $seen_peg{ $_ }++ }
                 map  { @{ $md5_to_pegs->{ $_ } || [] } }
                 @md5s;

    wantarray ? @pegIDs : \@pegIDs;
}


#-------------------------------------------------------------------------------
#  Alignment data with pegIDs:
#
#   \@seqs               = peg_alignment_by_ID( [@servers,] $alignID );
# ( \@seqs, \%metadata ) = peg_alignment_by_ID( [@servers,] $alignID );
#
#       $metadata{ $pegID } = [ $peg_length, $trim_beg, $trim_end, $location_string ]
#       @servers can be $ALITREserverO and/or $SAPserverO.
#-------------------------------------------------------------------------------
sub peg_alignment_by_ID
{
    #  I will take neither, either or both objects in any order.
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $SAPserverO    = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    #  And then I will supply the missing ones.
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";
    $SAPserverO    ||= SAPserver->new()    or die "Could not get a new SAPserver\n";

    ( my $alignID = $_[0] || '' ) =~ s/^aln//;
    my ( $md5_alignment, $md5_metadata ) = md5_alignment_by_ID( $ALITREserverO, $alignID );
    $md5_alignment && $md5_metadata or return wantarray ? () : undef;
    my ( $peg_metadata, $md5ID_to_fidIDs_map ) = map_md5_to_fid( $SAPserverO, $md5_metadata );
    my $peg_alignment = md5_align_to_fid_align( $md5_alignment, $md5ID_to_fidIDs_map );

    wantarray ? ( $peg_alignment, $peg_metadata ) : $peg_alignment;
}


#-------------------------------------------------------------------------------
#  Alignment sequence metadata with pegIDs:
#
#   \%metadata = peg_alignment_metadata( [@servers,] $alignID );
#
#       $metadata{ $pegID } = [ $peg_length, $trim_beg, $trim_end, $location_string ]
#       @servers can be $ALITREserverO and/or $SAPserverO.
#-------------------------------------------------------------------------------
sub peg_alignment_metadata
{
    #  I will take neither, either or both objects in any order.
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $SAPserverO    = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    #  And then I will supply the missing ones.
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";
    $SAPserverO    ||= SAPserver->new()    or die "Could not get a new SAPserver\n";

    ( my $alignID = $_[0] || '' ) =~ s/^aln//;
    my $md5_metadata = md5_alignment_metadata( $ALITREserverO, $alignID ) || {};
    my ( $peg_metadata, undef ) = map_md5_to_fid( $SAPserverO, $md5_metadata );

    $peg_metadata;
}


#-------------------------------------------------------------------------------
#  role-based alignments:
#-------------------------------------------------------------------------------
#  Alignments with a given role
#
#    @alignID_count = aligns_with_role( [@servers,] $role );
#   \@alignID_count = aligns_with_role( [@servers,] $role );
#
#    Servers can be FIG, ALITREserver and/or SAPserver. If FIG is supplied,
#    the ids are from the local SEED, otherwise they are from the SAPserver.
#-------------------------------------------------------------------------------
sub aligns_with_role
{
    my $fig;
    #  I will take neither, either or both objects in any order.
    $fig           = shift if eval { $_[0]->isa( 'FIG' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $fig           = shift if eval { $_[0]->isa( 'FIG' ) };
    $SAPserverO    = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $fig           = shift if eval { $_[0]->isa( 'FIG' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $fig           = shift if eval { $_[0]->isa( 'FIG' ) };
    #  And then I will supply the missing ones.
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";

    my $role = shift || '';
    length( $role ) or return wantarray ? () : [];

    my @md5s;
    if ( $fig )
    {
        my @pegIDs = $fig->prots_for_role( $role );
        my %md5s = map { $_ => 1 } values %{ $fig->md5_of_peg_bulk( \@pegIDs ) || {} };
        @md5s = keys %md5s;
    }
    else
    {
        $SAPserverO ||= SAPserver->new() or return wantarray ? () : [];
        my $pegIDs = ( $SAPserverO->occ_of_role( -roles => [$role] ) || {} )->{ $role } || [];
        my %md5s = map { $_ => 1 } values %{ pegs_to_md5( $SAPserverO, @$pegIDs ) || {} };
        @md5s = keys %md5s;
    }

    my %aligncnt;
    my $alignH = $ALITREserverO->aligns_with_md5ID( -ids => \@md5s ) || {};
    foreach my $md5 ( keys %$alignH )
    {
        foreach ( @{ $alignH->{ $md5 } || [] } ) { $aligncnt{ $_ }++ }
    }

    my @cnts = sort { $b->[1] <=> $a->[1] || $a->[0] cmp $b->[0] }
               map  { [ $_, $aligncnt{$_} ] }
               keys %aligncnt;

    wantarray ? @cnts : \@cnts;
}


#===============================================================================
#  Trees
#===============================================================================
#
#    @treeIDs = all_treeIDs( [$ALITREserverO] );
#   \@treeIDs = all_treeIDs( [$ALITREserverO] );
#
#-------------------------------------------------------------------------------
sub all_treeIDs
{
    $ALITREserverO   = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new()
        or die "Could not get a new ALITREserver\n";

    my $treeIDs = $ALITREserverO->all_treeIDs() || [];

    wantarray ? @$treeIDs : $treeIDs;
}


#-------------------------------------------------------------------------------
#  md5 based trees:
#-------------------------------------------------------------------------------
#
#    @treeIDs = trees_with_md5ID( [$ALITREserverO,] $md5 );
#   \@treeIDs = trees_with_md5ID( [$ALITREserverO,] $md5 );
#
#-------------------------------------------------------------------------------
sub trees_with_md5ID
{
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new()
        or die "Could not get a new ALITREserver\n";

    my $md5     = $_[0] || '';
    my $treeH   = $ALITREserverO->trees_with_md5ID( -ids => [ $md5 ] ) || {};
    my $treeIDs = $treeH->{ $md5 } || [];

    wantarray ? @$treeIDs : $treeIDs;
}

#-------------------------------------------------------------------------------
#
#    @md5IDs = md5IDs_in_tree( [$ALITREserverO,] $treeID );
#   \@md5IDs = md5IDs_in_tree( [$ALITREserverO,] $treeID );
#
#-------------------------------------------------------------------------------
sub md5IDs_in_tree
{
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new()
        or die "Could not get a new ALITREserver\n";

    ( my $treeID = $_[0] || '' ) =~ s/^tree//;
    my $md5H   = $ALITREserverO->md5IDs_in_tree( -ids => [ $treeID ] ) || [];
    my $md5IDs = $md5H->{ $treeID } || [];

    wantarray ? @$md5IDs : $md5IDs;
}

#-------------------------------------------------------------------------------
#  Run through a list of MD5 IDs, and produce a hash table that maps
#  each MD5 ID to a list of trees that contain it.
#
#   \%trees_of_md5ID = md5IDs_to_trees( [$ALITREserverO,] @md5IDs );
#
#-------------------------------------------------------------------------------
sub md5IDs_to_trees
{
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new()
        or die "Could not get a new ALITREserver\n";

    $ALITREserverO->trees_with_md5ID( -ids => \@_ ) || {};
}

#-------------------------------------------------------------------------------
#  Run through a list of peg IDs, and produce a hash table that maps
#  each pegID to a list of trees that contain it.
#
#   \%trees_of_pegID = pegIDs_to_trees( [@servers,] @pegIDs );
#
#    @servers can be $ALITREserverO and/or $SAPserverO.
#-------------------------------------------------------------------------------
sub pegIDs_to_trees
{
    #  I will take neither, either or both objects in any order.
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $SAPserverO    = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    #  And then I will supply the missing ones.
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";
    $SAPserverO    ||= SAPserver->new()    or die "Could not get a new SAPserver\n";

    my @pegIDs = @_;
    my $md5s_of_pegIDs = pegs_to_md5( $SAPserverO, @pegIDs );
    my @md5IDs = values %$md5s_of_pegIDs;
    my $trees_of_md5ID = md5IDs_to_trees( $ALITREserverO, @md5IDs );
    my %trees_of_pegID = map { $_ => $trees_of_md5ID->{ $md5s_of_pegIDs->{ $_ } } } @pegIDs;

    \%trees_of_pegID;
}

#-------------------------------------------------------------------------------
#
#      $tree              = md5_tree_by_ID( [$ALITREserverO,] $treeID );
#    ( $tree, $metadata ) = md5_tree_by_ID( [$ALITREserverO,] $treeID );
#
#-------------------------------------------------------------------------------
sub md5_tree_by_ID
{
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new()
        or die "Could not get a new ALITREserver\n";

    ( my $treeID = $_[0] || '' ) =~ s/^tree//;
    #***************************************************************************
    #  This function returns the tree and metadata pair for each treeID.
    #***************************************************************************
    my $treeH  = $ALITREserverO->md5_tree_by_ID( -ids => [ $treeID ] ) || {};
    my $treeL  = $treeH && $treeH->{ $treeID } ? $treeH->{ $treeID } : [];

    wantarray ? @$treeL : $treeL->[0];
}


#-------------------------------------------------------------------------------
#  peg based trees:
#-------------------------------------------------------------------------------
#
#    @treeIDs = trees_with_pegID( [@servers,] $pegID );
#   \@treeIDs = trees_with_pegID( [@servers,] $pegID );
#
#-------------------------------------------------------------------------------
sub trees_with_pegID
{
    #  I will take neither, either or both objects in any order.
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $SAPserverO    = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    #  And then I will supply the missing ones.
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";
    $SAPserverO    ||= SAPserver->new()    or die "Could not get a new SAPserver\n";

    trees_with_md5ID( $ALITREserverO, peg_to_md5( $SAPserverO, $_[0] || '' ) );
}

#-------------------------------------------------------------------------------
#
#    @pegIDs = pegIDs_in_tree( [@servers,] $treeID );
#   \@pegIDs = pegIDs_in_tree( [@servers,] $treeID );
#
#-------------------------------------------------------------------------------
sub pegIDs_in_tree
{
    #  I will take neither, either or both objects in any order.
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $SAPserverO    = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    #  And then I will supply the missing ones.
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";
    $SAPserverO    ||= SAPserver->new()    or die "Could not get a new SAPserver\n";

    ( my $treeID = $_[0] || '' ) =~ s/^tree//;
    my @md5s = md5IDs_in_tree( $ALITREserverO, $treeID );
    #  Do a batch translation
    my $md5_to_pegs = md5s_to_pegs( $SAPserverO, @md5s );
    my %seen_peg;
    my @pegIDs = grep { ! $seen_peg{ $_ }++ }
                 map  { @{ $md5_to_pegs->{ $_ } || [] } }
                 @md5s;

    wantarray ? @pegIDs : \@pegIDs;
}

#-------------------------------------------------------------------------------
#
#      $tree              = peg_tree_by_ID( [$servers,] $treeID );
#    ( $tree, $metadata ) = md5_tree_by_ID( [$servers,] $treeID );
#
#-------------------------------------------------------------------------------
sub peg_tree_by_ID
{
    #  I will take neither, either or both objects in any order.
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $SAPserverO    = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    #  And then I will supply the missing ones.
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";
    $SAPserverO    ||= SAPserver->new()    or die "Could not get a new SAPserver\n";

    ( my $treeID = $_[0] || '' ) =~ s/^tree//;
    my ( $md5_tree, $md5_metadata ) = md5_tree_by_ID( $ALITREserverO, $treeID );
    $md5_tree && $md5_metadata or return wantarray ? () : undef;
    
    my ( $peg_metadata, $md5ID_to_pegIDs_map ) = map_md5_to_fid( $SAPserverO, $md5_metadata );
    
    my $peg_tree = md5_tree_to_fid_tree( $md5_tree, $md5ID_to_pegIDs_map );

    wantarray ? ( $peg_tree, $peg_metadata ) : $peg_tree;
}


#===============================================================================
#  Functions for interconverting alignments and trees that md5-based ids and
#  fid-based ids.  Because the md5 id is based on the sequences, multiple
#  fids can have the same md5 id.  These are reduced to a single instance on
#  conversion to md5, and expanded to all known corresponding fids on conversion
#  back to fids.
#
#      (\@md5_align, \%md5_metadata) = fid_align_to_md5_align($sap,  \@fid_align, \%fid_metadata, $relaxed );
#      (\@fid_align, \%fid_metadata) = md5_align_to_fid_align($sap,  \@md5_align, \%md5_metadata, $relaxed );
#       $md5_tree  = fid_tree_to_md5_tree($sap,  $fid_tree, $relaxed );
#       $fid_tree  = md5_tree_to_fid_tree($sap,  $md5_tree, $relaxed );
#
#  sap              Sapling database object
#  @fid_align       An alignment, as fid_definition_sequence triples.
#  @md5_align       An alignment, as md5_definition_sequence triples.
#  %md5_metadata    hash mapping sequence IDs to MD5s with additional metadata
#                   (sequence id => [md5, len, beg, end, locations])
#  %fid_metadata    hash mapping feature IDs to sequence IDs with additional metadata
#                   (unique id => [fid, len, beg, end, locations])
#  $fid_tree        A gjonewick tree structure with fid ids.
#  $md5_tree        A gjonewick tree structure with md5 ids.
#  $relaxed         If set to a true value, untranslatable ids are passed through,
#                       rather than deleted.
#===============================================================================

## Run through the metadata associated with a FID-based alignment or tree and
## produce the md5 metadata and a map that translates FID IDs to MD5 IDs.

sub map_fid_to_md5
{
    my ( $sap, $fid_metadata, $relaxed ) = @_;
    $fid_metadata && ref( $fid_metadata ) eq 'HASH'
        or return ();

    my %md5_metadata;
    my %md5_covered; # maps each md5 to a list of the covered locations [beg, end]
    my %fidID_to_md5ID_map;

    foreach my $id ( keys %$fid_metadata )
    {
        my $fidMeta = $fid_metadata->{$id};
        my ($fid, $len, $beg, $end, $locations) = @$fidMeta;
        my $md5 = peg_to_md5( $sap, $fid );
        $md5 = $fid if ! $md5 && $relaxed;
        next if ! $md5;
        my $cover = $md5_covered{$md5} || [];
        my $found;
        foreach my $b_e (@$cover) {
            if (overlaps(@$b_e, $beg, $end)) {
                $found = 1;
                last;
            }
        }
        next if $found;
        my $md5ID = @$cover ? "$md5-" . scalar @$cover : $md5;
        push @{$md5_covered{$md5}}, [$beg, $end];
        $md5_metadata{$md5ID} = [$md5, $len, $beg, $end, $locations];
        $fidID_to_md5ID_map{$id} = $md5ID;
    }

    return (\%md5_metadata, \%fidID_to_md5ID_map);
}


sub fid_align_to_md5_align
{
    my ( $fid_align, $fidID_to_md5ID_map ) = @_;
    $fid_align && ref( $fid_align ) eq 'ARRAY' &&
        $fidID_to_md5ID_map && ref( $fidID_to_md5ID_map ) eq 'HASH'
        or return ();

    my @md5_align;

    foreach ( @$fid_align )
    {
        my $id = $_->[0];
        my $md5ID = $fidID_to_md5ID_map->{$id};
        next if ! $md5ID;
        push @md5_align, [ $md5ID, $_->[1], $_->[2] ];
    }

    return \@md5_align;
}


sub map_md5_to_fid
{
    my ( $sap, $md5_metadata, $relaxed ) = @_;
    $md5_metadata && ref( $md5_metadata ) eq 'HASH'
        or return ();

    my %fid_metadata;
    my %fids_seen;
    my %md5ID_to_fidIDs_map;

    my $md5_to_pegs = md5s_to_pegs( $sap, keys %$md5_metadata );

    foreach my $md5ID ( keys %$md5_metadata )
    {
        my $md5Metadata = $md5_metadata->{$md5ID};
        my ($md5, $len, $beg, $end, $location) = @$md5Metadata;
        my @fids = @{ $md5_to_pegs->{ $md5 } || [ ] };
        @fids = ( $md5 ) if ! @fids && $relaxed;
        foreach my $fid ( @fids )
        {
            my $fidID = $fid;
            if ($fids_seen{$fid}++) {
                $fidID = "$fid-" . $fids_seen{$fid};
            }
            $fid_metadata{$fidID} = [$fid, $len, $beg, $end, $location];
            push @{$md5ID_to_fidIDs_map{$md5ID}}, $fidID;
        }
    }

    return (\%fid_metadata, \%md5ID_to_fidIDs_map);
}


sub md5_align_to_fid_align
{
    my ( $md5_align, $md5ID_to_fidIDs_map ) = @_;
    $md5_align && ref( $md5_align ) eq 'ARRAY' && $md5ID_to_fidIDs_map &&
        ref( $md5ID_to_fidIDs_map ) eq 'HASH'
        or return ();

    my @fid_align;
    my %fid_metadata;

    foreach ( @$md5_align )
    {
        my $md5ID  = $_->[0];
        my @fidIDs = @{ $md5ID_to_fidIDs_map->{ $md5ID } || [] };
        foreach my $fidID ( @fidIDs )
        {
            push @fid_align, [ $fidID, $_->[1], $_->[2] ];
        }
    }

    return \@fid_align;
}


sub fid_tree_to_md5_tree
{
    my ( $fid_tree, $fidID_to_md5ID_map ) = @_;
    $fid_tree && ref( $fid_tree ) eq 'ARRAY' &&
        $fidID_to_md5ID_map && ref( $fidID_to_md5ID_map ) eq 'HASH'
        or return undef;

    gjonewicklib::newick_relabel_tips( gjonewicklib::newick_subtree( $fid_tree, keys %$fidID_to_md5ID_map ), $fidID_to_md5ID_map );
}


sub md5_tree_to_fid_tree
{
    my ( $md5_tree, $md5ID_to_fidIDs_map ) = @_;
    $md5_tree && ref( $md5_tree ) eq 'ARRAY' &&
        $md5ID_to_fidIDs_map && ref( $md5ID_to_fidIDs_map ) eq 'HASH'
        or return ();

    my @tips = gjonewicklib::newick_tip_list( $md5_tree );
    @tips or return undef;
    
    my $prune = 0;
    foreach my $md5ID ( @tips )
    {
        $prune = 1 if (! $md5ID_to_fidIDs_map->{$md5ID});
    }

    $md5_tree = gjonewicklib::newick_subtree( $md5_tree, [ keys %$md5ID_to_fidIDs_map ] ) if $prune;
    return expand_duplicate_tips( gjonewicklib::copy_newick_tree( $md5_tree ), $md5ID_to_fidIDs_map );
}


sub md5_align_and_tree_to_fid_version
{
    my ($sap, $md5_align, $md5_tree, $md5_metadata, $relaxed) = @_;
    my ($fid_metadata, $md5ID_to_fidIDs_map) = map_md5_to_fid($sap, $md5_metadata, $relaxed);
    my $fid_tree = md5_tree_to_fid_tree($md5_tree, $md5ID_to_fidIDs_map);
    my $fid_align = md5_align_to_fid_align($md5_align, $md5ID_to_fidIDs_map);
    return ($fid_align, $fid_tree, $fid_metadata);
}


sub fid_align_and_tree_to_md5_version
{
    my ($sap, $fid_align, $fid_tree, $fid_metadata, $relaxed) = @_;
    my ($md5_metadata, $fidID_to_md5ID_map) = map_fid_to_md5($sap, $fid_metadata, $relaxed);
    my $md5_tree = fid_tree_to_md5_tree($fid_tree, $fidID_to_md5ID_map);
    my $md5_align = fid_align_to_md5_align($fid_align, $fidID_to_md5ID_map);
    return ($md5_align, $md5_tree, $md5_metadata);
}


#-------------------------------------------------------------------------------
#
#  Return TRUE if the overlap between two sets of coordinates is sufficient to
#  treat them as essentially the same.
#
#   $okFlag = overlaps($beg1, $end1, $beg2, $end2);
#
#-------------------------------------------------------------------------------
sub overlaps
{
    # Get the parameters.
    my ( $beg1, $end1, $beg2, $end2 ) = @_;
    # Compute the number of overlapping residues.
    my $over = my_min( $end1, $end2 ) - my_max( $beg1, $beg2 ) + 1;
    # Return TRUE if the overlap is 80% of the shorter length.
    return $over >= 0.8 * my_min($end1 - $beg1 + 1, $end2 - $beg2 + 1);
}

sub my_min { $_[0] < $_[1] ? $_[0] : $_[1] }
sub my_max { $_[0] > $_[1] ? $_[0] : $_[1] }


#-------------------------------------------------------------------------------
#  Use a hash to relabel, and potentially expand the tips in a newick tree.
#
#  $node = expand_duplicate_tips( $node, \%new_names )
#
#-------------------------------------------------------------------------------
sub expand_duplicate_tips
{
    my ( $node, $new_names ) = @_;

    my @desc = gjonewicklib::newick_desc_list( $node );

    if ( @desc )
    {
        foreach ( @desc ) { expand_duplicate_tips( $_, $new_names ) }
    }
    else
    {
        my $new;
        if ( gjonewicklib::node_has_lbl( $node )
          && defined( $new = $new_names->{ gjonewicklib::newick_lbl( $node ) } )
           )
        {
            my @new = @$new;
            if ( @new == 1 )
            {
                gjonewicklib::set_newick_lbl( $node, $new[0] );
            }
            elsif ( @new > 1 )
            {
                gjonewicklib::set_newick_desc_ref( $node, [ map { [ [], $_, 0 ] } @new ] );
                gjonewicklib::set_newick_lbl( $node, undef );
            }
        }
    }

    $node;
}


#-------------------------------------------------------------------------------
#  Support for counting distinct roles in alignments and trees:
#
#    @role_count_pairs = roles_in_align( [$SAPserverO,] $alignID );
#   \@role_count_pairs = roles_in_align( [$SAPserverO,] $alignID );
#
#    $role  = majority_role_in_align( [$SAPserverO,] $alignID );
#
#    @role_count_pairs = roles_in_tree( [$SAPserverO,] $treeID );
#   \@role_count_pairs = roles_in_tree( [$SAPserverO,] $treeID );
#
#    $role  = majority_role_in_tree( [$SAPserverO,] $treeID );
#
#-------------------------------------------------------------------------------

sub roles_in_align
{
    #  I will take neither, either or both objects in any order.
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $SAPserverO    = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };

    #  And then I will supply the missing ones.
    $SAPserverO    ||= SAPserver->new()    or die "Could not get a new SAPserver\n";
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";

    ( my $alignID = $_[0] || '' ) =~ s/^aln//;
    my @md5s = md5IDs_in_align( $ALITREserverO, $alignID );
    my $md5_function = md5s_to_functions( $SAPserverO, @md5s ) || {};

    # my %cnt; 
    # for my $function ( map { $md5_function->{ $_ } } @md5s )
    # {
    #     next unless defined $function && $function =~ /\S/;
    #     $cnt{ $_ }++ for SeedUtils::roles_of_function( $function );
    # }

    my %func_cnt;
    for ( @md5s ) { $func_cnt{ $md5_function->{ $_ } }++ };

    my %cnt; 
    for my $function ( keys %func_cnt )
    {
        next unless defined $function && $function =~ /\S/;
        $cnt{ $_ } += $func_cnt{ $function } for SeedUtils::roles_of_function( $function );
    }

    my @pairs = sort { $b->[1] <=> $a->[1] || lc $a->[0] cmp lc $b->[0] }
                map { [ $_, $cnt{ $_ } ] }
                keys %cnt;

    wantarray ? @pairs : \@pairs;
}


sub roles_in_tree
{
    #  I will take neither, either or both objects in any order.
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $SAPserverO    = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };

    #  And then I will supply the missing ones.
    $SAPserverO    ||= SAPserver->new()    or die "Could not get a new SAPserver\n";
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";

    ( my $treeID = $_[0] || '' ) =~ s/^tree//;
    my @md5s = md5IDs_in_tree( $ALITREserverO, $treeID );
    my $md5_function = md5s_to_functions( $SAPserverO, @md5s ) || {};

    my %cnt; 
    for my $function ( map { $md5_function->{ $_ } } @md5s )
    {
        next unless defined $function && $function =~ /\S/;
        $cnt{ $_ }++ for SeedUtils::roles_of_function( $function );
    }

    my @pairs = sort { $b->[1] <=> $a->[1] || lc $a->[0] cmp lc $b->[0] }
                map { [ $_, $cnt{ $_ } ] }
                keys %cnt;

    wantarray ? @pairs : \@pairs;
}


sub majority_role_in_align
{
    my ($first_pair) = roles_in_align( @_ );
    $first_pair ? $first_pair->[0] : undef;
}


sub majority_role_in_tree
{
    my ($first_pair) = roles_in_tree( @_ );
    $first_pair ? $first_pair->[0] : undef;
}


my %peg_to_function;
sub pegs_to_function
{
    $SAPserverO = shift if eval { $_[0]->isa( 'SAPserver' ) };
    
    #  Remove the pegIDs already in the cache
    my @pegID = grep { ! exists $peg_to_function{ $_ } } @_;
    if ( @pegID )
    {
        #  Get the remaining pegIDs from the Sapling:
        $SAPserverO ||= SAPserver->new() or die "Could not get a new SAPserver\n";
        my $functionH = $SAPserverO->ids_to_functions( -ids => \@pegID ) || {};
        foreach ( @pegID )
        {
            my $func = $functionH->{ $_ };
            $peg_to_function{ $_ } = defined $func ? $func : '';
        }
    }
    
    #  Return the whole cache
    wantarray ? %peg_to_function : \%peg_to_function;
}


#-------------------------------------------------------------------------------
#  Support for a function lookup for md5 IDs:
#
#    \%md5_function = md5s_to_functions( [$SAPserverO,] @md5s );
#
#  When an md5ID maps to multiple pegIDs, only one pegID's function is counted.
#-------------------------------------------------------------------------------

my %md5_function;
sub md5s_to_functions
{
    $SAPserverO = shift if eval { $_[0]->isa( 'SAPserver' ) };

    #  Remove the md5IDs already in the cache
    my @new_md5s = grep { ! exists $md5_function{ $_ } } @_;

    if ( @new_md5s )
    {
        #  Batch lookup of md5 ids returns hash of all currently known translations
        $SAPserverO ||= SAPserver->new() or die "Could not get a new SAPserver\n";
        my $md5_to_pegs = md5s_to_pegs( $SAPserverO, @new_md5s );

        #  Just get first peg (this is dangerous, but will become more consistent in future) 
        my %md5_to_one_peg = map { $md5_to_pegs->{ $_ } ? ( $_ => $md5_to_pegs->{ $_ }->[0] ) : () } @new_md5s;
        my @new_pegs       = values %md5_to_one_peg;
        my $new_peg_func   = pegs_to_function( $SAPserverO, @new_pegs );

        #  Add new functions to those known, converting undef to blank
        foreach ( @new_md5s )
        {
            my $func = $new_peg_func->{ $md5_to_one_peg{ $_ } };
            $md5_function{ $_ } = defined $func ? $func : '';
        }
    }

    #  Return hash of all known
    wantarray ? %md5_function : \%md5_function;
}


#-------------------------------------------------------------------------------
#  Support for getting projections between MD5 IDs:
#
#   ( \@sets, \%metadata ) = min_paralog_trees( [$SAPserverO,] $treeID, $opts );
#     \@list_of_tips       = min_paralog_trees( [$SAPserverO,] $treeID, $opts );
#
#  Options:
#
#     fract_cover  => $fract (D=0.75)    # min fraction of median protein length
#     genome_excl  => \@gids             # genomes to exclude 
#     genome_incl  => \@gids             # genomes to include
#     min_tip_cnt  => $threshold (D=3)   # min size of subtrees to keep
#
#   @sets = [ [ subtree_1_tip_1, subtree_1_tip_2 ... ], [ subtree_2_tip1, ... ] ... ]
#   $metadata{ $tip } = [ peg, len, trim_beg, trim_end, location_string ]
#
#   @list_of_tips = [ [ tip, subtree_index, peg, len, trim_beg, trim_end, location_string ], ... ]
#
#-------------------------------------------------------------------------------

sub min_paralog_trees
{
    $SAPserverO = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $SAPserverO ||= SAPserver->new() or die "Could not get a new SAPserver\n";

    my ($treeID, $opts) = @_;
    
    my $min_tip_cnt = $opts->{min_tip_cnt} || 3;
    my $genome_excl = $opts->{genome_excl};
    my $genome_incl = $opts->{genome_incl};
    my $fract_cover = $opts->{fract_cover} || 0.75;

    my (%gid_excl, %gid_incl);
    %gid_excl = map { $_ => 1 } @$genome_excl if $genome_excl && ref $genome_excl eq 'ARRAY';
    %gid_incl = map { $_ => 1 } @$genome_incl if $genome_incl && ref $genome_incl eq 'ARRAY';

    my ($tree, $metadata) = peg_tree_by_ID($treeID);
    if (!$tree) {
        wantarray ? return ([], {}) : return [];
    }

    my @lens =  map { my $md = $metadata->{$_};
                      [ $_, $md->[3] - $md->[2] ] } keys %$metadata;
    
    my $median = ( sort { $a <=> $b } map { $_->[1] } @lens )[ int(@lens/2) ];
    my @keep = grep { my $gid = SeedUtils::genome_of($_);
                      ( ! $genome_incl || $gid_incl{ $gid } ) && ! $gid_excl{ $gid } } 
               map { $_->[1] >= $fract_cover * $median ? $_->[0] : () } @lens;

    if (@keep < 3) {
        wantarray ? return ([], {}) : return [];
    }
    
    $tree = gjonewicklib::reroot_newick_to_midpoint_w($tree);
    $tree = gjonewicklib::newick_subtree($tree, \@keep);

    my ($sets) = min_paralog_tree_1($tree);
    $sets = [ sort { @$b <=> @$a } grep { @$_ >= $min_tip_cnt } @$sets ];

    if (wantarray) {
        return ($sets, $metadata);
    } else {
        my @list;
        my $i;
        for my $set (@$sets) {
            ++$i;
            for my $tip (@$set) {
                push @list, [ $tip, $i, @{$metadata->{$tip}} ];
            }
        }
        return \@list;
    }
}

sub min_paralog_tree_1 
{
    my ($node) = @_;
    my $sets = [ ];
    my $gids = [ ];

    my (@list, %gidcnt);
    my $merge = 1;
    my @desc = gjonewicklib::newick_desc_list( $node );
    
    if ( @desc ) {
        for my $n ( @desc ) {
            my ( $sets2, $gids2 ) = min_paralog_tree_1( $n );
            $merge = 0 if @$sets2 > 1;
            if ( $merge && $gids2 ) {
                for ( @$gids2 ) {
                    if ( $gidcnt{ $_ }++ ) {
                        # print STDERR "$_\n";
                        $merge = 0; last;
                    }
                }
            }
            for (@$sets2) { push @$sets, $_; }
        }
        if ($merge) {
            $gids = [ keys %gidcnt ];
            my @merged_set;
            for my $set (@$sets) { 
                for (@$set) { push @merged_set, $_; }
            }
            $sets = [ \@merged_set ];
        }

    } else {
        my $tip   = gjonewicklib::newick_lbl( $node );
        my ($gid) = ($tip =~ /(\d+\.\d+)/);
        if ($gid) {
            $sets = [ [ $tip ] ];
            $gids = [ $gid ];
        }
    }

    ( $sets, $gids );
}


#-------------------------------------------------------------------------------
#  Support for getting projections between MD5 IDs:
#
#   \%md5_projections = get_md5_projections( [$SAPserverO,] @md5s [, \%opts] );
#
#  Options:
#
#     minScore  => $threshold (D=0)   # only get projections with scores greater than threshold
#     details   => bool       (D=0)   # get detailed projections if set, see below
#
#  details = 0:
#
#     { $md5_1 => [ $md5_1_a, $md5_1_b, ... ],
#       $md5_2 => [ $md5_2_a, $md5_2_b, ... ],
#                  ... };
#  
#  details = 1:
#
#     { $md5_1 => [ [ $md5_1_a, $context_1_a, $ident_1_a, $score_1_a ],
#                   [ $md5_1_b, $context_1_b, $ident_1_b, $score_1_b ],
#                   ... ],
#       $md5_2 => [ [ $md5_2_a, $context_2_a, $ident_2_a, $score_2_a ],
#                   [ $md5_2_b, $context_2_b, $ident_2_b, $score_2_b ],
#                   ... ],
#                 ... };
#
#-------------------------------------------------------------------------------

sub get_md5_projections
{
    #  I will take neither, either or both objects in any order.
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $SAPserverO    = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };

    #  And then I will supply the missing ones.
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";

    my $opts = ref $_[-1] eq 'HASH' ? pop @_ : {};
    my @md5s = @_ or die "Empty list of MD5 IDs\n";;

    my $minScore = $opts->{ minScore } || $opts->{ min_score } || 0;
    my $details  = $opts->{ details }  || $opts->{ full }      || 0;

    $ALITREserverO->get_projections( -ids => \@md5s, -minScore => $minScore, -details => $details );
}


#==============================================================================
#  Get data on alignments with a protein
#
#   @align_data = data_on_aligns_with_prot( [$ALITREserverO,] [$SAPserverO,] $fid [, \%opts] );
#   @align_data = data_on_aligns_with_prot( [$ALITREserverO,] [$SAPserverO,] $md5 [, \%opts] );
#
#   or with $opts->{fid} or $opts->{md5}
#
#   @align_data = data_on_aligns_with_prot( [$ALITREserverO,] [$SAPserverO,]      [  \%opts] );
#
#   Alignment data:
#
#   [ $alignID, $coverage, \@roles ]
#
#        $coverage = [ $beg, $end, $len ]
#       \@roles    = [ [ $role1, $cnt1 ], [ $role2, $cnt2 ], ... ]
#
#   Options:
#
#       coverage => $frac  # minimum fraction of protein residues included
#       fid      => $fid   # feature id, if not an argument
#       md5      => $md5   # protein sequence md5, if not an argument
#       noroles  => $bool  # omit roles and counts (faster)
#       roles    => $nrole # limit on number of roles returned
#       sap      => $sapO  # Sapling server object
#
#==============================================================================

sub data_on_aligns_with_prot
{
    #  I will take neither, either or both objects in any order.
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $SAPserverO    = shift if eval { $_[0]->isa( 'SAPserver' ) };
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };

    my $opts = $_[-1] && ref( $_[-1] ) eq 'HASH' ? pop
             : $_[ 0] && ref( $_[ 0] ) eq 'HASH' ? shift
             :                                     {};

    #  Fill in the missing server objects
    $SAPserverO    ||= $opts->{sap}
                   ||= SAPserver->new()    or die "Could not get a new SAPserver\n";
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";

    my $roles = defined $opts->{roles} ? $opts->{roles}
              : $opts->{noroles}       ? 0
              :                          1e9;

    my $min_cov = $opts->{ coverage } || 0;

    my $md5 = $_[0]                     ? ( valid_md5( $_[0] ) ? shift : peg_to_md5( $SAPserverO, $_[0] ) )
            : valid_md5( $opts->{md5} ) ? $opts->{md5}
            : $opts->{fid}              ? peg_to_md5( $SAPserverO, $opts->{fid} )
            :                             '';

    my %alignID_coverages = map { $_->[0] => $_ }
                            alignment_coverages_of_md5( $md5 );

    my @align_ids = keys %alignID_coverages;

    #  Compile data on available alignments

    my @aligns = ();
    foreach my $id ( @align_ids )
    {
        #  Work out the sequence coverage report:
        my ( undef, $md5, $len, $covers ) = @{ $alignID_coverages{$id} || [] };
        my $coverage;
        my $n_cov = 0;
        my $min_n = $min_cov * ( $len || 0 );
        if ( $len && $covers && @$covers )
        {
            $coverage = [ $covers->[0]->[0], $covers->[-1]->[1], $len ];
            foreach ( @$covers ) { $n_cov += $_->[1] - $_->[0] + 1 }
            next if $n_cov < $min_n;
        }
        else
        {
            $coverage = [];
        }

        #  Compile the role data:

        my @roles = $roles ? roles_in_align( $SAPserverO, $id ) : ();
        splice @roles, $roles if @roles > $roles;

        push @aligns, [ [ $id, $coverage, @roles ? \@roles : () ], $n_cov ];
    }

    @aligns = map  { $_->[0] }
              sort { $b->[1] <=> $a->[1] || $a->[0]->[0] cmp $b->[0]->[0] }
              @aligns;

    wantarray ? @aligns : \@aligns;
}


sub valid_md5 { $_[0] && $_[0] =~ /^[0-9A-Za-z]{32}$/ }

#-------------------------------------------------------------------------------
#  PSSM and rpsblast db
#
#    db_name = pssm_from_alignID( [$ALITREserverO,] $alignID, \%options )
#
#  The IDs of sequences in the database are their MD5s.
#
#  Options (all options are related to BlastInterface::alignment_to_pssm()):
#
#    ignore_msa_master => $bool     # ignore the master sequence when psiblast creates PSSM (D = 0)
#    ignoreMaster      => $bool     # ignore the master sequence when psiblast creates PSSM (D = 0)
#    max_sim           => $fract    # maximum identity of sequences in profile (D = no_limit)
#    min_sim           => [ $min_fract_ident, @IDs ]
#                                  # eliminate divergent sequences (D = undef)
#    msa_master_id     => $id      # ID of the sequence in in MSA for psiblast to use as a master
#    msa_master_idx    => 1-based index of the sequence in MSA for psiblast to use as a master
#    out_pssm          => output PSSM filename or handle (D = stdout)
#    outPSSM           => output PSSM filename or handle (D = stdout)
#    pseudo_master     => $bool    # create a consensus master sequence with all columns in alignment
#    pseudoMaster      => $bool    # create a consensus master sequence with all columns in alignment
#    title             => $title   # title of the PSSM (D = alignID)
#
#-------------------------------------------------------------------------------

sub pssm_from_alignID
{
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";

    eval { require BlastInterface; }
        or warn "AlignsAndTreesServer::pssm_from_alignID: failed in require BlastInterface\n"
            and return undef;

    my ( $alignID, $opts ) = @_;

    my $align;
    $align = md5_alignment_by_ID( $ALITREserverO, $alignID )
        and @$align
            or warn "AlignsAndTreesServer::pssm_from_alignID: failed to retrieve alignment '$alignID'\n"
                and return undef;

    BlastInterface::alignment_to_pssm( $align, $opts );
}

# 
# 
# options: title

sub rps_db_from_alignIDs
{
    $ALITREserverO = shift if eval { $_[0]->isa( 'ALITREserver' ) };
    $ALITREserverO ||= ALITREserver->new() or die "Could not get a new ALITREserver\n";

    eval { require BlastInterface; }
        or warn "AlignsAndTreesServer::rps_db_from_alignIDs: failed in require BlastInterface\n"
            and return undef;

    my ( $alignIDs, $db, $opts ) = @_;

    $alignIDs && ref( $alignIDs ) eq 'ARRAY' && @$alignIDs
        or warn "AlignsAndTreesServer::rps_db_from_alignIDs: call without alignment IDs\n"
            and return undef;

    $db or warn "AlignsAndTreesServer::rps_db_from_alignIDs: undefined database file name\n"
        and return undef;

    $opts = {} unless $opts && ref($opts) eq 'HASH';
    
    my %tmp_dir_opts = ( base => 'rps_db_tmp', %$opts );
    my ( $tmp_dir, $save_dir ) = SeedAware::temporary_directory( \%tmp_dir_opts );

    my @files;
    my %pssm_opts = %$opts;
    foreach my $alignID (@$alignIDs)
    {
        $alignID =~ s/^(\d{8})/aln$1/;
        $pssm_opts{ title } = $alignID;
        $pssm_opts{ out_pssm } = "$tmp_dir/$alignID.pssm";
        my $file = pssm_from_alignID( $alignID, \%pssm_opts )
            or warn "AlignsAndTreesServer::rps_db_from_alignIDs: failed to build PSSM for '$alignID'\n"
                and next;
        push @files, $file;
    }
    @files
        or warn "AlignsAndTreesServer::rps_db_from_alignIDs: failed to build any PSSM\n"
            and return undef;

    my $title = $opts->{ title };
    if ( ! $title )
    {
        $title = $db; 
        $title =~ s|.*/||;
    }
    my %rps_opts = ( title => $title );
    BlastInterface::build_rps_db( \@files, $db, \%rps_opts );
}

1;
