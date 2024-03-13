#
# Copyright (c) 2003-2011 University of Chicago and Fellowship
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

package AlignsAndTrees;

#===============================================================================
#  perl functions for loading and accessing Alignments and Trees based on md5
#
#  Usage:  use AlignsAndTrees;
#
#    @alignIDs = all_alignIDs();
#   \@alignIDs = all_alignIDs();
#
#    @alignIDs = aligns_with_md5ID( $md5 );
#   \@alignIDs = aligns_with_md5ID( $md5 );
#
#   \%md5_alignIDs = md5IDs_to_aligns( $sap, @md5IDs );
#
#    @md5IDs   = md5IDs_in_align( $sap, $alignID );
#   \@md5IDs   = md5IDs_in_align( $sap, $alignID );
#
#   \@seqs                       = md5_alignment_by_ID( $sap, $alignID );
# ( \@seqs, \%md5_row_metadata ) = md5_alignment_by_ID( $sap, $alignID );
#   \%md5_row_metadata           = md5_alignment_metadata( $sap, $alignID );
#
#       $md5_row_metadata{ $seqID } = [ $md5, $peg_length, $trim_beg, $trim_end, $location_string ]
#
#   \%md5_row_metadata = alignment_metadata_by_md5( $sap, $alignID, @md5IDs );
#
#        $md5_row_metadata{ $seqID } = [ $md5ID, $peg_length, $trim_beg, $trim_end, $location_string ]
#
#    @metadata = alignments_metadata_by_md5( $sap, @md5IDs );
#   \@metadata = alignments_metadata_by_md5( $sap, @md5IDs );
#
#        @metadata = ( [ $alignID, $seqID, $md5, $peg_length, $trim_beg, $trim_end, $location_string ], ... )
#
#    @treeIDs  = all_treeIDs( );
#   \@treeIDs  = all_treeIDs( );
#
#    @treeIDs  = trees_with_md5ID( $md5 );
#   \@treeIDs  = trees_with_md5ID( $md5 );
#
#   \%md5_treeIDs = md5IDs_to_trees( $sap, @md5IDs );
#
#    @md5IDs   = md5IDs_in_tree( $treeID );
#   \@md5IDs   = md5IDs_in_tree( $treeID );
#
#    $tree              = md5_tree_by_ID( $treeID );
#  ( $tree, \%metadata) = md5_tree_by_ID( $treeID );
#   \%metadata          = md5_tree_metadata( $treeID );
#
#    @fids = md5_to_pegs($sap, $md5);
#   \@fids = md5_to_pegs($sap, $md5);
#
#    $md5  = peg_to_md5($sap, $peg);
#
#  ( \@md5_align, \%md5_metadata ) = fid_align_to_md5_align( $sap, \@fid_align, \%fid_metadata, $relaxed );
#  ( \@fid_align, \%fid_metadata ) = md5_align_to_fid_align( $sap, \@md5_align, \%md5_metadata, $relaxed );
#
#    $md5_tree = fid_tree_to_md5_tree( $sap,  $fid_tree, $relaxed );
#    $fid_tree = md5_tree_to_fid_tree( $sap,  $md5_tree, $relaxed );
#
#  ( \@fid_align, $fid_tree, \%fid_metadata) = md5_align_and_tree_to_fid_version( \@md5_align, $md5_tree, \%md5_metadata)
#  ( \@md5_align, $md5_tree, \%md5_metadata) = fid_align_and_tree_to_md5_version( \@fid_align, $fid_tree, \%fid_metadata)
#
#   $alignID = load_md5_alignment_and_tree( $sap, $md5_align, $md5_tree, \@seq_metadata, \@alignment_tree_metadata );
#   delete_md5_alignment_and_tree( $sap, $alignID );
#
#   @metadata = alignment_tree_metadata( $sap, $alignID );
#  \@metadata = alignment_tree_metadata( $sap, $alignID );
#
#   @metadata = ( [ alignment-method, alignment-parameters, alignment-properties,
#                   tree-method, tree-parameters, tree-properties ], ... )
#
#-------------------------------------------------------------------------------
#  Some internal subroutines:
#
#   Alignments and trees are in files in a data directory.  This locates the
#   directory so that they can be read and written.
#
#   $data_dir = locate_data_dir();
#
#   Some basic sanity checks on the ids in an alignment, tree and associated
#   per sequence metadata.
#
#   $okFlag = validate_alignment_and_tree_ids( $alignment, $tree, $sequenceMetadata);
#
#   In deciding if two sequences with the same MD5 are really the saome thing,
#   we ask if they overlap by >80% of the length of the shorter sequence.
#
#   $boolean = overlaps($beg1, $end1, $beg2, $end2);
#
#   An MD5 tree tip can be expanded to a large number of tree tips.  This
#   routine uses the new_names hash to relabel each tree tip to a list of
#   tips, creating a multifurcation of zero-length branches when expanding
#   one to several.
#      
#   $node = expand_duplicate_tips( $node, \%new_names )
#
#===============================================================================

use gjoseqlib    qw( read_fasta );
use gjonewicklib qw( read_newick_tree );
# use FIG_Config;
use Tracer;
use Storable qw( nfreeze thaw );
use SeedUtils;
use strict;

my $data_dir = undef;
sub locate_data_dir
{
    if ( ! $data_dir )
    {
        if ( $ENV{ ATNG } && -d $ENV{ ATNG } ) {
            $data_dir = $ENV{ ATNG };
        } else {
            require FIG_Config;           # Only necessary for a few functions
            $data_dir = "$FIG_Config::fig/ATNG";
        }

        if ( ! $data_dir || ! -d $data_dir ) {
            die "Could not locate directory of alignments and trees.\n";
        }
    }

    $data_dir;
}

#-------------------------------------------------------------------------------
#
# Load an MD5 alignment and tree.
#
#   $alignID = load_md5_alignment_and_tree($sap, $alignStruct, $treeStruct,
#                                          $sequenceMetadata, $alignmentMetadata);
#
# where $sap is the Sapling database object.
#       $alignStruct describes the alignment
#       $treeStruct describes the tree
#       $sequenceMetadata describes the relationship of each protein sequence to the alignment
#       $alignmentMetadata contains information about how the alignment and tree were computed
#
#   $treeStruct is a gjonewick tree structure
#   $alignStruct is a reference to an array of protein FASTA triples (id, comment, sequence)
#   $sequenceMetadata is a reference to an array of tuples containing the description of
#               how each protein sequence is used by the alignment (sequence-id, md5, len, begin, end, locations)
#   $alignmentMetadata is reference to a list of strings (alignment-method,
#               alignment-parameters, alignment-properties, tree-method,
#               tree-parameters, tree-properties)
#
#-------------------------------------------------------------------------------
sub load_md5_alignment_and_tree {
    # Get the parameters.
    my ($sap, $alignStruct, $treeStruct, $sequenceMetadata, $alignmentMetadata) = @_;
    # Insure the alignment is valid.
    if (! validate_alignment_and_tree_ids($alignStruct, $treeStruct, $sequenceMetadata)) {
        die "Invalid alignment data.\n";
    }
    # Insure we have a place to store the alignment flat files.
    my $dataDir = locate_data_dir();
    if (! $dataDir) {
        die "No alignment storage directory.\n";
    }
    # First, we compute the ID by looking in the database to see the current largest
    # ID number.
    my ($newID) = $sap->GetFlat('AlignmentTree', "ORDER BY AlignmentTree(id) DESC LIMIT 1",
                                [], "id");
    if (! $newID) {
        $newID = "00000001";
    } else {
        $newID++;
    }
    # Insert the alignment-tree root record.
    $sap->InsertObject('AlignmentTree',
                       id                   => $newID,
                       alignment_method     => $alignmentMetadata->[0],
                       alignment_parameters => $alignmentMetadata->[1],
                       alignment_properties => $alignmentMetadata->[2],
                       tree_method          => $alignmentMetadata->[3],
                       tree_parameters      => $alignmentMetadata->[4],
                       tree_properties      => $alignmentMetadata->[5]);
    # Now we store the alignment file in the data directory.
    my $alignFile = "$dataDir/ali$newID.fa";
    gjoseqlib::print_alignment_as_fasta($alignFile, $alignStruct);
    # Next we store the tree file in the same directory.
    my $treeFile = "$dataDir/tree$newID.nwk";
    gjonewicklib::writeNewickTree($treeStruct, $treeFile);
    # Then we store a frozen version of the tree in a subdirectory.
    my $frozenFile = "$dataDir/FROZEN/tree$newID.frozen";
    SeedUtils::verify_dir("$dataDir/FROZEN");
    open(F, ">$frozenFile") or die "Could not open $frozenFile";
    print F nfreeze($treeStruct);
    close(F);
    # Finally, we create the realtionships that connect the protein sequences to
    # the alignment in the database.
    for my $seqTuple (@$sequenceMetadata) {
        # Get the components of this sequence's tuple.
        my ($seqID, $md5, $len, $begin, $end, $locations) = @$seqTuple;
        # Insert a relationship record for this sequence/alignment pair.
        $sap->InsertObject('Aligns', from_link => $newID, sequence_id => $seqID,
                           to_link => $md5, len => $len, begin => $begin,
                           end => $end, properties => $locations);
    }
    # Return the alignment/tree ID.
    return $newID;
}

#-------------------------------------------------------------------------------
#
# Delete an md5 alignment and tree.
#
#   delete_md5_alignment_and_tree($sap, $alignID);
#
# where $sap is the Sapling database object
#       $alignID is the ID of the alignment/tree to delete
#
#-------------------------------------------------------------------------------
sub delete_md5_alignment_and_tree {
    # Get the parameters.
    my ($sap, $alignID) = @_;
    # Get the data directory.
    my $dataDir = locate_data_dir();
    if (! $dataDir) {
        die "No alignment/tree data directory present.\n";
    }
    # Delete the alignment file if it exists.
    my $alignFile = "$dataDir/ali$alignID.fa";
    if (-f $alignFile) {
        unlink $alignFile;
    }
    # Delete the tree file if it exists.
    my $treeFile = "$dataDir/tree$alignID.nwk";
    if (-f $treeFile) {
        unlink $treeFile;
    }
    # Delete the frozen tree file if it exists.
    my $frozenFile = "$dataDir/FROZEN/tree$alignID.frozen";
    if (-f $frozenFile) {
        unlink $frozenFile;
    }
    # Delete the alignment from the database.
    $sap->Delete(AlignmentTree => $alignID);
}

#-------------------------------------------------------------------------------
#
# Verify that the sequence IDs are consistent in the load data for an md5 alignment
# and tree.
#
#   $okFlag = validate_alignment_and_tree_ids($alignStruct, $treeStruct,
#                                          $sequenceMetadata);
#
# where $alignStruct describes the alignment (if undef, will not be checked)
#       $treeStruct describes the tree (if undef, will not be checked)
#       $sequenceMetadata describes the relationship of each protein sequence to the alignment
#
#   $treeStruct is a gjonewick tree structure
#   $alignStruct is a reference to an array of protein FASTA triples (id, comment, sequence) 
#   $sequenceMetadata is a reference to an array of tuples containing the description of
#               how each protein sequence is used by the alignment (sequence-id, md5, len, begin, end, locations)
#               or a reference to a hash mapping each sequence ID to the description
#               (sequence-id => [md5, len, begin, end, locations])
#-------------------------------------------------------------------------------
sub validate_alignment_and_tree_ids {
    # Get the parameters.
    my ($alignStruct, $treeStruct, $sequenceMetadata) = @_;
    # Assume we're valid unless we determine otherwise.
    my $retVal = 1;
    # Get a hash of the sequence IDs from the sequence metadata. Each of the other
    # objects must have exactly this set of sequence IDs.
    my $sequences;
    if (! $sequenceMetadata) {
        warn "No sequence metadata passed into alignment validation.\n";
        return 0;
    } elsif (ref $sequenceMetadata eq 'HASH') {
        # Hash passed in. Just copy it.
        $sequences = $sequenceMetadata;
    } elsif (ref $sequenceMetadata eq 'ARRAY') {
        # Array passed in. Create a hash from the sequence elements.
        $sequences = { map { $_->[0] => 1 } @$sequenceMetadata };
        # Verify that there are no duplicates in the sequence metadata.
        if (scalar(keys %$sequences) != scalar @$sequenceMetadata) {
            warn "Duplicate sequence ID found in sequence metadata for alignment.\n";
            return 0;
        }
    } else {
        warn "Unrecognized sequence data structure passed into alignment validator.\n";
        return 0;
    }
    my $sequenceCount = scalar keys %$sequences;
    # The sequence metadata appears valid.
    if ($alignStruct) {
        # We have an alignment. Now loop through the alignment structure
        # to insure its sequence IDs match. We'll use this hash to make sure there are
        # no duplicates.
        my %alignSeen;
        # Verify that the alignment is a valid array.
        if (ref $alignStruct ne 'ARRAY') {
            $retVal = 0;
            warn "Alignment is not an array.\n";
        } else {
            for my $alignRow (@$alignStruct) {
                my $sequence = $alignRow->[0];
                # Insure this sequence exists.
                if (! $sequences->{$sequence}) {
                    # Here it was not in the metadata.
                    $retVal = 0;
                    warn "Invalid sequence ID $sequence found in alignment FASTA.\n";
                } elsif ($alignSeen{$sequence}) {
                    # Here we've seen it twice.
                    $retVal = 0;
                    warn "Duplicate sequence ID $sequence found in alignment FASTA.\n";
                } else {
                    # Here the sequence is valid.
                    $alignSeen{$sequence} = 1;
                }
            }
            # Insure everything was found.
            if (scalar(keys %alignSeen) != $sequenceCount) {
                $retVal = 0;
                warn "Some sequence IDs were missing from the alignment FASTA.\n";
            }
        }
    }
    # Only check the tree structure if it exists.
    if ($treeStruct) {
        # Now we validate the tree. We'll use this hash to make sure there are
        # no duplicates in the tree.
        my %treeSeen;
        # Verify that the tree is a valid newick tree.
        if (ref $treeStruct ne 'ARRAY') {
            $retVal = 0;
            warn "Alignment tree is not an array.\n";
        } else {
            # Get the list of sequences in the tree.
            for my $sequence (gjonewicklib::newick_tip_list($treeStruct)) {
                # Insure this sequence exists.
                if (! $sequences->{$sequence}) {
                    # Here it was not in the metadata.
                    $retVal = 0;
                    warn "Invalid sequence ID $sequence found in alignment tree.\n";
                } elsif ($treeSeen{$sequence}) {
                    # Here we've seen it twice.
                    $retVal = 0;
                    warn "Duplicate sequence ID $sequence found in alignment tree.\n";
                } else {
                    # Here the sequence is valid.
                    $treeSeen{$sequence} = 1;
                }
            }
            # Insure everything was found.
            if (scalar(keys %treeSeen) != $sequenceCount) {
                $retVal = 0;
                warn "Some sequence IDs were missing from the alignment tree.\n";
            }
        }
    }
    # Return the determination indicator.
    return $retVal;
}

#-------------------------------------------------------------------------------
#
#    @alignIDs = all_alignIDs($sap);
#   \@alignIDs = all_alignIDs($sap);
#
#
# where $sap is the Sapling database object.
#-------------------------------------------------------------------------------
sub all_alignIDs
{
    my ( $sap ) = @_;
    my @ids = $sap ? $sap->GetFlat('AlignmentTree', "", [], "id") : ();
    wantarray ? @ids : \@ids;
}

#-------------------------------------------------------------------------------
#
#    @alignIDs = aligns_with_md5ID( $sap, $md5 );
#   \@alignIDs = aligns_with_md5ID( $sap, $md5 );
#
# where $sap is the Sapling database object.
#       $md5 is an MD5 protein sequence ID whose alignments are desired
#-------------------------------------------------------------------------------
sub aligns_with_md5ID
{
    my ( $sap, $md5 ) = @_;
    my @ids;
    if ( $sap && $md5 ) {
        @ids = $sap->GetFlat('Aligns', 'Aligns(to-link) = ?', [$md5],
                             'from-link');
    }
    wantarray ? @ids : \@ids;
}

#-------------------------------------------------------------------------------
#
#   \%md5_alignIDs = md5IDs_to_aligns( $sap, @md5IDs );
#
#     $md5_alignIDs{ $md5ID } = \@alignIDs_for_md5ID
#
# where $sap is the Sapling database object.
#       @md5IDs is an list of MD5 protein sequence IDs whose alignments are desired
#-------------------------------------------------------------------------------
sub md5IDs_to_aligns
{
    my ( $sap, @md5IDs ) = @_;
    $sap && @md5IDs or return {};

    # Declare the return hash.
    my $retVal = {};
    # Loop through the incoming model IDs.
    for my $md5 (@md5IDs) {
        # Get the list of reactions for this model.
        my @ids = $sap->GetFlat('Aligns', "Aligns(to-link) = ?", [$md5], 'from-link');
        # Store them in the return hash.
        $retVal->{$md5} = \@ids;
    }

    # Return the results.
    return $retVal;
}

#-------------------------------------------------------------------------------
#
#   @metadata = alignment_tree_metadata($sap, $alignID);
#  \@metadata = alignment_tree_metadata($sap, $alignID);
#
# where $sap is the Sapling database object.
#       $alignID is an alignment whose metadata is desired.
# returns (alignment-method, alignment-parameters, alignment-properties,
#          tree-method, tree-parameters, tree-properties)
#
#-------------------------------------------------------------------------------
sub alignment_tree_metadata
{
    # Get the parameters.
    my ($sap, $alignID) = @_;
    my $fields = [];
    if ($alignID) {
        # Retrieve the metadata from the main alignment/tree record.
        ($fields) = $sap->GetAll("AlignmentTree",
                                 "AlignmentTree(id) = ?",
                                 [$alignID],
                                 [qw(alignment-method
                                     alignment-parameters
                                     alignment-properties
                                     tree-method
                                     tree-parameters
                                     tree-properties)]);
    }
    # Return the result in the user-specified manner.
    wantarray ? @$fields : $fields;
}

#-------------------------------------------------------------------------------
#
#    @md5IDs = md5IDs_in_align( $sap, $alignID );
#   \@md5IDs = md5IDs_in_align( $sap, $alignID );
#
# where $sap is the Sapling database object.
#       $alignID is an alignment whose protein list is desired.
#-------------------------------------------------------------------------------
sub md5IDs_in_align
{
    my ( $sap, $alignID ) = @_;
    $sap && $alignID or return wantarray ? () : [];
    my %seen;
    my @md5IDs = grep { ! $seen{$_}++ }
                 $sap->GetFlat('Aligns', 'Aligns(from-link) = ?', [$alignID], 'to-link');
    wantarray ? @md5IDs : \@md5IDs;
}

#-------------------------------------------------------------------------------
#
#   \@seqs               = md5_alignment_by_ID( $sap, $alignID );
# ( \@seqs, \%metadata ) = md5_alignment_by_ID( $sap, $alignID );
#           \%metadata   = md5_alignment_metadata( $sap, $alignID );
#
#       $metadata{ $seqID } = [ $md5, $peg_length, $trim_beg, $trim_end, $location_string ]
#
# where $sap is the Sapling database object.
#       $alignID is an alignment whose MD5 relationship data is desired.
#-------------------------------------------------------------------------------
sub md5_alignment_by_ID
{
    my ( $sap, $alignID ) = @_;
    $sap && $alignID or return ();
    my @align;
    if ( $data_dir ||= locate_data_dir() )
    {
        my $file = "$data_dir/ali$alignID.fa";
        @align = map { $_->[1] = ''; $_ } gjoseqlib::read_fasta( $file ) if -f $file;
    }

    wantarray ? ( \@align, md5_alignment_metadata( $sap, $alignID ) ) : \@align;
}

#-------------------------------------------------------------------------------
#
#   \%md5_row_metadata = md5_alignment_metadata( $sap, $alignID );
#
#        $md5_row_metadata{ $seqID } = [ $md5ID, $peg_length, $trim_beg, $trim_end, $location_string ]
#
# where $TreeServerO is an alignment and tree server object.
#       $alignID is an alignment whose MD5 relationship data is desired.
#
#       $metadata{ $sedID } = [ $md5, $peg_length, $trim_beg, $trim_end, $location_string ]
#-------------------------------------------------------------------------------
sub md5_alignment_metadata
{
    my ( $sap, $alignID ) = @_;
    my %metadata;
    if ( $sap && $alignID )
    {
        %metadata = map { my ($md5, @data) = @$_; ( $md5 => \@data ) }
                    $sap->GetAll('Aligns', 'Aligns(from-link) = ?', [$alignID],
                        [qw(sequence-id to-link len begin end properties)]);
    }
    \%metadata;
}

#-------------------------------------------------------------------------------
#
#    %metadata = alignment_metadata_by_md5( $sap, $alignID, @md5IDs );
#   \%metadata = alignment_metadata_by_md5( $sap, $alignID, @md5IDs );
#
# where  $sap is the Sapling database object.
#        $alignID is an alignment whose MD5 relationship data is desired.
#       \@md5IDs is a list of the md5IDs for which the data are desired.
#
#       $metadata{ $seqID } = [ $md5, $peg_length, $trim_beg, $trim_end, $location_string ]
#-------------------------------------------------------------------------------
sub alignment_metadata_by_md5
{
    my ( $sap, $alignID, @md5IDs ) = @_;
    my %metadata;

    if ( $sap && $alignID && @md5IDs )
    {
        %metadata = map { my ( $seqID, @data ) = @$_; $seqID => \@data }
                    $sap->GetAll('Aligns',
                                 'Aligns(from-link) = ? AND Aligns(to-link) IN ( ' . join( ', ', qw(?) x @md5IDs ) . ' )',
                                 [ $alignID, @md5IDs ],
                                 [qw( sequence-id to-link len begin end properties )]);
    }
    \%metadata;
}

#-------------------------------------------------------------------------------
#
#    @metadata = alignments_metadata_by_md5( $sap, @md5IDs );
#   \@metadata = alignments_metadata_by_md5( $sap, @md5IDs );
#
# where $sap is the Sapling database object.
#       \@md5IDs is a list of the md5IDs for which the data are desired.
#
#       @metadata = ( [ $alignID, $seqID, $md5, $peg_length, $trim_beg, $trim_end, $location_string ], ... )
#-------------------------------------------------------------------------------
sub alignments_metadata_by_md5
{
    my ( $sap, @md5IDs ) = @_;
    my @metadata = ();

    if ( $sap && @md5IDs )
    {
        @metadata = $sap->GetAll('Aligns',
                                 'Aligns(to-link) IN ( ' . join( ', ', qw(?) x @md5IDs ) . ' )',
                                 \@md5IDs,
                                 [qw( from-link sequence-id to-link len begin end properties )]);
    }

    wantarray ? @metadata : \@metadata;
}

#-------------------------------------------------------------------------------
#
#    @treeIDs = all_treeIDs( $sap );
#   \@treeIDs = all_treeIDs( $sap );
#
#-------------------------------------------------------------------------------
sub all_treeIDs
{
    return all_alignIDs(@_);
}

#-------------------------------------------------------------------------------
#
#    @treeIDs = trees_with_md5ID( $sap, $md5 );
#   \@treeIDs = trees_with_md5ID( $sap, $md5 );
#
#-------------------------------------------------------------------------------
sub trees_with_md5ID
{
    return aligns_with_md5ID(@_);
}


#-------------------------------------------------------------------------------
#
#   \%md5_treeIDs = md5IDs_to_trees( $sap, @md5IDs );
#
#-------------------------------------------------------------------------------
sub md5IDs_to_trees
{
    return md5IDs_to_aligns(@_);
}


#-------------------------------------------------------------------------------
#
#    @md5IDs = md5IDs_in_tree( $sap, $treeID );
#   \@md5IDs = md5IDs_in_tree( $sap, $treeID );
#
#-------------------------------------------------------------------------------
sub md5IDs_in_tree
{
    return md5IDs_in_align(@_);
}

#-------------------------------------------------------------------------------
#
#   $tree               = md5_tree_by_ID( $sap, $treeID );
# ( $tree, \%metadata ) = md5_tree_by_ID( $sap, $treeID );
#          \%metadata   = md5_tree_metadata( $sap, $treeID );
#
#       $metadata{ $seqID } = [ $md5, $peg_length, $trim_beg, $trim_end, $location_string ]
#
# where $sap is the Sapling database object.
#       $treeID is a tree whose MD5 relationship data is desired.
#-------------------------------------------------------------------------------
sub md5_tree_by_ID
{
    my ( $sap, $treeID ) = @_;
    my ( $file, $frozen );
    if ( $treeID && ( $data_dir ||= locate_data_dir() ) )
    {
        $file   = "$data_dir/tree$treeID.nwk";
        $frozen = "$data_dir/FROZEN/tree$treeID.frozen";
    }

    my $tree = $frozen && -f $frozen ? thaw_tree_file( $frozen ) :
               $file   && -f $file   ? gjonewicklib::read_newick_tree( $file ) : undef;

    wantarray ? ( $tree, md5_tree_metadata( $sap, $treeID ) ) : $tree;
}

sub thaw_tree_file 
{
    my ( $file ) = @_;
    return undef unless $file && -s $file;

    my $ref_frozen_tree = gjoseqlib::slurp( $file );
    my $tree = thaw($$ref_frozen_tree);

    return $tree;
}


sub md5_tree_metadata
{
    return md5_alignment_metadata(@_);
}

#-------------------------------------------------------------------------------
#
#    @fids = md5_to_pegs($sap, $md5);
#   \@fids = md5_to_pegs($sap, $md5);
#
# where $sap is the Sapling database object.
#       $md5 is an MD5 protein ID
#
# This method returns all the pegs that produce the indicated protein.
#
#-------------------------------------------------------------------------------
sub md5_to_pegs {
    # Get the parameters.
    my ($sap, $md5) = @_;
    # Read the features for the indicated protein.
    my @retVal = $sap->GetFlat("IsProteinFor", 'IsProteinFor(from-link) = ?', 
                            [$md5], 'to-link');
    # Return the result in the caller-specified manner.
    wantarray ? @retVal : \@retVal;
}


#-------------------------------------------------------------------------------
#
#   $md5 = peg_to_md5($sap, $peg);
#
# Return the MD5 ID of the protein produced by the indicated peg.
#
# where $sap is the Sapling database object.
#       $peg is a feature ID.
#
#-------------------------------------------------------------------------------
sub peg_to_md5 {
    # Get the parameters.
    my ($sap, $peg) = @_;
    # Read the protein for the indicated PEG.
    my ($retVal) = $sap->GetFlat("Produces", "Produces(from-link) = ?",
                                 [$peg], 'to-link');
    # Return the result.
    $retVal;
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
    my %fid_id_to_md5_id_map;

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
        $fid_id_to_md5_id_map{$id} = $md5ID;
    }

    return (\%md5_metadata, \%fid_id_to_md5_id_map);
}

sub fid_align_to_md5_align
{
    my ( $fid_align, $fid_id_to_md5_id_map ) = @_;
    $fid_align && ref( $fid_align ) eq 'ARRAY' &&
        $fid_id_to_md5_id_map && ref( $fid_id_to_md5_id_map ) eq 'HASH'
        or return ();

    my @md5_align;

    foreach ( @$fid_align )
    {
        my $id = $_->[0];
        my $md5ID = $fid_id_to_md5_id_map->{$id};
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
    my %md5_id_to_fid_ids_map;

    foreach my $md5ID ( keys %$md5_metadata )
    {
        my $md5Metadata = $md5_metadata->{$md5ID};
        my ($md5, $len, $beg, $end, $location) = @$md5Metadata;
        my @fids = md5_to_pegs( $sap, $md5 );
        @fids = ( $md5 ) if ! @fids && $relaxed;
        foreach my $fid ( @fids )
        {
            my $fidID = $fid;
            if ($fids_seen{$fid}++) {
                $fidID = "$fid-" . $fids_seen{$fid};
            }
            $fid_metadata{$fidID} = [$fid, $len, $beg, $end, $location];
            push @{$md5_id_to_fid_ids_map{$md5ID}}, $fidID;
        }
    }

    return (\%fid_metadata, \%md5_id_to_fid_ids_map);
}

sub md5_align_to_fid_align
{
    my ( $md5_align, $md5_id_to_fid_ids_map ) = @_;
    $md5_align && ref( $md5_align ) eq 'ARRAY' && $md5_id_to_fid_ids_map &&
        ref( $md5_id_to_fid_ids_map ) eq 'HASH'
        or return ();

    my @fid_align;
    my %fid_metadata;

    foreach ( @$md5_align )
    {
        my $md5ID  = $_->[0];
        my @fidIDs = @{$md5_id_to_fid_ids_map->{$md5ID}};
        foreach my $fidID ( @fidIDs )
        {
            push @fid_align, [ $fidID, $_->[1], $_->[2] ];
        }
    }

    return \@fid_align;
}


sub fid_tree_to_md5_tree
{
    my ( $fid_tree, $fid_id_to_md5_id_map ) = @_;
    $fid_tree && ref( $fid_tree ) eq 'ARRAY' &&
        $fid_id_to_md5_id_map && ref( $fid_id_to_md5_id_map ) eq 'HASH'
        or return undef;

    gjonewicklib::newick_relabel_tips( gjonewicklib::newick_subtree( $fid_tree, keys %$fid_id_to_md5_id_map ), $fid_id_to_md5_id_map );
}


sub md5_tree_to_fid_tree
{
    my ( $md5_tree, $md5_id_to_fid_ids_map ) = @_;
    $md5_tree && ref( $md5_tree ) eq 'ARRAY' &&
        $md5_id_to_fid_ids_map && ref( $md5_id_to_fid_ids_map ) eq 'HASH'
        or return ();

    my @tips = gjonewicklib::newick_tip_list( $md5_tree );
    @tips or return undef;

    my $prune = 0;
    foreach my $md5ID ( @tips )
    {
        $prune = 1 if (! $md5_id_to_fid_ids_map->{$md5ID});
    }

    $md5_tree = gjonewicklib::newick_subtree( $md5_tree, [ keys %$md5_id_to_fid_ids_map ] ) if $prune;
    return expand_duplicate_tips( gjonewicklib::copy_newick_tree( $md5_tree ), $md5_id_to_fid_ids_map );
}


sub md5_align_and_tree_to_fid_version {
    my ($sap, $md5_align, $md5_tree, $md5_metadata, $relaxed) = @_;
    my ($fid_metadata, $md5_id_to_fid_ids_map) = map_md5_to_fid($sap, $md5_metadata, $relaxed);
    my $fid_tree = md5_tree_to_fid_tree($md5_tree, $md5_id_to_fid_ids_map);
    my $fid_align = md5_align_to_fid_align($md5_align, $md5_id_to_fid_ids_map);
    return ($fid_align, $fid_tree, $fid_metadata);
}


sub fid_align_and_tree_to_md5_version {
    my ($sap, $fid_align, $fid_tree, $fid_metadata, $relaxed) = @_;
    my ($md5_metadata, $fid_id_to_md5_id_map) = map_fid_to_md5($sap, $fid_metadata, $relaxed);
    my $md5_tree = fid_tree_to_md5_tree($fid_tree, $fid_id_to_md5_id_map);
    my $md5_align = fid_align_to_md5_align($fid_align, $fid_id_to_md5_id_map);
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
sub overlaps {
    # Get the parameters.
    my ($beg1, $end1, $beg2, $end2) = @_;
    # Compute the number of overlapping residues.
    my $over = min($end1, $end2) - max($beg1, $beg2) + 1;
    # Return TRUE if the overlap is 80% of the shorter length.
    return $over >= 0.8 * min($end1 - $beg1 + 1, $end2 - $beg2 + 1);
}

sub min { $_[0] < $_[1] ? $_[0] : $_[1] }
sub max { $_[0] > $_[1] ? $_[0] : $_[1] }


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


1;
