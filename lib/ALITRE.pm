#
# Copyright (c) 2003-2006 University of Chicago and Fellowship
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
package ALITRE;

use strict;
use AlignsAndTrees;
use Tracer;
use SeedUtils;
use ServerThing;

=head1 Alignment and Tree Server Function Object

This file contains the functions and utilities used by the Alignment and Tree Server
(B<alitre_server.cgi>). The various methods listed in the sections below represent
function calls direct to the server. These all have a signature similar to the
following.

    my $results = $altObject->function_name($args);

where C<$altObject> is an object created by this module, 
C<$args> is a parameter structure, and C<function_name> is the server function name.
The output $results is a scalar, generally a hash reference, but sometimes a
string or a list reference.

=head2 Constructor

Use

    my $altObject = ALITREserver->new();

to create a new alignment/tree server function object. The server function object
is used to invoke the L</Primary Methods> listed below. See L<ALITREserver> for
more information on how to create this object and the options available.

=cut

#
# Actually, if you are using ALITRE.pm, you should do ALITRE->new(), not ALITREserver->new()
# That comment above is for the benefit of the pod doc stuff on how to use ALITREserver
# that is generated from this file.
#

sub new {
 my ( $class, $sap ) = @_;

 # Create the sapling object.
 if ( !defined $sap ) {
  $sap = ERDB::GetDatabase('Sapling');
 }

 # Create the server object.
 my $retVal = { db => $sap };

 # Bless and return it.
 bless $retVal, $class;
 return $retVal;
}

=head1 Primary Methods

=head2 Server Utility Methods

You will not use the methods in this section very often. Currently, the only one
present (L</methods>) is used by the server framework for maintenance and control
purposes.

=head3 methods

    my $methodList =        $altObject->methods();

Return a reference to a list of the methods allowed on this object.

=cut

use constant METHODS => [
 qw(
   alignment_metadata_by_md5
   alignments_metadata_by_md5
   alignment_tree_metadata
   aligns_with_md5ID
   all_alignIDs
   all_treeIDs
   expand_duplicate_tips
   fid_align_and_tree_to_md5_version
   fid_align_to_md5_align
   fid_tree_to_md5_tree
   get_projections
   map_fid_to_md5
   map_md5_to_fid
   md5IDs_in_align
   md5IDs_in_tree
   md5_align_and_tree_by_ID
   md5_align_and_tree_to_fid_version
   md5_align_to_fid_align
   md5_alignment_by_ID
   md5_alignment_metadata
   md5_tree_by_ID
   md5_tree_to_fid_tree
   trees_with_md5ID
   )
];

sub methods {

 # Get the parameters.
 my ($self) = @_;

 # Return the result.
 return METHODS;
}

=head2 Client Methods

=head3 alignment_tree_metadata

    my $alignHash =         $altObject->alignment_tree_metadata({
                                -ids => [$alt1, $alt2, ...]
                            });

Return the construction metadata for the alignment and tree in each specified
alignment/tree pair. The construction metadata describes how the alignment
and tree were built from the raw data.

=over 4

=item parameters

The parameter should be a reference to a hash with the following keys:

=over 8

=item -ids

Reference to a list of alignment/tree IDs.

=back

=item RETURN

Returns a reference to a hash mapping each incoming alignment ID to a 6-tuple
of metadata information, including (0) the name of the method used to build
the alignment, (1) the parameters passed to the alignment method, (2) the
properties of the alignment process, (3) the name of the method used to build
the tree, (4) the parameters passed to the tree method, and (5) the properties
of the tree-building process.

    $alignHash => { $alt1 => [$almethod1, $alparms1, $alprops1,
                              $trmethod1, $trparms1, $trprops1],
                    $alt2 => [$almethod2, $alparms2, $alprops2,
                              $trmethod2, $trparms2, $trprops2],
                    ... };

=back

=cut

sub alignment_tree_metadata {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the Sapling database.
 my $sap = $self->{db};

 # Create the return hash.
 my $retVal = {};

 # Get the alignment IDs.
 my $ids = ServerThing::GetIdList( -ids => $args );

 # Loop through the IDs, extracting the metadata.
 for my $id (@$ids) {
  $retVal->{$id} = AlignsAndTrees::alignment_tree_metadata( $sap, $id );
 }

 # Return the result hash.
 return $retVal;
}

=head3 aligns_with_md5ID

    my $protHash =          $altObject->aligns_with_md5ID({
                                -ids => [$prot1, $prot2, ...]
                            });

Return a list of the alignment/tree pairs containing each of the specified proteins.

=over 4

=item parameters

The parameter should be a reference to a hash with the following keys:

=over 8

=item -ids

Reference to a list of MD5 protein IDs.

=back

=item RETURN

Returns a reference to a hash mapping each incoming protein ID to a list of the
IDs for the alignments containing that protein.

    $protHash = { $prot1 => [$alt1a, $alt1b, ...],
                  $prot2 => [$alt2a, $alt2b, ...],
                  ... };

=back

=cut

sub aligns_with_md5ID {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the Sapling database.
 my $sap = $self->{db};

 # Create the return hash.
 my $retVal = {};

 # Get the list of incoming IDs.
 my $ids = ServerThing::GetIdList( -ids => $args );

 # Loop through the protein IDs, finding the alignments.
 for my $id (@$ids) {
  $retVal->{$id} = AlignsAndTrees::aligns_with_md5ID( $sap, $id );
 }

 # Return the result hash.
 return $retVal;
}

=head3 all_alignIDs

    my $idList =            $altObject->all_alignIDs();

Return a list of all the alignment IDs in the database.

=over 4

=item RETURN

Returns a reference to a list of alignment IDs for all the alignments in the database.

    $idList = [$alt1, $alt2, ...];

=back

=cut

sub all_alignIDs {

 # Get the parameters.
 my ($self) = @_;

 # Get the sapling database.
 my $sap = $self->{db};

 # Get the list of IDs.
 my $retVal = AlignsAndTrees::all_alignIDs($sap);

 # Return it to the caller.
 return $retVal;
}

=head3 all_treeIDs

    my $idList =            $altObject->all_treeIDs();

Return a list of all the tree IDs in the database.

=over 4

=item RETURN

Returns a reference to a list of IDs for all the trees in the database. (Note:
as currently construed, this is the same as a list of all the alignment IDs, since
each alignment has exactly one associated tree and it has the same ID.)

    $idList = [$alt1, $alt2, ...];

=back

=cut

sub all_treeIDs {
 return all_alignIDs(@_);
}

=head3 expand_duplicate_tips

    my $newTree =       $altObject->expand_duplicate_tips({
                            -tree => $actualTree,
                            -map => { $oldName1 => [$newName1a, $newName1b, ...],
                                      $oldName2 => [$newName2a, $newName2b, ...],
                                      ... }
                        });

Rename and possibly expand the tips of the specified tree data structure using the
specified mapping.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -tree

Reference to a list that encodes a newick phylogenetic tree.

=item -map

Reference to a hash that maps node names to new node names. Each new node name is
a reference to a list of names. If the list is a singleton, the mapping is a simple
renaming. If the list contains muliple entries, the node will be expanded into duplicates.

=back

=item RETURN

Returns a new version of the tree with the specified renamings performed.

=back

=cut

sub expand_duplicate_tips {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the incoming tree.
 my $tree = $args->{-tree};
 if ( !$tree ) {
  Confess("Missing -tree parameter to expand_duplicate_tips.");
 } elsif ( ref $tree ne 'ARRAY' ) {
  Confess("Invalid -tree parameter for expand_duplicate_tips.");
 }

 # Get the incoming name map.
 my $map = $args->{-map};
 if ( !$map ) {
  Confess("Missing -map parameter to expand_duplicate_tips.");
 } elsif ( ref $map ne 'HASH' ) {
  Confess("Invalid -map parameter to expand_duplicate_tips.");
 }

 # Perform the expansion. Note that the expansion is actually done in place,
 # and the method simply returns the incoming argument.
 my $retVal = AlignsAndTrees::expand_duplicate_tips( $tree, $map );

 # Return the modified tree.
 return $retVal;
}

=head3 fid_align_and_tree_to_md5_version

    my $md5Tuple =          $altObject->fid_align_and_tree_to_md5_version({
                                -align => $fid_align,
                                -tree => $fid_tree,
                                -meta => $fid_meta,
                                -relaxed => 1
                            });

Convert a PEG-based alignment/tree pair to an MD5-based alignment/tree pair. Each
PEG identifier in the alignment and tree will be converted to the corresponding MD5
protein identifier. This may cause some nodes in the tree and items in the alignment
to be collapsed into a single instance, since multiple PEGs can produce the same
protein.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys:

=over 8

=item -align

The PEG-based alignment to convert.

=item -tree

The corresponding phylogenetic tree.

=item -meta

Reference to a hash mapping each feature ID in the alignment and tree to a
description of which part of the resulting protein was used.

=item -relaxed (optional)

If TRUE, then incoming feature IDs that are not found in the database will be
left untranslated in the output. Otherwise, such IDs will cause an error. The
default is FALSE.

=back

=item RETURN

Returns a reference to a 3-tuple containing (0) the MD5 version of the incoming
alignment, (1) the MD5 version of the incoming tree, and (2) a reference to a
hash describing which portion of each protein was used in the alignment.

    $md5Tuple = [$md5_align, $md5_tree, $md5_metadata];

=back

=cut

sub fid_align_and_tree_to_md5_version {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the Sapling database.
 my $sap = $self->{db};

 # Get the parameters.
 my $align = $args->{-align}
   || Confess("No alignment specified in fid_align_and_tree_to_fid_version.");
 my $tree = $args->{-tree}
   || Confess("No tree specified in fid_align_and_tree_to_fid_version.");
 my $meta = $args->{-meta}
   || Confess("No metadata specified in fid_align_and_tree_to_fid_version.");
 my $relaxed = $args->{-relaxed} || 0;

 # Convert the alignment and tree.
 my ( $newAlign, $newTree, $newMeta ) =
   AlignsAndTrees::fid_align_and_tree_to_md5_version( $sap, $align, $tree,
  $meta, $relaxed );

 # Return the results.
 return [ $newAlign, $newTree, $newMeta ];

}

=head3 fid_align_to_md5_align

    my $md5align =          $altObject->fid_align_to_md5_align({
                                -align => $fid_align,
                                -map => $fid_to_md5_map
                            });

Use a map produced by L</map_fid_to_md5> to convert a PEG-based tree to an MD5-based tree.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys:

=over 8

=item -align

The PEG-based alignment to be converted.

=item -map

A hash tha maps each feature ID in the alignment to the corresponding MD5 protein
ID.

=back

=item RETURN

Returns a new version of the alignment with the feature IDs replaced by MD5 protein
IDs using the data in the map.

=back

=cut

sub fid_align_to_md5_align {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the incoming alignment and map.
 my $align = $args->{-align}
   || Confess("No alignment specified in fid_align_to_md5_align.");
 my $map = $args->{-map}
   || Confess("No map specified in fid_align_to_md5_align,");

 # Perform the conversion.
 my $retVal = AlignsAndTrees::fid_align_to_md5_align( $align, $map );

 # Return the result.
 return $retVal;
}

=head3 fid_tree_to_md5_tree

    my $md5tree =          $altObject->fid_tree_to_md5_tree({
                                -tree => $fid_tree,
                                -map => $fid_to_md5_map
                            });

Use a map produced by L</map_fid_to_md5> to convert a PEG-based tree to an MD5-based tree.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys:

=over 8

=item -tree

The PEG-based tree to be converted.

=item -map

A hash tha maps each feature ID in the tree to the corresponding MD5 protein
ID.

=back

=item RETURN

Returns a new version of the tree with the feature IDs replaced by MD5 protein
IDs using the data in the map.

=back

=cut

sub fid_tree_to_md5_tree {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the incoming tree and map.
 my $tree = $args->{-tree}
   || Confess("No tree specified in fid_tree_to_md5_tree.");
 my $map = $args->{-map}
   || Confess("No map specified in fid_tree_to_md5_tree,");

 # Perform the conversion.
 my $retVal = AlignsAndTrees::fid_tree_to_md5_tree( $tree, $map );

 # Return the result.
 return $retVal;
}

=head3 get_projections

    my $protHash =          $altObject->get_projections({
                                -ids => [$prot1, $prot2, ...],
                                -minScore => 0,
                                -details => 1
                            });

Get all the proteins that are clear bidirectional best hits (projections) of the 
specified proteins. The call can specify whether or not to include the details of 
the projection (score, context, percent identity) and also whether to only include
projections with a certain minimum score.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of MD5 protein IDs.

=item -minScore (optional)

Minimum score for projections to be included. Only projections with the specified
score or greater will be included in the output. The score ranges from 0 to 1. The
default is C<0>, which includes everything.

=item -details (optional)

If TRUE, then for each projected protein, the gene-context count, percent identity,
and score will be included in the output. The default is FALSE (only the protein IDs
are returned).

=back

=item RETURN

Returns a reference to a hash mapping each incoming protein ID to a list of the
proteins that are clear bidirectional best hits. If C<-details> is FALSE, each
protein will be represented by its MD5 protein ID. If C<-details> is TRUE, each
protein will be represented by a 4-tuple consisting of (0) the MD5 protein ID,
(1) the number of homologous genes in the immediate context of the two proteins
(up to a maximum of 10), (2) the percent match between the two protein sequences,
and (3) the score of the projections (0 to 1).

=over 8

=item -details FALSE

    $protHash = { $prot1 => [$prot1a, $prot1b, ...],
                  $prot2 => [$prot2a, $prot2b, ...],
                  ... };

=item -details TRUE

    $protHash = { $prot1 => [[$prot1a, $context1a, $match1a, $score1a],
                             [$prot1b, $context1b, $match1b, $score1b], 
                             ... ],
                  $prot2 => [[$prot2a, $context2a, $match2a, $score2a],
                             [$prot2b, $context2b, $match2b, $score2b], 
                             ... ],
                  ... };

=back

=back

=cut

sub get_projections {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of protein IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Compute the minimum score.
    my $minScore = $args->{-minScore} || 0;
    # Compute the detail level.
    my $details = ($args->{-details} ? 1 : 0);
    # Create the return variable.
    my $retVal = {};
    # Loop through the proteins.
    for my $id (@$ids) {
        # Compute the filter parameters.
        my $parms = [$id, $minScore];
        # We'll put our query results in here. The query itself depends on the
        # details mode.
        my @projections;
        # The query has to be tried twice: once from each direction.
        for my $rel (qw(ProjectsOnto IsProjectedOnto)) {
            if ($details) {
                # Here we have a detailed projection request.
                push @projections, $sap->GetAll($rel,
                    "$rel(from-link) = ? AND $rel(score) >= ?",
                    $parms, "to-link gene-context percent-identity score");
            } else {
                # Here we want only the projected protein ID.
                push @projections, $sap->GetFlat($rel,
                    "$rel(from-link) = ? AND $rel(score) >= ?", $parms, "to-link");
            }
        }
        # Store the results found in the return hash.
        $retVal->{$id} = \@projections;
    }
    # Return the hash of results.
    return $retVal;
}

=head3 map_fid_to_md5

    my $md5Tuple =          $altObject->map_fid_to_md5({
                                -meta => $fid_metadata,
                                -relaxed => 0
                            });

Analyze the metadata for a PEG-basedalignment/tree pair and compute the metadata for
the corresponding MD5-based data structions along with a mapping from the PEG IDs
to MD5 IDs.

=over 4

=item parameter

Reference to a hash with the following keys:

=over 8

=item -meta

Reference to a hash mapping each FIG feature ID in an alignment/tree pair to
information describing which part of each feature's protein was used.

=item -relaxed (optional)

If TRUE, then incoming feature IDs that are not found in the database will be
left untranslated in the output. Otherwise, such IDs will cause an error. The
default is FALSE.

=back

=item RETURN

Returns a reference to a 2-tuple containing (0) the MD5-based metadata hash
creating from the incoming hash and (1) a hash mapping each incoming feature ID
to the corresponding MD5 protein ID.

    $md5Tuple => [$md5_metadata, { $fida => $md5a, $fidb => $md5b, ... }];

=back

=cut

sub map_fid_to_md5 {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the Sapling database.
 my $sap = $self->{db};

 # Get the metadata structure.
 my $meta = $args->{-meta}
   || Confess("No metadata structure passed to map_fid_to_md5.");

 # Compute the relax flag.
 my $relaxed = $args->{-relaxed} || 0;

 # Perform the conversion.
 my ( $newMeta, $map ) =
   AlignsAndTrees::map_fid_to_md5( $sap, $meta, $relaxed );

 # Return the result.
 return [ $newMeta, $map ];
}

=head3 map_md5_to_fid

    my $fidTuple =          $altObject->map_md5_to_fid({
                                -meta => $md5_metadata,
                                -relaxed => 0
                            });

Analyze the metadata for an MD5 alignment/tree pair and compute the metadata for
the corresponding PEG-based data structions along with a mapping from the MD5 IDs
to the PEG IDs.

=over 4

=item parameter

Reference to a hash with the following keys:

=over 8

=item -meta

Reference to a hash mapping each MD5 protein ID in an alignment/tree pair to
information describing which part of each protein was used.

=item -relaxed (optional)

If TRUE, then incoming MD5 IDs that are not found in the database will be
left untranslated in the output. Otherwise, such IDs will cause an error. The
default is FALSE.

=back

=item RETURN

Returns a reference to a 2-tuple containing (0) the PEG-based metadata hash
created from the incoming hash and (1) a hash mapping each incoming MD5 protein
ID to a list of corresponding FIG feature IDs.

    $fidTuple => [$fid_metadata, { $md5a => [$fida1, $fida2, ...],
                                   $md5b => [$fidb1, $fidb2, ...],
                                   ... }];

=back

=cut

sub map_md5_to_fid {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the Sapling database.
 my $sap = $self->{db};

 # Get the metadata structure.
 my $meta = $args->{-meta}
   || Confess("No metadata structure passed to map_md5_to_fid.");

 # Compute the relax flag.
 my $relaxed = $args->{-relaxed} || 0;

 # Perform the conversion.
 my ( $newMeta, $map ) =
   AlignsAndTrees::map_md5_to_fid( $sap, $meta, $relaxed );

 # Return the result.
 return [ $newMeta, $map ];
}

=head3 md5IDs_in_align

    my $altHash =           $altObject->md5IDs_in_align({
                                -ids => [$alt1, $alt2, ...]
                            });

For each incoming alignment ID, return a list of the MD5 protein IDs for the proteins
found in the alignment.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys:

=over 8

=item -ids

Reference to a list of alignment IDs.

=back

=item RETURN

Returns a reference to a hash mapping each incoming alignment ID to a list of
the proteins found in the alignment. Each protein is represented by an MD5 protein
ID.

    $altHash = { $alta => [$md5a1, $md5a2, ... ],
                 $altb => [$md5b1, $md5b2, ... ],
                 ... };

=back

=cut

sub md5IDs_in_align {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the list of alignment IDs.
 my $ids = ServerThing::GetIdList( -ids => $args );

 # Get the sapling database.
 my $sap = $self->{db};

 # Declare the return hash.
 my $retVal = {};

 # Loop through the incoming IDs.
 for my $id (@$ids) {

  # Get the MD5s for this alignment.
  my $md5List = AlignsAndTrees::md5IDs_in_align( $sap, $id );

  # Store them in the return hash.
  $retVal->{$id} = $md5List;
 }

 # Return the result hash.
 return $retVal;
}

=head3 md5IDs_in_tree

    my $altHash =           $altObject->md5IDs_in_tree({
                                -ids => [$alt1, $alt2, ...]
                            });

For each incoming tree ID, return a list of the MD5 protein IDs for the proteins
found in the tree.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys:

=over 8

=item -ids

Reference to a list of tree IDs.

=back

=item RETURN

Returns a reference to a hash mapping each incoming tree ID to a list of
the proteins found in the tree. Each protein is represented by an MD5 protein
ID.

    $altHash = { $alta => [$md5a1, $md5a2, ... ],
                 $altb => [$md5b1, $md5b2, ... ],
                 ... };

=back

=cut

sub md5IDs_in_tree {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the list of tree IDs.
 my $ids = ServerThing::GetIdList( -ids => $args );

 # Get the sapling database.
 my $sap = $self->{db};

 # Declare the return hash.
 my $retVal = {};

 # Loop through the incoming IDs.
 for my $id (@$ids) {

  # Get the MD5s for this tree.
  my $md5List = AlignsAndTrees::md5IDs_in_tree( $sap, $id );

  # Store them in the return hash.
  $retVal->{$id} = $md5List;
 }

 # Return the result hash.
 return $retVal;
}

=head3 md5_align_and_tree_by_ID

    my $tupleHash =         $altObject->md5_align_and_tree_by_ID({
                                -ids => [$alt1, $alt2, ...]
                            });

Return the alignment and tree for each specified ID. The return hash will contain
a 3-tuple for each tree ID consisting of the alignment, the tree, and the metadata
describing the proteins involved.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys:

=over 8

=item -ids

Reference to a list of alignment/tree pair IDs.

=back

=item RETURN

Returns a reference to a hash mapping each incoming ID to a 3-tuple containing (0) the
identified MD5 protein alignment, (1) the associated phylogenetic tree, and (2) a
hash describing what portion of each protein was used in the alignment.

    $tupleHash = { $alt1 => [$md5_align1, $md5_tree1, $md5_metadata1],
                   $alt2 => [$md5_align2, $md5_tree2, $md5_metadata2],
                   ... };

=back

=cut

sub md5_align_and_tree_by_ID {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the Sapling database.
 my $sap = $self->{db};

 # Declare the return hash.
 my $retVal = {};

 # Get the list of incoming IDs.
 my $ids = ServerThing::GetIdList( -ids => $args );

 # Loop through the list.
 for my $id (@$ids) {

  # Get the alignment and the metadata.
  my ( $align, $meta ) = AlignsAndTrees::md5_alignment_by_ID( $sap, $id );

  # Get the tree as well.
  my $tree = AlignsAndTrees::md5_tree_by_ID( $sap, $id );

  # Return all three items.
  $retVal->{$id} = [ $align, $tree, $meta ];
 }

 # Return the result hash.
 return $retVal;
}

=head3 md5_align_and_tree_to_fid_version

    my $fidTuple =          $altObject->md5_align_and_tree_to_fid_version({
                                -align => $md5_align,
                                -tree => $md5_tree,
                                -meta => $md5_metadata,
                                -relaxed => 1
                            });

Convert an MD5 alignment/tree pair to a PEG-based alignment-tree pair. Each protein in
the alignment or tree will be translated to a corresponding FIG feature ID. In some
cases, this may cause a single protein to be replicated to include all the features
that produce that protein.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -align

Reference to the MD5 alignment to be converted.

=item -tree

Reference to the corresponding phylogenetic tree.

=item -meta

Reference to a hash mapping each MD5 protein ID in the alignment and tree to a
description of what section of the protein was used.

=item -relaxed (optional)

If TRUE, then incoming feature IDs that are not found in the database will be
left untranslated in the output. Otherwise, such IDs will cause an error. The
default is FALSE.

=back

=item RETURN

Returns a reference to a 3-tuple containing (0) a PEG-based version of the
incoming alignment, (1) a PEG-based version of the incoming tree, and (2) a
reference to a hash mapping each feature ID in the new alignment and tree to
a description of what section of the feature's protein was used.

    $fidTuple = [$fid_align, $fid_tree, $fid_metadata];

=back

=cut

sub md5_align_and_tree_to_fid_version {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the sapling database.
 my $sap = $self->{db};

 # Get the parameters.
 my $align = $args->{-align}
   || Confess("No alignment specified in md5_align_and_tree_to_fid_version.");
 my $tree = $args->{-tree}
   || Confess("No tree specified in md5_align_and_tree_to_fid_version.");
 my $meta = $args->{-meta}
   || Confess("No metadata specified in md5_align_and_tree_to_fid_version.");
 my $relaxed = $args->{-relaxed} || 0;

 # Convert the alignment and tree.
 my ( $newAlign, $newTree, $newMeta ) =
   AlignsAndTrees::md5_align_and_tree_to_fid_version( $sap, $align, $tree,
  $meta, $relaxed );

 # Return the results.
 return [ $newAlign, $newTree, $newMeta ];
}

=head md5_align_to_fid_align

    my $fidAlign =          $altObject->md5_align_to_fid_align({
                                -align => $md5_align,
                                -map => $md5_to_fid_map
                            });

Use a map produced by L</map_fid_to_md5> to convert an MD5-based alignment to a PEG-based
alignment. Since a single protein may be generated by multiple features, this could
result in alignment entries being replicated in the result.

=over 4

=item parameter

The parameter should be a reference to a hash containing the following keys.

=over 8

=item -align

The MD5-based alignment to be converted.

=item -map

Reference to a hash mapping each MD5 protein ID to a list of the corresponding FIG
feature IDs.

=back

=item RETURN

Returns a new version of the alignment with the MD5 protein IDs replaced by FIG
feature IDs.

=back

=cut

sub md5_align_to_fid_align {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the alignment and the map.
 my $align = $args->{-align}
   || Confess("No alignment specified in md5_align_to_fid_align.");
 my $map = $args->{-map}
   || Confess("No map specified in md5_align_to_fid_align.");

 # Perform the conversion.
 my $retVal = AlignsAndTrees::md5_align_to_fid_align( $align, $map );

 # Return the result.
 return $retVal;
}

=head3 md5_alignment_by_ID

    my $altHash =           $altObject->md5_alignment_by_ID({
                                -ids => [$alt1, $alt2, ...]
                            });

Return the alignments with the specified IDs. The return hash will contain a
2-tuple for each alignment ID consisting of the alignment itself followed by
the metadata describing the proteins in the alignment.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys:

=over 8

=item -ids

Reference to a list of alignment IDs.

=back

=item RETURN

Returns a reference to a hash mapping each incoming ID to an MD5 alignment.

    $altHash = { $alt1 => $md5_align1, $alt2 => $md5_align2, ... };

=back

=cut

sub md5_alignment_by_ID {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the Sapling database.
 my $sap = $self->{db};

 # Get the list of incoming IDs.
 my $ids = ServerThing::GetIdList( -ids => $args );

 # Declare the return hash.
 my $retVal = {};

 # Loop through the incoming IDs.
 for my $id (@$ids) {

  # Get the tree and metadata for this ID.
  $retVal->{$id} = [ AlignsAndTrees::md5_alignment_by_ID( $sap, $id ) ];
 }

 # Return the result hash.
 return $retVal;

}


=head3 md5_alignment_metadata

    my $alignHash =         $altObject->md5_alignment_metadata({
                                -ids => [$alignID1, alignID2, ...]
                            });

Return the sequence metadata for the alignment.

=over 4

=item parameters

The parameter should be a reference to a hash with the following keys:

=over 8

=item -ids

Reference to a list of alignment/tree IDs.

=back

=item RETURN

Returns a reference to a hash mapping each incoming alignment ID to a 5-tuple
of metadata information.

    $metadataHash => { $alignID1 => [ $md5, $peg_length, $trim_beg, $trim_end, $location_string ], 
                       $alignID2 => [ $md5, $peg_length, $trim_beg, $trim_end, $location_string ], 
                       ... };

=back

=cut

sub md5_alignment_metadata {
    my ($self, $args) = @_;
    my $sap = $self->{db};
    my $ids = ServerThing::GetIdList( -ids => $args );
    my $retVal = {};
    
    for my $id (@$ids) {
        $retVal->{$id} = AlignsAndTrees::md5_alignment_metadata( $sap, $id );
    }

    return $retVal;
}


=head3 alignment_metadata_by_md5

    my $alignHash = $altObject->alignment_metadata_by_md5({
                                -ids => [$alignID, md5ID1, md5ID2, ...]
                            });

Return the sequence metadata for the alignment.

=over 4

=item parameters

The parameter should be a reference to a hash with the following keys:

=over 8

=item -ids

Reference to a list of alignment whose first element is an align ID followed by MD5 IDs.

=back

=item RETURN

Returns a reference to a hash mapping each alignmen row ID to a 5-tuple
of metadata information.

    $metadataHash => { $seqID1 => [ $md5, $peg_length, $trim_beg, $trim_end, $location_string ], 
                       $seqID2 => [ $md5, $peg_length, $trim_beg, $trim_end, $location_string ], 
                       ... };

=back

=cut

sub alignment_metadata_by_md5 {
    my ($self, $args) = @_;
    my $sap = $self->{db};
    my $ids = ServerThing::GetIdList( -ids => $args );

    my ($alignID, @md5IDs) = @$ids;

    my $metadata = AlignsAndTrees::alignment_metadata_by_md5( $sap, $alignID, @md5IDs );

    return $metadata;
}


=head3 alignments_metadata_by_md5

    my $metadataHash = $altObject->alignments_metadata_by_md5({
                                -ids => [md5ID1, md5ID2, ...]
                            });

Return the sequence metadata for the alignment.

=over 4

=item parameters

The parameter should be a reference to a hash with the following keys:

=over 8

=item -ids

Reference to a list of MD5 IDs.

=back

=item RETURN

Returns a reference to a list of a 7-tuple of metadata information.

    $metadataList => [ [ $alignID1, $seqID1, $md5, $peg_length, $trim_beg, $trim_end, $location_string ], 
                       [ $alignID2, $seqID2, $md5, $peg_length, $trim_beg, $trim_end, $location_string ], 
                       ... ] ];

=back

=cut

sub alignments_metadata_by_md5 {
    my ($self, $args) = @_;
    my $sap = $self->{db};
    my $ids = ServerThing::GetIdList( -ids => $args );
    
    my @md5IDs = @$ids;

    my $metadata = AlignsAndTrees::alignments_metadata_by_md5( $sap, @md5IDs );
    
    return $metadata;
}


=head3 md5_tree_by_ID

    my $tupleHash =         $altObject->md5_tree_by_ID({
                                -ids => [$alt1, $alt2, ...]
                            });

Return the trees with the specified IDs. The return hash will contain a 2-tuple
for each tree ID consisting of the tree itself followed by the metadata describing
the proteins in the tree.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys:

=over 8

=item -ids

Reference to a list of tree IDs.

=back

=item RETURN

Returns a reference to a hash that maps each incoming tree ID to a 2-tuple consisting of
(0) a data structure containing the identified phylogenetic tree represented as a
newick-format list, and (1) a hash containing the metadata for the leaves of the tree.

    $tupleHash = { $tree1 => [$md5_tree1, $md5_metadata1],
                   $tree2 => [$md5_tree2, $md5_metadata2],
                   ... };

=back

=cut

sub md5_tree_by_ID {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the sapling database.
 my $sap = $self->{db};

 # Get the incoming IDs.
 my $ids = ServerThing::GetIdList( -ids => $args );

 # Declare the return hash.
 my $retVal = {};

 # Loop through the incoming IDs.
 for my $id (@$ids) {

  # Get the tree and metadata for this ID.
  $retVal->{$id} = [ AlignsAndTrees::md5_tree_by_ID( $sap, $id ) ];
 }

 # Return the result hash.
 return $retVal;
}

=head md5_tree_to_fid_tree

    my $fidtree =          $altObject->md5_tree_to_fid_tree({
                                -tree => $md5_tree,
                                -map => $md5_to_fid_map
                            });

Use a map produced by L</map_fid_to_md5> to convert an MD5-based tree to a PEG-based
tree. Since a single protein may be generated by multiple features, this could
result in tree nodes being replicated in the result.

=over 4

=item parameter

The parameter should be a reference to a hash containing the following keys.

=over 8

=item -tree

The MD5-based tree to be converted.

=item -map

Reference to a hash mapping each MD5 protein ID to a list of the corresponding FIG
feature IDs.

=back

=item RETURN

Returns a new version of the tree with the MD5 protein IDs replaced by FIG
feature IDs.

=back

=cut

sub md5_tree_to_fid_tree {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the tree and the map.
 my $tree = $args->{-tree}
   || Confess("No tree specified in md5_tree_to_fid_tree.");
 my $map = $args->{-map}
   || Confess("No map specified in md5_tree_to_fid_tree.");

 # Perform the conversion.
 my $retVal = AlignsAndTrees::md5_tree_to_fid_tree( $tree, $map );

 # Return the result.
 return $retVal;
}

=head3 trees_with_md5ID

    my $protHash =          $altObject->trees_with_md5ID({
                                -ids => [$prot1, $prot2, ...]
                            });

Return a list of the alignment/tree pairs containing each of the specified proteins.

=over 4

=item parameters

The parameter should be a reference to a hash with the following keys:

=over 8

=item -ids

Reference to a list of MD5 protein IDs.

=back

=item RETURN

Returns a reference to a hash mapping each incoming protein ID to a list of the
IDs for the trees containing that protein.

    $protHash = { $prot1 => [$alt1a, $alt1b, ...],
                  $prot2 => [$alt2a, $alt2b, ...],
                  ... };

=back

=cut

sub trees_with_md5ID {

 # Get the parameters.
 my ( $self, $args ) = @_;

 # Get the Sapling database.
 my $sap = $self->{db};

 # Create the return hash.
 my $retVal = {};

 # Get the list of incoming IDs.
 my $ids = ServerThing::GetIdList( -ids => $args );

 # Loop through the protein IDs, finding the trees.
 for my $id (@$ids) {
  $retVal->{$id} = AlignsAndTrees::trees_with_md5ID( $sap, $id );
 }

 # Return the result hash.
 return $retVal;
}

1;
