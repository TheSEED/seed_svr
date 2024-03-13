package StandaloneTrees;

use strict;
use FIG_Config;
use DBI;
use Data::Dumper;
use SeedUtils;
use FIG;
use gjoseqlib;
use gjonewicklib;
use File::Slurp;
use Storable qw(nfreeze thaw);

use base 'Class::Accessor';

StandaloneTrees->mk_accessors(qw(fig dbh data_dir));

sub new
{
    my($class, $fig) = @_;

    $fig //= new FIG;

    my $dbh = DBI->connect("DBI:mysql:database=pubseed_sapling_14;host=localhost", "seed");
#    my $dbh = DBI->connect("DBI:mysql:database=pubseed_sapling_14;host=arborvitae.cels.anl.gov", "seed");
    $dbh or die "Cannot connect to pubseed_sapling_14";

    my $data_dir = "/vol/public-pseed/FIGdisk/FIG/ATNG";

    my $self =
    {
	dbh => $dbh,
	fig => $fig,
	data_dir => $data_dir,
    };

    return bless $self, $class;
}

sub data_path
{
    my($self, $f) = @_;
    return $self->data_dir . "/$f";
}

sub read_alignment
{
    my($self, $id) = @_;
    my $file = $self->data_path("ali$id.fa");
    my @align = map { $_->[1] = ''; $_ } gjoseqlib::read_fasta($file);
    return @align;
}

sub aligns_with_pegID
{
    my($self, $id) = @_;

    my $md5 = $self->fig->md5_of_peg($id);

    my $res = $self->dbh->selectcol_arrayref(qq(SELECT from_link
						FROM Aligns
						WHERE to_link = ?), undef, $md5);
    return @$res;
}

sub peg_alignment_by_ID
{
    my($self, $alignID) = @_;

    my ( $md5_alignment, $md5_metadata ) = $self->md5_alignment_by_ID($alignID );
    $md5_alignment && $md5_metadata or return wantarray ? () : undef;
    my ( $peg_metadata, $md5ID_to_fidIDs_map ) = $self->map_md5_to_fid($md5_metadata );
    my $peg_alignment = md5_align_to_fid_align( $md5_alignment, $md5ID_to_fidIDs_map );

    wantarray ? ( $peg_alignment, $peg_metadata ) : $peg_alignment;
}


sub md5_alignment_by_ID
{
    my($self, $alignID) = @_;

    my @align = $self->read_alignment($alignID);
    
    wantarray ? (\@align, $self->md5_alignment_metadata($alignID)) : \@align;
}

sub md5_alignment_metadata
{
    my($self, $alignID) = @_;

    my $res = $self->dbh->selectall_arrayref(qq(SELECT sequence_id, to_link, len, begin, end, properties
						FROM Aligns
						WHERE from_link = ?), undef, $alignID);

    my $md = {};
    for my $ent (@$res)
    {
	my($seq, $to, $len, $begin, $end, $properties) = @$ent;
	$md->{$seq} = [$to, $len, $begin, $end, $properties];
    }
    return $md;
}

sub peg_tree_by_ID
{
    my($self, $treeID) = @_;
    $treeID  =~ s/^tree//;
    my ( $md5_tree, $md5_metadata ) = $self->md5_tree_by_ID($treeID );
    $md5_tree && $md5_metadata or return wantarray ? () : undef;
    
    my ( $peg_metadata, $md5ID_to_pegIDs_map ) = $self->map_md5_to_fid($md5_metadata );
    
    my $peg_tree = md5_tree_to_fid_tree( $md5_tree, $md5ID_to_pegIDs_map );

    wantarray ? ( $peg_tree, $peg_metadata ) : $peg_tree;
}



sub md5_tree_by_ID
{
    my($self, $treeID) = @_;

    my $file = $self->data_path("tree$treeID.nwk");
    my $frozen = $self->data_path("FROZEN/tree$treeID.frozen");

    my $tree;
    if (-f $frozen)
    {
	$tree = thaw(scalar read_file($frozen));
    }
    else
    {
	$tree = gjonewicklib::read_newick_tree($file);
    }
    wantarray ? ($tree, $self->md5_alignment_metadata($treeID)) : $tree;
	
}

my %md5_to_pegs;

sub md5s_to_pegs
{
    my($self, @md5) = @_;

    return wantarray ? () : {} if ! @md5;

    #  Remove the md5 ids already in the cache

    @md5 = grep { ! exists $md5_to_pegs{ $_ } } @md5;

    #  Get the remaining md5 ids from the Sapling:

    if ( @md5 )
    {
	for my $md5 (@md5)
	{
	    my @pegs = $self->fig->pegs_with_md5($md5);
	    push(@{$md5_to_pegs{$md5}}, $_) foreach @pegs;
	}
    }

    #  Return the whole cache

    wantarray ? %md5_to_pegs : \%md5_to_pegs;
}

my %peg_to_md5;
sub pegs_to_md5
{
    my($self, @pegID) = @_;

    return wantarray ? () : {} if ! @pegID;

    #  Remove the pegIDs already in the cache

    @pegID = grep { ! exists $peg_to_md5{ $_ } } @pegID;

    #  Get the remaining pegIDs from the Sapling:

    if ( @pegID )
    {
	my $ret = $self->fig->md5_of_peg_bulk(\@pegID);
	for my $peg (@pegID)
	{
	    $peg_to_md5{$peg} = $ret->{$peg};
	}
    }

    #  Return the whole cache

    wantarray ? %peg_to_md5 : \%peg_to_md5;
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
    my($self, @md5s) = @_;
    #  Remove the md5IDs already in the cache
    my @new_md5s = grep { ! exists $md5_function{ $_ } } @md5s;

    if ( @new_md5s )
    {
        #  Batch lookup of md5 ids returns hash of all currently known translations

        my $md5_to_pegs = $self->md5s_to_pegs( @new_md5s );

        #  Just get first peg (this is dangerous, but will become more consistent in future) 
        my %md5_to_one_peg = map { $md5_to_pegs->{ $_ } ? ( $_ => $md5_to_pegs->{ $_ }->[0] ) : () } @new_md5s;
        my @new_pegs       = values %md5_to_one_peg;
        my $new_peg_func   = $self->fig->function_of_bulk(\@new_pegs);

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


sub map_md5_to_fid
{
    my ( $self, $md5_metadata, $relaxed ) = @_;
    $md5_metadata && ref( $md5_metadata ) eq 'HASH'
        or return ();

    my %fid_metadata;
    my %fids_seen;
    my %md5ID_to_fidIDs_map;

    my $md5_to_pegs = $self->md5s_to_pegs(keys %$md5_metadata );

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

sub peg_to_md5
{
    my($self, $peg) = @_;
    return $self->fig->md5_of_peg($peg);
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
    my($self, @md5IDs) = @_;

    if (@md5IDs == 1 && ref($md5IDs[0]))
    {
	@md5IDs = @{$md5IDs[0]};
    }
	

    my $qs = join(",", map { "?" } @md5IDs);
    my $res = $self->dbh->selectall_arrayref(qq(SELECT from_link, sequence_id, to_link, len, begin, end, properties
						FROM Aligns
						WHERE to_link IN ($qs)), undef, @md5IDs);
    wantarray ? @$res : $res;
}

sub md5IDs_in_align
{
    my($self, $alignID) = @_;

    my %ids;    

    my $res = $self->dbh->selectall_arrayref(qq(SELECT to_link
						FROM Aligns
						WHERE from_link = ?), undef, $alignID);
    $ids{$_->[0]} = 1 foreach @$res;
    my @ids = keys %ids;
    wantarray ? @ids : \@ids;
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
    my($self, $md5, @alignIDs) = @_;

    $md5 or return wantarray ? () : [];

    my %keep = map { $_ => 1 } @alignIDs;

    # [ $alignID, $seqID, $md5, $peg_length, $trim_beg, $trim_end, $location_string ]
    my %metadata;
    foreach ( @{ $self->alignments_metadata_by_md5( [$md5] ) || [] } )
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
    my($self, $alignID) = @_;

    $alignID =~ s/^aln//;

    my @md5s = $self->md5IDs_in_align( $alignID );
    my $md5_function = $self->md5s_to_functions(  @md5s ) || {};

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
    my $self = shift;
    my $opts = ref $_[-1] eq 'HASH' ? pop @_ : {};
    my @md5s = @_ or die "Empty list of MD5 IDs\n";;

    my $minScore = $opts->{ minScore } || $opts->{ min_score } || 0;
    my $details  = $opts->{ details }  || $opts->{ full }      || 0;

    #    $ALITREserverO->get_projections( -ids => \@md5s, -minScore => $minScore, -details => $details );

    my %ret;
    for my $id (@md5s)
    {
	my $res1 = $self->dbh->selectall_arrayref(qq(SELECT from_link, gene_context, percent_identity, score
						     FROM ProjectsOnto
						     WHERE to_link = ? AND
						     score > ?), undef, $id, $minScore);
	my $res2 = $self->dbh->selectall_arrayref(qq(SELECT to_link, gene_context, percent_identity, score
						     FROM ProjectsOnto
						     WHERE from_link = ? AND
						     score > ?), undef, $id, $minScore);
	my @proj;
	for my $ent (@$res1, @$res2)
	{
	    my($match, $gene_context, $iden, $score) = @$ent;
	    if ($details)
	    {
		push(@proj, $ent);
	    }
	    else
	    {
		push(@proj, $match);
	    }
	}
	$ret{$id} = \@proj;
    }
    return \%ret;
}


1;
