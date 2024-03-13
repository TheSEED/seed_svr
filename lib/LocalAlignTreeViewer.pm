package LocalAlignTreeViewer;
#package CGIAlignTreeViewer;

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

use strict;
use HTML                 qw( fid_link
                             java_buttons
                           );
use SAPserver;
use SeedUtils            qw( genome_of );
use SeedAware            qw( open_tmp_file );
use AlignsAndTreesServer qw( aligns_with_pegID
                             get_md5_projections
                             md5s_to_pegs
                             peg_alignment_by_ID
                             peg_alignment_metadata
                             peg_to_md5
                             peg_tree_by_ID
                             pegs_to_md5
                             roles_in_align
                           );
use FIGgjo               qw( colorize_roles );
#  gd_tree_0 is invoked via a "require gd_tree", only when requested.
use gjoalign2html        qw( alignment_2_html_table
                             color_alignment_by_consensus
                             color_alignment_by_residue
                             repad_alignment
                           );
use gjonewicklib         qw( aesthetic_newick_tree
                             formatNewickTree
                             newick_relabel_nodes
                             newick_subtree
                             reroot_newick_to_midpoint_w
                             text_plot_newick
                           );

use StandaloneTrees;

use Data::Dumper;
use Time::HiRes qw( time );

sub run
{
    my( $fig, $cgi, $sapObject, $user, $url, $hidden_form_var ) = @_;

    my $atree = StandaloneTrees->new($fig);

    return help() if $cgi->param( 'page_help' );

    my $action        = $cgi->param( 'request' )       || ''; # assign & change focus are the only requests
    my $ali_tree_id   = $cgi->param( 'ali_tree_id' )   || '';
    my @ali_tree_ids  = $cgi->multi_param( 'at_ids' );
    my $align_format  = $cgi->param( 'align_format' )  || 'default';    # default || fasta || clustal
    my $align_id      = $cgi->param( 'align_id' );
    my $au            = $cgi->param( 'assign_using' );
    my $assign_using  = ( $au =~ /^Sap/i ) || ( ! $fig ) ? 'Sapling' : 'SEED';
    my @checked       = $cgi->multi_param( 'checked' );
    my $color_aln_by  = $cgi->param( 'color_aln_by' )  || 'residue';   # consensus || residue || none
    my $fid           = $cgi->param( 'fid' )           || '';
    my $from          = $cgi->param( 'from' )          || '';          # peg assignment to propagate
    my $rep_pegs      = $cgi->param( 'rep_pegs' )      || 'all';       # all || roles || dlit || paralog
  # my $show_aliases  = $cgi->param( 'show_aliases' )  || '';
    my $show_align    = $cgi->param( 'show_align' );
    my $show_tree     = $cgi->param( 'show_tree' );
    my $condense      = $cgi->param( 'condense' )      || 'default'; # default || none
    my $tree_format   = $cgi->param( 'tree_format' )   || 'default'; # default || newick || png
    my $tree_id       = $cgi->param( 'tree_id' );

    # The html will be assembled here.

    my @html = ();

    #------------------------------------------------------------------------------
    #  Convert the cgi paramater values to a local summary of the work to be done
    #------------------------------------------------------------------------------

    #  Let's see if we can work out missing values from other data:

    $fid         ||= $checked[0] if @checked == 1;
    $ali_tree_id ||= $align_id || $tree_id || '';
    $ali_tree_id   = '' if $action =~ /ali.* tree.* with.* prot/i;  #  Forced update of list
    if ( ( ! $ali_tree_id ) && ( ! @ali_tree_ids ) && $fid )
    {
        @ali_tree_ids = $atree->aligns_with_pegID( $fid );
    }
    $ali_tree_id ||= $ali_tree_ids[0] if @ali_tree_ids == 1;

    #  Move alignment and tree selection information into one id and two booleans

    $show_align ||= $align_id;
    $show_tree  ||= $tree_id;

    #------------------------------------------------------------------------------
    #  We have the analysis paramaters.  Put them in a local hash so they can be passed to
    #  subroutines.
    #------------------------------------------------------------------------------

    my $data = {};

    $data->{atree} = $atree;

    $data->{ fig }  =  $fig;
    $data->{ sap }  =  $sapObject;
    $data->{ cgi }  =  $cgi;
    $data->{ html } = \@html;
    $data->{ user } =  $user;
    $data->{ url }  =  $url;

    $data->{ action }        =  $action;
    $data->{ ali_tree_id }   =  $ali_tree_id;
    $data->{ ali_tree_ids }  = \@ali_tree_ids;
    $data->{ align_format }  =  $align_format;
    $data->{ assign_using }  =  $assign_using;
    $data->{ can_assign }    =  $user && ( $assign_using =~ /SEED/i );
    $data->{ checked }       = \@checked;
    $data->{ color_aln_by }  =  $color_aln_by;
    $data->{ fid }           =  $fid;
    my $form_name = 'alignment_and_tree';
    $data->{ form_name }     =  $form_name;
    $data->{ from }          =  $from;
    $data->{ rep_pegs }      =  $rep_pegs;
  # $data->{ show_aliases }  =  $show_aliases;
    $data->{ show_align }    =  $show_align;
    $data->{ show_tree }     =  $show_tree;
    $data->{ condense }      =  $condense;
    $data->{ tree_format }   =  $tree_format;

    #------------------------------------------------------------------------------
    #  Start the page:
    #------------------------------------------------------------------------------

    my( $this_html, $title ) = page_head_html( $data );

    #------------------------------------------------------------------------------
    #  Deal with assignments:
    #------------------------------------------------------------------------------

    make_assignments( $data )  if ( $data->{ action } =~ /^Assign/i );

    #------------------------------------------------------------------------------
    #  Change the focus peg:
    #------------------------------------------------------------------------------

    $fid = $from  if ( $data->{ action } =~ /focus/i && $from );

    #------------------------------------------------------------------------------
    #  Special form for help information:
    #------------------------------------------------------------------------------

    push @html, $cgi->start_form( -method => 'post',
				  -action => $url,
				  -name   => 'alignment_and_tree_help',
				  -target => 'Alignment and Tree Help'
				);
    push @html, $hidden_form_var;
    push @html, $cgi->submit( -name => 'page_help', -value => 'Alignment and Tree Help' ), ' (in new window)';
    push @html, $cgi->end_form, $cgi->br;

    #------------------------------------------------------------------------------
    #  Start the form:
    #------------------------------------------------------------------------------

    push @html, $cgi->start_form( -method => 'post',
				  -action => $url,
				  -name   => $form_name
				);
    push @html, $hidden_form_var;

    #------------------------------------------------------------------------------
    #  Alignment and tree format controls:
    #------------------------------------------------------------------------------

    add_general_options( $data );

    #------------------------------------------------------------------------------
    #  Collect all of the necessary alignment and/or tree data:
    #------------------------------------------------------------------------------

    if ( $data->{ ali_tree_id } && ( $data->{ show_align } || $data->{ show_tree } ) )
    {
	compile_alignment_and_tree_data( $data );
    }

    #------------------------------------------------------------------------------
    #  Alignment dispaly
    #------------------------------------------------------------------------------

    show_alignment( $data )  if ( $data->{ ali_tree_id } && $data->{ show_align } );

    #------------------------------------------------------------------------------
    #  Tree display
    #------------------------------------------------------------------------------

    show_tree( $data )  if ( $data->{ ali_tree_id } && $data->{ show_tree } );

    #------------------------------------------------------------------------------
    #  Select alignments and trees with given fid
    #------------------------------------------------------------------------------

    show_alignments_and_trees_with_fid( $data )  if ( ! $data->{ ali_tree_id } );

    #------------------------------------------------------------------------------
    #  Finish form and body content
    #------------------------------------------------------------------------------

    push @html, $cgi->end_form, '';

    #------------------------------------------------------------------------------
    #  Report the output
    #------------------------------------------------------------------------------

    return ( join( "\n", @html, '' ), $title );
}



#==============================================================================
#  Only subroutines below
#==============================================================================
#  This is a sufficient set of escaping for text in HTML (function and alias):
#
#     $html = html_esc( $text )
#------------------------------------------------------------------------------

sub html_esc { local $_ = $_[0]; s/\&/&amp;/g; s/\>/&gt;/g; s/\</&lt;/g; $_ }


#===============================================================================
#  Start the HTML
#===============================================================================

sub page_head_html
{
    my ( $data ) = @_;
    my $html = $data->{ html } || [];

    my $ali_tree_id = $data->{ ali_tree_id };
    my $fid         = $data->{ fid };
    my $show_align  = $data->{ show_align };
    my $show_tree   = $data->{ show_tree };

    my $title;
    if ( $show_align && $ali_tree_id )
    {
        if ( $show_tree ) { $title = "Protein Alignment $ali_tree_id" }
        else              { $title = "Protein Alignment and Tree $ali_tree_id" }
    }
    elsif ( $show_tree && $ali_tree_id )
    {
        $title = "Protein Tree $ali_tree_id";
    }
    else
    {
        if ( $fid ) { $title = "Protein Alignment and Tree Selector for '$fid'" }
        else        { $title = "Protein Alignment and Tree Selector" }
    }

    #  This stuff is because different browsers render the contents differently.
    #  $height and $lsize are essential for the CSS here-document below.

    my $agent  = $ENV{ HTTP_USER_AGENT } || '';
    my $height = $agent =~ /Safari/i  ? '110%'
               : $agent =~ /Firefox/i ? '100%'
               :                        '100%';
    my $lsize  = $agent =~ /Safari/i  ? '160%'
               : $agent =~ /Firefox/i ? '130%'
               :                        '140%';

    push @$html, <<"End_of_CSS";
<!--
  --  Different browsers handle the layout of the line-drawing characters
  --  differently. This is An attempt to make it look okay.
  -->
<STYLE Type="text/css">
  /* Support for HTML printer graphics tree */
  DIV.tree {
    border-spacing:  0px;
    font-size:      100%;
    line-height: $height;
    white-space:  nowrap;
  }
  DIV.tree A {
    color:          black;
    text-decoration: none;
  }
  DIV.tree PRE {
    padding:      0px;
    margin:       0px;
    font-size: $lsize;
    display:   inline;
  }
  DIV.tree INPUT {
    padding: 0px;
    margin:  0px;
    height: 10px;    /* ignored by Firefox */
    width:  10px;    /* ignored by Firefox */
  }
  DIV.tree SPAN.w {  /* used for tree white space */
    color: white;
  }
</STYLE>

End_of_CSS

    push @$html, <<"End_of_JavaScript";
<SCRIPT Src="./Html/css/FIG.js" Type="text/javascript"></SCRIPT>

<SCRIPT Language="JavaScript">
//
// Tree Tip Selection Support (includin Undo and Redo)
//

function setAll( myForm, state )
{
    // Adjust the undo stack
    if ( ! myForm.undo ) { myForm.undo = new Array(); myForm.nUndo = 0 }
    while ( myForm.undo.length > myForm.nUndo ) { myForm.undo.pop() }
    // Create the new entry
    var history = new Array();
    myForm.undo.push( history ); myForm.nUndo++;
    // Do the requested selection
    if ( state == undefined ) { state = true }
    for ( var i = 0; i < myForm.checked.length; i++ )
    {
        var box = myForm.checked[i];
        if ( box.checked != state )
        {
            box.checked = state;
            history.push( [box,state] );
        }
    }
    setUndo( myForm, true );
    setRedo( myForm, false );
}

function setSubtree( myForm, nodeNum, state )
{
    if ( ! myForm.nodes[nodeNum] ) { return }
    // Adjust the undo stack
    if ( ! myForm.undo ) { myForm.undo = new Array(); myForm.nUndo = 0 }
    while ( myForm.undo.length > myForm.nUndo ) { myForm.undo.pop() }
    // Create the new entry
    var history = new Array();
    myForm.undo.push( history ); myForm.nUndo++;
    // Do the requested selection
    if ( state == undefined ) { state = true }
    setSubtree2( myForm, nodeNum, state, history );
    setUndo( myForm, true );
    setRedo( myForm, false );
}

function setSubtree2( myForm, nodeNum, state, history )
{
    var node = myForm.nodes[nodeNum];
    var desc = node[0];
    if ( desc ) { for ( var i = 0; i < desc.length; i++ ) { setSubtree2( myForm, desc[i], state, history ) } }
    // Set the tip checkbox
    if ( ! ( node[1] == null ) ) { setTip( myForm, node[1], state, history ) }
}

function setTip( myForm, tipNum, state, history )
{
    var box = myForm.checked[tipNum];
    if ( ! box ) { return }
    if ( box.checked == state ) { return }
    box.checked = state;
    history.push( [box,state] );
}

function toggleTip( myForm, tipBox )
{
    if ( ! myForm ) { return }
    var state = tipBox.checked;

    // Adjust the undo stack
    if ( ! myForm.undo ) { myForm.undo = new Array(); myForm.nUndo = 0 }
    while ( myForm.undo.length > myForm.nUndo ) { myForm.undo.pop() }
    // Create the new entry
    myForm.undo.push( [[tipBox,state]] );  // A list with a single change
    myForm.nUndo++;
    setUndo( myForm, true );
    setRedo( myForm, false );
}

function undoSelect( myForm )
{
    if ( ! myForm.undo || ! myForm.nUndo ) { return }
    var pairs = myForm.undo[--myForm.nUndo];
    for ( var i = 0; i < pairs.length; i++ )
    {
        var pair = pairs[i];
        pair[0].checked = ! pair[1];
    }
    setUndo( myForm, myForm.nUndo > 0 );
    setRedo( myForm, true );
}

function redoSelect( myForm )
{
    if ( ! myForm.undo || myForm.undo.length <= myForm.nUndo ) { return }
    var pairs = myForm.undo[myForm.nUndo++];
    for ( var i = 0; i < pairs.length; i++ )
    {
        var pair = pairs[i];
        pair[0].checked = pair[1];
    }
    setUndo( myForm, true );
    setRedo( myForm, myForm.undo.length > myForm.nUndo );
}

function setUndo( myForm, state )
{
   for ( var i = 0; i < myForm.UndoBtn.length; i++ ) { myForm.UndoBtn[i].disabled = ! state }
}

function setRedo( myForm, state )
{
   for ( var i = 0; i < myForm.RedoBtn.length; i++ ) { myForm.RedoBtn[i].disabled = ! state }
}
</SCRIPT>

End_of_JavaScript

    return( $html, $title );
}


#===============================================================================
#  Make requested assignments.
#===============================================================================

sub make_assignments
{
    my ( $data ) = @_;

    my $fig  = $data->{ fig };
    my $sap  = $data->{ sap };
    my $cgi  = $data->{ cgi };
    my $html = $data->{ html };
    my $user = $data->{ user };
    my $from = $data->{ from };

    my $func;
    if ( defined( $from ) && ( $func = $fig->function_of( $from, $user ) ) && @{ $data->{ checked } } )
    {
        $func =~ s/\s+\#[^\#].*$//;       #  Remove single hash comments

        if (  $data->{ assign_using } =~ m/SEED/i && $fig )
        {
            my @pegs = @{ $data->{ checked } };
            my $assign_opts = { annotation => "Assignment projected from $from based on tree proximity"
                              };
            my ( $nsucc, $nfail, $nmoot ) = $fig->assign_function( \@pegs, $user, $func, $assign_opts );
            push @$html, $cgi->h3( "$nsucc assignments changed." )    if $nsucc;
            push @$html, $cgi->h3( "$nmoot assignments were moot." )  if $nmoot;
            push @$html, $cgi->h3( "$nfail assignments failed." )     if $nfail;
        }
        elsif ( 0 )
        {
            #  We currently do not have assignment outside the SEED environment
        }
    }
    else
    {
        push @$html, $cgi->h3( 'Cannot assign with no radio button selected.' )  if ! defined( $from );
        push @$html, $cgi->h3( 'Cannot assign a null function.' )                if $from && ! $func;
        push @$html, $cgi->h3( 'Cannot assign with no check boxes selected.' )   if ! @{ $data->{ checked } };
    }
}


#===============================================================================
#  Push the general page options into the html.
#===============================================================================

sub add_general_options
{
    my ( $data ) = @_;
    my $cgi  = $data->{ cgi };
    my $html = $data->{ html } || [];
    my $user = $data->{ user };

    if ( @{ $data->{ checked } } && ! $data->{ show_tree } )
    {
        push @$html, $cgi->hidden( -name => 'checked', -value => $data->{ checked } );
    }

    if ( $data->{ ali_tree_id } )
    {
        push @$html, $cgi->hidden( -name => 'ali_tree_id', -value => $data->{ ali_tree_id } );
    }

    if ( $user )
    {
        push @$html, $cgi->hidden( -name => 'user', -value => $user );
    }
    else
    {
        push @$html, 'SEED user: ',
                      $cgi->textfield( -name     => 'user',
                                       -value    => '',
                                       -size     => 32,
                                       -override => 1
                                     ),
                      $cgi->br;
    }
    
    if ( $data->{ ali_tree_id } || $data->{ fid } )
    {
        push @$html, 'Focus protein ID: ';
    }
    else
    {
        push @$html, $cgi->h2( 'Enter a SEED protein id: ' );
    }

    push @$html, $cgi->textfield( -name => "fid", -size => 32, -value => $data->{ fid } ),
                 $cgi->submit( -name => 'request', -value => 'List all alignments and trees with this protein' ),
                 $cgi->br;

    if ( ! ( $data->{ show_align } || $data->{ show_tree } ) )
    {
        push @$html, '<SPAN Style="color: #CC0000">',
                     $cgi->h2( 'Neither alignment nor tree are selected below.  Please select at least one.' ),
                     '</SPAN>';
    }

    push @$html, $cgi->checkbox( -name     => 'show_align',
                                 -label    => 'Show alignment',
                                 -override => 1,
                                 -checked  => $data->{ show_align }
                               ),
                 '&nbsp;',
                 $cgi->checkbox( -name     => 'show_tree',
                                 -label    => 'Show tree',
                                 -override => 1,
                                 -checked  => $data->{ show_tree }
                               ),
                 $cgi->br, $cgi->br;

    push @$html, $data->{ can_assign } ? 'Use for functions and assignments: '
                                       : 'Use for functions: ';
    push @$html, ( map { "&nbsp;$_" }
                   $cgi->radio_group( -name     => 'assign_using',
                                      -override => 1,
                                      -values   => [ 'Sapling', 'SEED' ],
                                      -default  => $data->{ assign_using }
                                    )
                 ),
                 $cgi->br;

    push @$html, 'Of identical sequences, show: ',
                 ( map { "&nbsp;$_" }
                   $cgi->radio_group( -name     => 'condense',
                                      -override => 1,
                                      -values   => [ 'default', 'all' ],
                                      -default  => $data->{ condense } || 'default'
                                    )
                 ),
                 $cgi->br;

    push @$html, $cgi->br,
                 'Color alignment by: ',
                 ( map { "&nbsp;$_" }
                   $cgi->radio_group( -name     => 'color_aln_by',
                                      -override => 1,
                                      -values   => [ 'consensus', 'residue', 'none' ],
                                      -default  => $data->{ color_aln_by }
                                    )
                 ),
                 $cgi->br;

    push @$html, 'Alignment format: ',
                 ( map { "&nbsp;$_" }
                   $cgi->radio_group( -name     => 'align_format',
                                      -override => 1,
                                      -values   => [ 'default', 'fasta', 'clustal' ],
                                      -default  => $data->{ align_format } || 'default'
                                    )
                 ),
                 $cgi->br, $cgi->br;

    push @$html, 'Tree format: ',
                 ( map { "&nbsp;$_" }
                   $cgi->radio_group( -name     => 'tree_format',
                                      -override => 1,
                                      -values   => [ 'default', 'newick', 'png' ],
                                      -default  => $data->{ tree_format } || 'default'
                                    )
                 ),
                 $cgi->br;

  # push @$html, $cgi->checkbox( -name     => 'show_aliases',
  #                              -label    => 'Show aliases in tree',
  #                              -override => 1,
  #                              -checked  => $data->{ show_aliases }
  #                            ),
  #              $cgi->br;

    push @$html, $cgi->br,
                 $cgi->submit( -name => 'request', -value => 'Update' ),
                 $cgi->br;

    return @$html if wantarray;
}


#------------------------------------------------------------------------------
#  Compile all necessary data for alignments and trees.
#  The per sequence metadata are:
#
#      [ $peg_id, $peg_length, $trim_beg, $trim_end, $location_string ]
#
#------------------------------------------------------------------------------

sub compile_alignment_and_tree_data
{
    my ( $data ) = @_;

    ( $data->{ ali_tree_id } && ( $data->{ show_align } || $data->{ show_tree } ) )
        or return 0;

    my $html = $data->{ html } || [];
    my $sap  = $data->{ sap };
    my $cgi  = $data->{ cgi };
    my $fig  = $data->{ fig };
    my $fid  = $data->{ fid };
    my $user = $data->{ user };

    my $align = [];
    my $tree  = undef;
    my $metaH = {};

    if ( $data->{ show_align } )
    {
        ( $align, $metaH ) = $data->{atree}->peg_alignment_by_ID( $data->{ ali_tree_id } );
    }
    if ( $data->{ show_tree } )
    {
        ( $tree, $metaH ) = $data->{atree}->peg_tree_by_ID( $data->{ ali_tree_id } );
    }

    $metaH && %$metaH
        or push @$html, $cgi->h2( "No data for alignment and tree '$data->{ali_tree_id}'." );

    my @uids = keys %$metaH;    # Ids of alignment line and tree tips
    my %fid_of_uid = map { $_ => $metaH->{$_}->[0] } @uids;

    my %peg_seen = {};
    my @fids = grep { ! $peg_seen{$_}++ } values %fid_of_uid;

    #--------------------------------------------------------------------------
    #  In case of SEED, remove ids that do not exist locally:
    #--------------------------------------------------------------------------

    if ( @fids && $data->{ assign_using } =~ /^SEED/i && $fig )
    {
        #  Find the fids that are actually present in this SEED.

        # my %have_fid = map { $fig->is_real_feature( $_ ) ? ( $_ => 1 ) : () } @fids;
        my %have_fid = map { $_ => 1 } $fig->is_real_feature_bulk( \@fids );

        @uids = grep { $have_fid{ $fid_of_uid{ $_ } } } @uids;
        @fids = grep { $have_fid{ $_ } } @fids;

        my %have_uid = map { $_ => 1 } @uids;

        @$align = grep { $have_uid{ $_->[0] } } @$align;
        $tree   = gjonewicklib::newick_subtree( $tree, @uids ) if $tree;
        %$metaH = map { $have_uid{ $_ } ? ( $_ => $metaH->{$_} ) : () }
                  keys %$metaH;
    }

    #--------------------------------------------------------------------------
    #  Find the current functions and organism names:
    #--------------------------------------------------------------------------

    my $fid_funcH = {};
    my $orgH      = {};
    if ( @fids && $data->{ assign_using } =~ /^SEED/i && $fig )
    {
        $fid_funcH = $fig->function_of_bulk( \@fids );
        foreach my $peg ( @fids ) { $orgH->{ $peg } = $fig->org_of( $peg ) }
    }
    elsif ( @fids )
    {
        $sap     ||= SAPserver->new();
        $fid_funcH = $sap->ids_to_functions( -ids => \@fids ) || {};
        $orgH      = $sap->ids_to_genomes( -ids => \@fids, -name => 1 ) || {};
    }

    #--------------------------------------------------------------------------
    #  Condense identical sequences with identical functions
    #
    #  Keep genome of focus peg;
    #  For groups of identical sequences
    #      Keep one genome of each assignment;
    #      Keep same genome across paralog sets;
    #  We can restrict the calculation to md5 groups with more than one fid;
    #--------------------------------------------------------------------------

    #  Because dlits also need this converstion, we will put it here
    my $pegs_to_md5 = $data->{atree}->pegs_to_md5(  @fids );

    my %keep_genome;
    my %represents;
    my $condense = $data->{ condense } || 'default';
    if ( $condense =~ /^def/i )
    {
        #  Summarize the data as [ fid, md5, gid, func ] for each fid:
        my @tuples = map { [ $_,
                             $pegs_to_md5->{ $_ },
                             SeedUtils::genome_of( $_ ),
                             $fid_funcH->{ $_ }
                           ]
                         }
                     @fids;

        #  Data on a fid
        my %tuple_by_fid = map { $_->[0] => $_ } @tuples;

        #  Data on a set of identical sequences
        my %tuples_by_md5;
        foreach ( @tuples ) { push @{ $tuples_by_md5{ $_->[1] } }, $_ }

        #  Data on sequences from a genome
        my %tuples_by_gid;
        foreach ( @tuples ) { push @{ $tuples_by_gid{ $_->[2] } }, $_ }

        #  Number of different roles from a genome (priority for showing)
        my %n_roles;
        my %roles;
        foreach ( values %tuples_by_gid )
        {
            %roles = map { $_->[3] => 1 } @$_;
            $n_roles{ $_->[0]->[2] } = scalar keys %roles;
        }

        #  Work out the genomes to keep:

        #  The genome of every singleton sequence must be kept (nothing else
        #  can represent its possible paralogs)
        foreach ( grep { @$_ == 1 } values %tuples_by_md5 )
        {
            $keep_genome{ $_->[0]->[2] } = 1;
        }

        #  The focus genome must be kept, and it top priority as rep
        $keep_genome{ SeedUtils::genome_of( $fid ) || '' } = 2;

        #  Now work through the remaining md5 groups
        foreach ( sort { @$b <=> @$a } grep { @$_ > 1 } values %tuples_by_md5 )
        {
            #  Group the sequences by function
            my %by_func;
            foreach my $tuple ( @$_ ) { push @{ $by_func{ $tuple->[3] } }, $tuple }

            foreach ( sort { @$b <=> @$a } values %by_func )
            {
                my @equivs = sort { $keep_genome{$b->[2]} cmp $keep_genome{$a->[2]}
                                 || $n_roles{$b->[2]} <=> $n_roles{$a->[2]}
                                  }
                             @$_;
                my $rep = shift @equivs;
                $keep_genome{ $rep->[2] } ||= 1;
                #  List of genomes to be represented:
                @equivs = grep { ! $keep_genome{ $_->[2] } } @equivs;
                #  List of fids represented
                $represents{ $rep->[0] } = [ map { $_->[0] } @equivs ];
            }
        }
        #  We should cleanse the value in represents for late additions to
        #  the keep list, but this is not critical, so put if off for later
        #  improvements.

        #  Okay, now we filter the tree prune the master lists of uids and
        #  fids, and prune the alignment and tree.
        if ( keys %represents )
        {
            my %keep_fid = map  { $_->[0] => 1 }
                           map  { @{ $tuples_by_gid{ $_ } } }
                           grep { $_ && $keep_genome{ $_ } && $tuples_by_gid{$_} }
                           keys %keep_genome;

            @uids = grep { $keep_fid{ $fid_of_uid{ $_ } } } @uids;
            @fids = grep { $keep_fid{ $_ } } @fids;

            my %keep_uid = map { $_ => 1 } @uids;

            @$align = grep { $keep_uid{ $_->[0] } } @$align;
            $tree   = gjonewicklib::newick_subtree( $tree, @uids ) if $tree;
        }
    }

    #--------------------------------------------------------------------------
    #  Aliases
    #--------------------------------------------------------------------------

    # my $aliasH = {};
    # if ( $data->{ show_aliases } ) { 0 }

    #--------------------------------------------------------------------------
    #  Build a hash from fid to dlit PMID list
    #--------------------------------------------------------------------------

    my $dlitH;
    if ( $data->{ assign_using } =~ /^SEED/i && $fig )
    {
        # Collect all available dlits by md5
        my %dlits_by_md5;
        foreach ( @{ $fig->all_dlits() } )
        {
            push @{ $dlits_by_md5{ $_->[1] } }, $_->[2];  # Just the PubMed ID
        }
        # Index by peg those that apply to the alignment and/or tree
        my %dlits;
        foreach ( @fids )
        {
            my $dlits = $dlits_by_md5{ $pegs_to_md5->{ $_ } };
            $dlits{ $_ } = $dlits if $dlits;
        }
        $dlitH = \%dlits;
    }
    else
    {
        $dlitH = $sap->dlits_for_ids( -ids => \@fids );
    }

    #--------------------------------------------------------------------------
    #  Projections from focus peg's md5 to other pegs (by their md5):
    #
    #  $proj->{ $md5a } = [ $md5b, $n_shared, $identity, $score ]
    #--------------------------------------------------------------------------

    #  Get the projections for the focus peg md5
    my $md5   = $data->{atree}->peg_to_md5( $fid ) || '';
    my $projH = $data->{atree}->get_md5_projections( $md5, { details => 1 } ) || {};
    my @projs = @{ $projH->{ $md5 } || [] };

    #  What are the md5s in the data?
    my %have_md5 = map { $_ ? ( $_ => 1 ) : () } map { $pegs_to_md5->{$_} } @fids;
    #  Expand the projection md5 values to pegs
    my @proj_md5s   = grep { $have_md5{$_} } map { $_->[0] } @projs;
    my $md5_to_pegs = $data->{atree}->md5s_to_pegs( $md5, @proj_md5s );

    my %projection;
    foreach my $proj ( @projs )
    {
        my $md5b = shift @$proj;
        my @pegs = @{ $md5_to_pegs->{ $md5b } || [] };
        foreach ( @pegs ) { $projection{ $_ } = $proj }
    }

    #  Projections to identical sequences

    foreach ( @{ $md5_to_pegs->{ $md5 } || [] } ) { $projection{$_} = [10,100,1] }

    #--------------------------------------------------------------------------
    #  Subsystems html for labels:
    #--------------------------------------------------------------------------

    my %subsyst = map { $_ => subsystems( $fig, $cgi, $sap, $_ ) }
                  @fids;

    #--------------------------------------------------------------------------
    #  Coverage color for each sequence:
    #--------------------------------------------------------------------------

    my %cov_color = map { $_ => region_color( ( @{ $metaH->{$_} } )[2,3,1] ) }
                    @uids;

    #--------------------------------------------------------------------------
    #  Put in data hash
    #--------------------------------------------------------------------------

    # $data->{ alias }       =  $aliasH;
    $data->{ align }       =  $align;
    $data->{ cov_color }   = \%cov_color;
    $data->{ dlits }       =  $dlitH;
    $data->{ fid_func }    =  $fid_funcH;
    $data->{ fid_of_uid }  = \%fid_of_uid;
    $data->{ fids }        = \@fids;
    $data->{ md5 }         =  $md5;
    $data->{ org }         =  $orgH;
    $data->{ pegs_to_md5 } =  $pegs_to_md5;
    $data->{ projects }    = \%projection;
    $data->{ represents }  = \%represents;
    $data->{ seq_meta }    =  $metaH;
    $data->{ subsyst }     = \%subsyst;
    $data->{ tree }        =  $tree;
    $data->{ uids }        = \@uids;

    return @$html if wantarray;
}


#==============================================================================
#  Show an alignment
#==============================================================================

sub show_alignment
{
    my ( $data ) = @_;
    my $html = $data->{ html } || [];
    my $cgi  = $data->{ cgi };

    ( $data->{ ali_tree_id } && $data->{ show_align } ) or return;

    my $align = $data->{ align };
    $align && @$align
        or push @$html, $cgi->h2( "No data for alignment '$data->{ali_tree_id}'." );

    #  This defines the ordering.
    my @seq_ids = map { $_->[0] } @$align;

    push @$html, $cgi->h2( "Alignment $data->{ali_tree_id}" ) . "\n";

    my $fid_of_uid = $data->{ fid_of_uid };
    my $fid_func   = $data->{ fid_func };
    my $org        = $data->{ org };

    if ( $align && @$align && ( $data->{ align_format } =~ /^fasta/i ) )
    {
        my ( $id, $peg );
        my %def = map { $id = $_->[0];
                        $peg = $fid_of_uid->{ $id };
                        $id => join( ' ', $id,
                                          ( $fid_func->{ $id } ? $fid_func->{$id} : () ),
                                          ( $org->{ $id }      ? "[$org->{$id}]"    : () )
                                   )
                      }
                  @$align;

        push @$html, join( "\n",
                          "<PRE>",
                          ( map { ( ">$def{$_->[0]}", $_->[2] =~ m/(.{1,60})/g ) } @$align ),
                          "</PRE>\n"
                        );
    }

    elsif ( $align && @$align && ( $data->{ align_format } =~ /^clustal/i ) )
    {
        push @$html, "<PRE>\n", &to_clustal( $align ), "</PRE>\n";
    }

    elsif ( $align && @$align )
    {
        my ( $align2, $legend );

        #  Color by residue type:

        if ( $data->{ color_aln_by } eq 'residue' )
        {
            my %param1 = ( align => $align, protein => 1 );
            $align2 = gjoalign2html::color_alignment_by_residue( \%param1 );
        }

        #  Color by consensus:

        elsif ( $data->{ color_aln_by } ne 'none' )
        {
            my %param1 = ( align => $align );
            ( $align2, $legend ) = gjoalign2html::color_alignment_by_consensus( \%param1 );
        }

        #  No colors:

        else
        {
            $align2 = gjoalign2html::repad_alignment( $align );
        }

        #  Add organism names:

        foreach ( @$align2 ) { $_->[1] = $org->{ $_->[0] || '' } }

        #  Build a tool tip with organism names and functions:

        my %tips = map { $_ => [ $_, join( $cgi->hr, $org->{ $_ }, $fid_func->{ $_ } ) ] }
                   map { $_->[0] }
                   @$align2;
        $tips{ 'Consen1' } = [ 'Consen1', 'Primary consensus residue' ];
        $tips{ 'Consen2' } = [ 'Consen2', 'Secondary consensus residue' ];

        my %param2 = ( align   => $align2,
                       tooltip => \%tips
                     );
        $param2{ legend } = $legend if $legend;

        push @$html, join( "\n",
                           scalar gjoalign2html::alignment_2_html_table( \%param2 ),
                           $cgi->br,
                         );
    }

    return @$html if wantarray;
}


#------------------------------------------------------------------------------
#  Clustal format alignment
#------------------------------------------------------------------------------
sub to_clustal
{
    my( $alignment ) = @_;

    my($tuple,$seq,$i);
    my $len_name = 0;
    foreach $tuple ( @$alignment )
    {
        my $sz = length( $tuple->[0] );
        $len_name = ($sz > $len_name) ? $sz : $len_name;
    }

    my @seq  = map { $_->[2] } @$alignment;
    my $seq1 = shift @seq;
    my $cons = "\377" x length($seq1);
    foreach $seq (@seq)
    {
        $seq  = ~($seq ^ $seq1);
        $seq  =~ tr/\377/\000/c;
        $cons &= $seq;
    }
    $cons =~ tr/\000/ /;
    $cons =~ tr/\377/*/;

    push(@$alignment,["","",$cons]);

    my @out = ();
    for ($i=0; ($i < length($seq1)); $i += 50)
    {
        foreach $tuple (@$alignment)
        {
            my($id,undef,$seq) = @$tuple;
            my $line = sprintf("\%-${len_name}s %s\n", $id, substr($seq,$i,50));
            push(@out,$line);
        }
        push(@out,"\n");
    }
    return join("","CLUSTAL W (1.8.3) multiple sequence alignment\n\n\n",@out);
}


#==============================================================================
#  Tree:
#==============================================================================

sub show_tree
{
    my ( $data ) = @_;

    my $html = $data->{ html } || [];
    my $user = $data->{ user };
    my $cgi  = $data->{ cgi };

    my $tree = $data->{ tree };
    if ( ! $tree )
    {
        push @$html, $cgi->h2( "No data for tree '$data->{ali_tree_id}'." );
        return wantarray ? @$html : ();
    }

    push @$html, $cgi->h2( "Tree $data->{ali_tree_id}" ) . "\n"  if $tree;

    # my $alias      = $data->{ alias }      || {};
    my $can_assign = $data->{ can_assign } ||  0;
    my $cov_color  = $data->{ cov_color }  || {};
    my $dlits      = $data->{ dlits }      || {};
    my $fid_func   = $data->{ fid_func }   || {};
    my $fid_of_uid = $data->{ fid_of_uid } || {};
    my $form_name  = $data->{ form_name };
    my $org        = $data->{ org }        || {};
    my $proj       = $data->{ projects }   || {};
    my $represents = $data->{ represents } || {};
    my $subsyst    = $data->{ subsyst }    || {};

    #------------------------------------------------------------------
    #  Newick tree
    #------------------------------------------------------------------
    if ( $tree && ( $data->{ tree_format } =~ /^newick/i ) )
    {
        push @$html, "<PRE>\n" . gjonewicklib::formatNewickTree( $tree ) . "</PRE>\n";
    }

    #------------------------------------------------------------------
    #  PNG tree
    #------------------------------------------------------------------
    elsif ( $tree && ( $data->{ tree_format } =~ /^png/i ) )
    {
        my $okay;
        eval { require gd_tree_0; $okay = 1 };
        my $fmt;
        if ( $okay && ( $fmt = ( gd_tree::gd_has_png() ? 'png'  :
                                 gd_tree::gd_has_jpg() ? 'jpeg' :
                                                         undef
                               ) ) )
        {
            #------------------------------------------------------------------
            #  Formulate the desired labels
            #------------------------------------------------------------------
            my %labels;
            foreach my $id ( @{ $data->{ uids } } )
            {
                my   $peg = $fid_of_uid->{ $id };
                my   @label;
                push @label, $id;
                push @label, $fid_func->{ $peg }          if $fid_func->{ $peg };
                push @label, "[$org->{$peg}]"             if $org->{ $peg };
                # push @label, html_esc( $alias->{ $peg } ) if $alias->{ $peg };
        
                $labels{ $id } = join( ' ', @label );
            }

            #------------------------------------------------------------------
            #  Relabel the tips, midpoint root, pretty it up and draw
            #  the tree as printer plot
            #
            #  Adjustable parameters on text_plot_newick:
            #
            #     @lines = text_plot_newick( $node, $width, $min_dx, $dy )
            #------------------------------------------------------------------
            my $tree2 = gjonewicklib::newick_relabel_nodes( $tree, \%labels );
            my $tree3 = gjonewicklib::reroot_newick_to_midpoint_w( $tree2 );
        
            $tree = gjonewicklib::aesthetic_newick_tree( $tree3 );
            my $options = { thickness =>  2,
                            dy        => 15,
                          };
            my $gd = gd_tree::gd_plot_newick( $tree, $options );

            my ( $TreeFH, $file ) = SeedAware::open_tmp_file( 'align_and_tree', $fmt, $FIG_Config::temp );
            binmode $TreeFH;
            print   $TreeFH $gd->$fmt;
            close   $TreeFH;
            chmod   0644, $file;

            my $name = $file;
            $name =~ s/^.*\///;
            my $url = "$FIG_Config::temp_url/$name";
            push @$html, "<BR />\n"
                       . "<IMG Src='$url' Border=0>\n"
                       . "<BR />\n";
        }
        else
        {
            push @$html, "<H3>Failed to convert tree to PNG.  Sorry.</H3>\n"
                       . "<H3>Please choose another format above.</H3>\n";
        }
    }

    #------------------------------------------------------------------
    #  Printer plot tree
    #------------------------------------------------------------------
    else
    {
        #------------------------------------------------------------------
        #  Midpoint root, and use aesthetic order.
        #  This determines checkbox order on page, and hence their index
        #  numbers (which is needed by the JavaScript).
        #------------------------------------------------------------------

        $tree = gjonewicklib::reroot_newick_to_midpoint_w( $tree );
        $tree = gjonewicklib::aesthetic_newick_tree( $tree );

        #----------------------------------------------------------------------
        #  Formulate the desired labels:
        #----------------------------------------------------------------------
        #  Build a function-to-color translation table based on frequency of
        #  function. Normally white is reserved for the current function, but
        #  there is none here. Assign colors until we run out, then go gray.
        #  Undefined function is not in %func_color, and so is not in
        #  %formatted_func
        #----------------------------------------------------------------------

        my %formatted_func = &FIGgjo::colorize_roles( $fid_func );

        #----------------------------------------------------------------------
        #  Tree tips are in alignment and tree uid.  We will fetch them in
        #  order from the tree, and build the label that we want.
        #----------------------------------------------------------------------

        my %labels;
        my $nbox = 0;
        my %boxes;
        foreach my $id ( gjonewicklib::newick_tip_list( $tree ) )
        {
            my $peg      = $fid_of_uid->{ $id };
            my $url      = HTML::fid_link( $cgi, $peg, 0, 1 );
            my $link     = $cgi->a( { -href => $url, -target => '_blank' }, $peg );
            my $reping   = representing( $cgi, $represents->{ $peg }, $org );
            my $func     = $fid_func->{ $peg };
            my $functext = $func ? $formatted_func{ $func } : '';
            my $orgname  = $org->{ $peg } ? html_esc( $org->{ $peg } ) : '';
            my $checkbox;
            if ( $can_assign && $orgname )
            {
                $checkbox = qq(<INPUT Type=checkbox Name=checked Value="$peg" OnClick="JavaScript:toggleTip(document.$form_name,this)" />&nbsp;);
                $boxes{ $id } = $nbox++;
            }
            my $proj_scr = $proj->{ $peg } ? $proj->{ $peg }->[2] : 0;
            my $hbar     = score_to_hbar( $proj_scr );
            my $subsyst  = $subsyst->{ $peg };
            my $cov_col  = $cov_color->{ $id };
            $cov_col     = '' if $cov_col eq '#FFFFFF';

            if ( $dlits->{$peg} && @{$dlits->{$peg}} )
            {
                $functext = qq(<SPAN Style='font-weight: bold;'>$functext</SPAN>) if $functext;
                $orgname  = qq(<SPAN Style='font-weight: bold;'>$orgname</SPAN>)  if $orgname;
            }

            my   @label;
            push @label, "$link&nbsp;";
            push @label, $reping                                               if $reping;
            push @label, qq(<SPAN Style="background-color: $cov_col;">)        if $cov_col;
            push @label, '&nbsp;';
            push @label, $checkbox                                             if $checkbox;
            push @label, qq(<INPUT Type=radio Name=from Value="$peg" />&nbsp;) if $func;
            push @label, "</SPAN>"                                             if $cov_col;
            push @label, $hbar;
            push @label, $subsyst                                              if $subsyst;
            push @label, $functext                                             if $functext;
            push @label, "[$orgname]"                                          if $orgname;
            # push @label, html_esc( $alias->{ $peg } )                          if $alias->{ $peg };

            $labels{ $id } = join( ' ', @label );
        }

        #------------------------------------------------------------------
        #  Draw the tree as printer plot.
        #------------------------------------------------------------------

        my $plot_options = { chars  => 'html',        # html-encoded unicode box set
                             format => 'chrlist_lbl', # line = [ \@symbs, $label ]
                             dy     =>  1,
                             min_dx =>  1,
                             width  => 64,
                             hash   =>  1             # get the layout info
                           };

        my ( $lines, $hash ) = gjonewicklib::text_plot_newick( $tree, $plot_options );

        my @nodeinfo = node_info( $tree, $hash, \%boxes );

        my @jsNodes;
        foreach my $info ( sort { $a->[1] <=> $b->[1] } @nodeinfo )
        {
            my ( $node, $nn, $x, $y, $desc, $box ) = @$info;
            my $line = $lines->[$y];

            #  Set internal node actions

            if ( @$desc )
            {
                my $chars = $line->[0];
                $chars->[$x] = qq(<A OnClick="JavaScript:setSubtree(document.$form_name,$nn,true)">$chars->[$x]</A>);
            }

            #  Put in final tip labels

            my $lbl = $line->[1];
            if ( defined( $lbl ) && length( $lbl ) && defined( $labels{$lbl} ) )
            {
                $line->[1] = $labels{$lbl};
            }

            push @jsNodes, "[ [ " . join( ', ', @$desc ) . "], $box ]";
        }

        my $jsTree = join( ",\n     ", @jsNodes );

        # Build the output lines

        my $btns = join ("\n", $cgi->br,
                               qq(<INPUT Type=button Name=Select Value="Select All"   OnClick="JavaScript:setAll(document.$form_name,true)" />),
                               qq(<INPUT Type=button Name=Select Value="Deselect All" OnClick="JavaScript:setAll(document.$form_name,false)" />),
                               qq(<INPUT Type=button Name=UndoBtn Value=Undo Disabled=disabled OnClick="JavaScript:undoSelect(document.$form_name)" />),
                               qq(<INPUT Type=button Name=RedoBtn Value=Redo Disabled=disabled OnClick="JavaScript:redoSelect(document.$form_name)" />),
                               $cgi->br,
                               ( $can_assign ? $cgi->submit( -name => 'request', -value => 'Assign' ) : () ),
                               $cgi->submit( -name => 'request', -value => 'Change focus peg' ),
                               $cgi->br,
                               ''
                        );

        push @$html, $btns;

        push @$html, qq(<DIV Class="tree">\n);
        foreach ( @$lines )
        {
            my $graphic = join( '', @{$_->[0]} );
            $graphic =~ s/((&nbsp;)+)/<SPAN Class=w>$1<\/SPAN>/g;
            $graphic =~ s/&nbsp;/&#9474;/g;
            my $lbl = defined( $_->[1] ) && length( $_->[1] ) ? " $_->[1]" : '';
            push @$html, "<PRE>$graphic</PRE>$lbl<BR />\n";
        }
        push @$html, "</DIV>\n";

        push @$html, $btns, <<"End_of_jsTree";
<SCRIPT Language="JavaScript">
//
//  The tree definition in form of descendent list and checkbox number.
//  Only one is used for each node.
//
document.$form_name.nodes
 = [ $jsTree
   ];
</SCRIPT>
End_of_jsTree

    }

    return @$html if wantarray;
}


#
#  Information about each tree node for building JavsScript tree definition,
#  and for putting links on internal tree nodes.
#
sub node_info
{
    my ( $node, $hash, $boxnum, $nodenum ) = @_;
    $nodenum ||= 0;
    my ( $x, $y ) = gjonewicklib::node_row_col( $node, $hash );
    my $desc = [];
    my $lbl  = gjonewicklib::newick_lbl( $node ) || '';
    my $box  = defined( $lbl ) && defined $boxnum->{$lbl} ? $boxnum->{$lbl} : 'null';
    my @nodes = ( [ $node, $nodenum, $x, $y, $desc, $box ] );
    foreach ( gjonewicklib::newick_desc_list( $node ) )
    {
        my $nextnode = $nodenum + @nodes;
        push @$desc, $nextnode;
        push @nodes, node_info( $_, $hash, $boxnum, $nextnode )
    }
    @nodes;
}


#==============================================================================
#  Select alignments and trees with given fid
#==============================================================================

sub show_alignments_and_trees_with_fid
{
    my ( $data ) = @_;
    my $html = $data->{ html } || [];
    my $sap  = $data->{ sap };
    my $cgi  = $data->{ cgi };
    my $fid  = $data->{ fid };

    if ( @{ $data->{ ali_tree_ids } } )
    {

        #  Get alignments and coverage data for the md5
        my $md5 = $data->{ md5 } ||= $data->{atree}->peg_to_md5($fid ) || '';
        my %alignID_coverages = map { $_->[0] => $_ }
                                $data->{tree}->alignment_coverages_of_md5( $md5 );
        my @ali_tree_ids = sort keys %alignID_coverages;

        #  Display the available alignments
        push @$html, $cgi->h2( "Select an Alignment and/or Tree" ),
                    '<TABLE>',
                    '<TR><TH Style="text-align:center">ID<BR />coverage</TH><TH>Count</TH><TH>Role</TH><TR>',
                    '<TABLEBODY>';
        foreach my $id ( @ali_tree_ids )
        {
            #  Work out the sequence coverage report:
            my $coverage = '';
            my ( undef, $md5, $len, $covers ) = @{ $alignID_coverages{$id} || [] };
            if ( $len && $covers && @$covers )
            {
                my @coverage;
                foreach ( @$covers )
                {
                    my ( $beg, $end ) = @$_[0,1];
                    my $color = region_color( $beg, $end, $len );
                    push @coverage, "<BR /><SPAN Style='background-color:$color;'>$beg-$end/$len</SPAN>";
                }
                $coverage = join( '', @coverage );
            }

            #  Compile the role data:
            my @role_data = $data->{atree}->roles_in_align($id );
            #  Limit to top 25 roles
            splice @role_data, 25 if @role_data > 25;

            #  Build the first table row for this alignment and tree:
            my $nrow = @role_data;
            my ( $role, $cnt ) = @{ shift @role_data };
            $role = html_esc( $role );
            push @$html, "<TR><TD RowSpan=$nrow Style='text-align:center'><INPUT Type=radio Name=ali_tree_id Value=$id />&nbsp;$id$coverage</TD>";
            push @$html, "    <TD Style='text-align:right'>$cnt</TD>";
            push @$html, "    <TD>$role</TD>";
            push @$html, "</TR>";

            #  Build the rest of the table rows for this alignment and tree:
            foreach ( @role_data )
            {
                ( $role, $cnt ) = @$_;
                $role = html_esc( $role );
                push @$html, "<TR>";
                push @$html, "    <TD Style='text-align:right'>$cnt</TD>";
                push @$html, "    <TD>$role</TD>";
                push @$html, "</TR>";
            }

            push @$html, '<TR><TD ColSpan=3><HR /></TD></TR>';
        }
        push @$html, '</TABLEBODY>',
                    '</TABLE>', $cgi->br,
                    $cgi->submit( -name => 'request', -value => 'Update' ),
                    $cgi->br;
    }
    elsif ( $fid )
    {
        @{ $data->{ ali_tree_ids } } = $data->{atree}->aligns_with_pegID( $fid );
        push @$html, "Sorry, no alignments with protein id '$fid'\n<BR /><BR />\n" if ! @{ $data->{ ali_tree_ids } };
    }

    return @$html if wantarray;
}


#-------------------------------------------------------------------------------
#  Vertical and horizontal histogram bar graphics characters:
#
#  9601-9608
#  9615-9608
#
#-------------------------------------------------------------------------------
sub score_to_vbar
{
    my ( $scr ) = @_;
    my $code = int(($scr ** 0.8) / 0.15) + 9601;
    return "&#$code;";
}


sub score_to_hbar
{
    my ( $scr ) = @_;
    return '&nbsp;&nbsp;' if $scr == 0;
    my $code = 9615 - int( 7.999 * ($scr ** 1.00) );
    return "&#$code;";
}



sub score_to_hbar2
{
    my ( $scr ) = @_;
    return '&nbsp;&nbsp;' if $scr == 0;
    $scr <= 0.5 ? "&#@{[9615-int(15.999*$scr)]};&nbsp;"
                : "&#9608;&#@{[9623-int(15.999*$scr)]};";
}


#-------------------------------------------------------------------------------
#  Background color as a function of sequence region covered by alignment
#-------------------------------------------------------------------------------

sub region_color
{
    my ( $b, $e, $n ) = @_;
    my ( $l, $r ) = ( $e > $b ) ? ( $b, $e ) : ( $e, $b );
    #  Allow 10 residues of split ends:
    my $split = 10;
    $l = $l          <= $split ?  1 : $l - $split;
    $r = $r + $split >= $n     ? $n : $r + $split;
    #  Base the hue on the position within the range of possible positions
    my $max = $n - ( $r - $l ) - 1;
    my $pos = ( $l - 1 - $max/2 ) / ( $max + 20 ) + 0.5;
    # my $pos = 0.5 * ( $l + $r - 2 ) / ( $n - 1 );
    my $hue = 5/6 * $pos - 1/12;
    my $sat = 1 - ($r-$l+1)/$n;
    my $br  = 1;
    &rgb2html( &hsb2rgb( $hue, $sat, $br ) );
}


sub hsb2rgb
{
    my ( $h, $s, $br ) = @_;
    $h = 6 * ($h - floor($h));
    if ( $s  > 1 ) { $s  = 1 } elsif ( $s  < 0 ) { $s  = 0 }
    if ( $br > 1 ) { $br = 1 } elsif ( $br < 0 ) { $br = 0 }
    my ( $r, $g, $b ) = ( $h <= 3 ) ? ( ( $h <= 1 ) ? ( 1,      $h,     0      )
                                      : ( $h <= 2 ) ? ( 2 - $h, 1,      0      )
                                      :               ( 0,      1,      $h - 2 )
                                      )
                                    : ( ( $h <= 4 ) ? ( 0,      4 - $h, 1      )
                                      : ( $h <= 5 ) ? ( $h - 4, 0,      1      )
                                      :               ( 1,      0,      6 - $h )
                                      );
    ( ( $r * $s + 1 - $s ) * $br,
      ( $g * $s + 1 - $s ) * $br,
      ( $b * $s + 1 - $s ) * $br
    )
}


sub rgb2html
{
    my ( $r, $g, $b ) = @_;
    if ( $r > 1 ) { $r = 1 } elsif ( $r < 0 ) { $r = 0 }
    if ( $g > 1 ) { $g = 1 } elsif ( $g < 0 ) { $g = 0 }
    if ( $b > 1 ) { $b = 1 } elsif ( $b < 0 ) { $b = 0 }
    sprintf("#%02x%02x%02x", int(255.999*$r), int(255.999*$g), int(255.999*$b) )
}


sub floor
{
    my $x = $_[0];
    defined( $x ) || return undef;
    ( ( $x >= 0 ) || ( int($x) == $x ) ) ? int( $x ) : -1 - int( - $x )
}


#-------------------------------------------------------------------------------
#  Show sequences represented by another sequence
#
#-------------------------------------------------------------------------------
sub representing
{
    my ( $cgi, $repsL, $orgsH ) = @_;
    return '' unless $cgi, $repsL && ref( $repsL ) eq 'ARRAY' && @$repsL;
    $orgsH ||= {};
    my @list = map  { qq(<P Style="text-indent: -2em; margin-left: 2em; margin-top: 0px; margin-bottom: 0px;">$_->[1] ($_->[0])) }
               sort { lc $a->[1] cmp lc $b->[1] }
               map  { [ $_, $orgsH->{ $_ } || '' ] }
               @$repsL;
    my $list = join( '', map { $_ } @list );
    text_with_tool_tip( $cgi, 'representing', 'Also representing', $list, "&#x271A;" );
}


#-------------------------------------------------------------------------------
#  Subsystems: with tooltip giving details
#
#    $html = subsystems( $fig, $cgi, $sap, $fid )
#-------------------------------------------------------------------------------

sub subsystems
{
    my ( $fig, $cgi, $sap, $fid ) = @_;
    return '' unless ( $fig || $sap ) && $cgi && $fid;

    my @in_sub = $fig ? $fig->peg_to_subsystems( $fid ) : ();
    return '' unless @in_sub;

    my $ss_list = join( $cgi->br, map { s/\_/ /g; $_ } sort { lc $a cmp lc $b } @in_sub );
    my $n_ss = @in_sub;
    text_with_tool_tip( $cgi, 'subsystems', 'Subsystems', $ss_list, "&#x2318;" );
}


sub text_with_tool_tip
{
    my ( $cgi, $id, $title, $tip, $link_text ) = @_;
    $cgi->a( { id          => $id,
               onMouseover => "javascript:if(!this.tooltip) this.tooltip=new Popup_Tooltip(this, '$title', '$tip', ''); this.tooltip.addHandler(); return false;"
             },
             $link_text
           )
}


sub help
{
    return ( <<'End_of_help_body', 'Help with Aligment and Tree Display' );
<H2>Help with Alignment and Tree Display</H2>
<H2>Controls</H2>
The details of the functions available depend upon the current status of the align/tree selection
process.  The options may include:

<UL>
<LI><B>SEED user</B> is available for entering a SEED user ID only when none is know.  In the
SeedViewer, this should be accomplished by using the login box in the upper-right corner of
the page.

<LI><B>Focus protein ID</B> can be used to identify a new focus protein for the current alignment
and tree (by clicking the <B>Update</B> button), or for selecting a new alignment and/or tree
containing the chosen protein (by clicking the <B>List all alignments and trees with this
protein</B> button).

<LI><B>Enter a SEED protein id</B> is used as a starting point for finding appropriate alignments
and trees when there in no current focus protein, alignment or tree.

<LI>The <B>Show alignment</B> and <B>Show tree</B> checkboxes control which data are displayed.
If neither is selected, a reminder message is displayed above the boxes.

<LI><B>Use for functions and assignments</B> determines the source and destination of the functions
displayed and the destination of requested assignment changes.  The ability to make assignments
depends on the identity of the user, and the selected data source.  If a given protein is not
present in the requested database, it is not displayed in the alignment or tree.

<LI><B>Of identical sequences, show</B> allows some control of the display of redundant sequences
in the alignment and tree.  By default, identical sequences with identical annotations are merged
into a single entry and the existence of hidden sequences are indicated by a &#x271A; following
the sequence identifier.  By selecting "all", the user can force the display of all sequences.
This can be critical for using the WWW browser Find function to locate a sequence from a particular
genome.

<LI><B>Color alignment by</B> allows control over the nature of the coloring shown in the displayed
sequence alignment.  The option only applies to the "default" alignment format.
</B> allows selection among the available alignment display formats.

<LI><B>Tree format</B> allows selection among the available tree display formats.

<LI>The <B>Update</B> button is used to refresh the display to reflect the current setting. Note
that this button does not change function assignments, even if sequences have been selected with
the checkboxes (only the <B>Assign</B> button makes assignments).

<LI>When present, the <B>Select an Alignment and/or Tree</B> heading is followed by a list of
alignments and trees for the current Focus protein.  The corresponding radio button is used to
select one of those available.
</UL>

<H2>Alignment</H2>
An alignment of the sequences, in one of the available formats.

<H2>Tree</H2>
The following descriptions apply to the "printer-plot" format of the tree.  Other formats
currently have fewer associated options.

<H3>Tree Action Buttons</H3>

<UL>
<LI><B>Select All</B> and <B>Deselect All</B> are for setting the state of all of the tree checkboxes.
Select All should be used with <B>extreme care</B> when changing a function assignment. Trees often
include related proteins with distinct functions, and multifunctional proteins are frequently
interspersed with proteins performing a single role (see <B>background color</b> below).

<LI><B>Undo</B> and <B>Redo</B> progressively undo (or redo) selections made on the tree.
Note that these cannot be used to undo function assignments that have already been made.

<LI>The <B>Assign</B> button is used to apply the function of the protein currently selected
by its radio button to all currently selected proteins. The success and failure status of the
function changes is reported at the top the page refreshes. The "attemped protein assignments
ignored" generally reflect proteins of identical sequence that are not in the SEED being accessed.

<LI>The <B>Change focus peg</B> button is used to update the protein from which function projection
scores are assessed.
</UL>

<H3>Tree Display</H3>
The tree is a printer-plot approximation of the computed tree.  It is arbitrarily rooted at the
approximate midpoint (though the definition used for midpoint can have multiple valid solutions,
so the rooting you get might not correspond to your intuition). Clicking on an internal node in
the tree selects all descendants of the node. Clicking on the checkbox of a tree tip selects the
individual sequence. The data associated with each tree tip can include the following:<BR />

<UL>
<LI>The <B>fig identifier</B> (fid) of the sequence.  It is a link to the corresponding protein in
the SEED.

<LI>&#x271A; is an optional element indicating that the given sequence is representing one or more
other proteins that are identical in sequence and have identical annotations (including any comments
that might be present).

<LI>A <B>check box</B> for selecting the given protein, most commonly for assigning a new function.
All identical sequences will be assigned the same new function (even if they are not selected, or
not even visible).  The check box is only present if the user is identified, and the data are from
a SEED that can be modified by the user.

<LI>A <B>radio button</B> for selecting either a protein with the function to be applied to the
checked sequences, or to define a new focus sequence for the projection data.  The radio button
is only present if the protein has an assigned function.

<LI>The <B>background color</B> of the check box and radio button reflect the amount of the given
protein that is included in the alignment used to compute the tree.  It is as though white light
was spread out in a spectrum from red to violet along the length of the protein, and the light could
pass through the portions of the sequence that are in the alignment, but not through those portions
that are not in the alignment.  The transmitted light is mixed back together and the resulting color
is shown.  Full length alignment from white light; alignment of only the beginning of the protein
makes red light; alignment of only the middle of the protein makes green light; alignment of only
the end of the protein makes blue light.  The shorter the region used in the alignment, the more
saturated the color.  If a protein is too colorful, care should be exercised in assigning a
function; it might be a multidomain (multifunctional) protein.

<LI>A <B>black bar</B>, when present, has a width proportional to the function projection score
from the the current focus protein to the protein in question.

<LI><B>&#x2318;</B> indicates that the protein is in a subsystem.  Mouse-over of the symbol
reveals the list of subsystems.

<LI>The <B>function currently assigned</B> to the protein.  Identical roles (components of the
function) are presented on the same background color.  The relative frequencies of the various
roles are indicated by a color progression (from most common or less common) that is roughly
brown, red, orange, yellow, green, blue, violet.  The least recurring functions are all on a
gray background.

<LI>The function is followed by (in square brackets) the <B>organism</B> from which the sequence was
obtained.

<LI>The function and organism are displayed in bold if the sequence is <B>directly linked to
literature</B> (a dlit) that provides insights into the function of the protein.
</UL>

End_of_help_body
}


1;
