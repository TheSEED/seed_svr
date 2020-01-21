package gjocolorlib;

#===============================================================================
#  Utilities for manipulating colors.
#
#  Based on component values from 0 to 1; undefined values are silently set
#  to 0.  Unless stated otherwise, RGB values are linear, and HTML color
#  values are sRGB.
#
#     @rgb  = ( $red, $green, $blue )             # gamma 1.0
#     @srgb = ( $red, $green, $blue )             # WWW std for pseudo gamma 2.2
#     @hsb  = ( $hue, $saturation, $brightness )  # hue = 0 is red
#     @hsy  = ( $hue, $saturation, $luma )        # luma is percieved brightness
#     @cmy  = ( $cyan, $magenta, $yellow )
#     @cmyk = ( $cyan, $magenta, $yellow, $black )
#     $html = '#xxxxxx', or named color
#
#  With very few exepctions, all exported functions take color components,
#  or a reference to an array of color components. Similarly, output is
#  an array of components, or a reference to an array of color components,
#  depending on context.
#
#  Exported functions:
#
#     @srgb = rgb2srgb( @rgb )
#     @rgb  = srgb2rgb( @srgb )
#
#     @rgb  = hsb2rgb( @hsb )
#     @rgb  = hsy2rgb( @hsy )
#     @hsb  = rgb2hsb( @rgb )
#     @hsy  = rgb2hsy( @rgb )
#
#     @rgb  = cmy2rgb( @cmy )
#     @rgb  = cmyk2rgb( @cmyk )
#     @cmy  = rgb2cym( @rgb )
#     @cmyk = rgb2cymk( @rgb )
#
#     $gray = rgb2gray(  @rgb )
#     @rgb  = gray2rgb( $gray )
#
#     @rgb  = html2rgb( $html )         #  name or hex sRGB to RGB
#     $html = rgb2html( @rgb )          #  linear RGB to sRGB html
#     $html = rgb2html_g10( @rgb )      #  gamma = 1.0, or @rgb is sRGB
#     $html = rgb2html_g18( @rgb )      #  gamma = 1.8
#     $html = rgb2html_g22( @rgb )      #  gamma = 2.2
#     $html = gray2html( $gray )        #  linear gray to sRGB html
#     $html = gray2html_g10( $gray )    #  gamma = 1.0
#     $html = gray2html_g18( $gray )    #  gamma = 1.8
#     $html = gray2html_g22( $gray )    #  gamma = 2.2
#     $name = rgb2name( @rgb )          #  CSS 3.0 and SVG 1.0 colors
#
#     @rgb  = blend_rgb_colors( \@color1, \@color2, ... )
#     $html = blend_html_colors( $html1, $html2, ... )
#
#  Internal functions for validated input values (no reference as input):
#
#     @srgb = linear2srgb( @rgb )
#     @rgb  = srgb2linear( @srgb )
#     @rgb  = hs2rgb( $hue, $saturation )  # brighness = 1
#     @rgb  = rgb2hsb0( @rgb )
#     $gray = rgb2gray0( @rgb )
#    \@rgb  = htmlhex2rgb( $html )      # returns reference due to its usage
#     $html = srgb2html0( @srgb )
#     $dist = rgb_distance( \@rgb1, \@rgb2 )
#
#  Internal data include:
#
#     @name2rgb      #  Array of CamalCase standard color names and RGB values
#                    #      used to match name to given RGB value
#     %lc_name2rgb   #  Hash from lowercase name to RGB values (includes some
#                    #      nonstandard names)
#
#===============================================================================

use strict;

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(
        rgb2srgb
        srgb2rgb
        hsb2rgb
        hsy2rgb
        cmy2rgb
        cmyk2rgb
        rgb2gray
        gray2rgb
        html2rgb
        rgb2html
        rgb2html_g10
        rgb2html_g18
        rgb2html_g22
        gray2html
        gray2html_g10
        gray2html_g18
        gray2html_g22
        rgb2name
        blend_rgb_colors
        blend_html_colors
        );
our @EXPORT_OK = qw(
        UI_Orange
        UI_Formal_Orange
        UI_Blue
        UI_Formal_Blue
        );


my $UI_Blue          = [ 125/255,  60/255, 125/255 ];
my $UI_Orange        = [ 244/255, 127/255,  36/255 ];
my $UI_Formal_Blue   = [ 110/255, 139/255, 191/255 ];
my $UI_Formal_Orange = [ 239/255, 138/255,  28/255 ];

#  Use floor to get proper modulo 1 values

sub floor
{
    my $x = $_[0] || 0;
    ( $x >= 0 ) || ( int($x) == $x ) ? int( $x ) : -1 - int( - $x )
}

sub min { ( $_[0] < $_[1] ) ? $_[0] : $_[1] }
sub max { ( $_[0] > $_[1] ) ? $_[0] : $_[1] }

#
#  Many values need to be clipped to a range from 0 to 1.  The first function
#  returns a value in the range.  It creates the value if necessary.  The
#  second function adjusts args in place.  It cannot create a value.
#
#  $zero2one = zero2one( $val )
#

sub zero2one { local $_ = $_[0]; ( ! defined || ( $_ <= 0 ) ) ? 0 : $_ > 1 ? 1 : $_ }

#
#  make_zero2one( $a, $b, ... )
#

sub make_zero2one   # in place
{
    foreach ( @_ )
    {
        if ( ! defined || $_ <= 0 ) { $_ = 0 } elsif ( $_ > 1 ) { $_ = 1 }
    }
}

#===============================================================================
#  Conversions between linear RGB and sRGB, the standard for default color
#  values.
#-------------------------------------------------------------------------------
#  The breakpoints for the conversions are different in the specification at
#  www.w3.org and the version on wikipedia.org.  The latter actually matches
#  values to more decimal places and is used here.
#
#      http://www.w3.org/Graphics/Color/sRGB.html
#      http://en.wikipedia.org/wiki/SRGB
#
#-------------------------------------------------------------------------------
#  Internal functions for doing the forward and reverse transforms from
#  validated triplets:
#
#    @sRGB = linear2srgb(  @RGB )
#     @RGB = srgb2linear( @sRGB )
#
#-------------------------------------------------------------------------------
sub linear2srgb
{
  # map { $_ <= 0.00314   ? 12.92*$_ : 1.055*($_**(1/2.4))-0.055 } @_   # w3.org spec
    map { $_ <= 0.0031308 ? 12.92*$_ : 1.055*($_**(1/2.4))-0.055 } @_   # wiki spec
}


sub srgb2linear
{
  # map { $_ <= 0.03928 ? $_/12.92 : (($_+0.055)/1.055)**2.4 } @_   # w3.org spec
    map { $_ <= 0.04045 ? $_/12.92 : (($_+0.055)/1.055)**2.4 } @_   # wiki spec
}


#-------------------------------------------------------------------------------
#  External functions for doing the forward and reverse transforms:
#
#   \@sRGB = rgb2srgb( \@RGB )
#   \@sRGB = rgb2srgb(  @RGB )
#    @sRGB = rgb2srgb( \@RGB )
#    @sRGB = rgb2srgb(  @RGB )
#
#    \@RGB = srgb2rgb( \@sRGB )
#    \@RGB = srgb2rgb(  @sRGB )
#     @RGB = srgb2rgb( \@sRGB )
#     @RGB = srgb2rgb(  @sRGB )
#
#-------------------------------------------------------------------------------
sub rgb2srgb
{
    my ( $r, $g, $b ) = $_[0] && ref( $_[0] ) eq 'ARRAY' ? @{$_[0]} : @_;
    make_zero2one( $r, $g, $b );
    my @srgb = linear2srgb( $r, $g, $b );
    wantarray ? @srgb : \@srgb;
}


sub srgb2rgb
{
    my ( $r, $g, $b ) = $_[0] && ref( $_[0] ) eq 'ARRAY' ? @{$_[0]} : @_;
    make_zero2one( $r, $g, $b );
    my @rgb = srgb2linear( $r, $g, $b );
    wantarray ? @rgb : \@rgb;
}


#-------------------------------------------------------------------------------
#  Look-up table for 0 - 255 index of sRGB (i.e., html hex values) to
#  linear RGB value in 0-1 range.
#-------------------------------------------------------------------------------

my @index2rgb = map { srgb2linear( $_/255 ) } ( 0 .. 255 );


#===============================================================================
#  Interconversions of RGB and HSB and HSY values.  HSY is same idea as HSB,
#  but with luma (the visual lightness, y, range 0 - 1).  Colors are darkened
#  or desaturated as necessary to match the luma value.
#
#     @rgb = hsb2rgb( @hsb )
#     @rgb = hsy2rgb( @hsy )
#     @hsb = rgb2hsb( @rgb )
#     @hsy = rgb2hsy( @rgb )
#
#-------------------------------------------------------------------------------
#
#  Internal function to get RGB from validated hue and saturation:
#

sub hs2rgb
{
    my ( $h, $s ) = @_;
    my $h6 = 6 * ( $h - floor($h) );   #  Hue is cyclic modulo 1
    my $m  = 1 - $s;
    map { $_ * $s + $m } ( $h6 <= 3 ) ? ( ( $h6 <= 1 ) ? (     1, $h6,     0 )
                                        : ( $h6 <= 2 ) ? ( 2-$h6,   1,     0 )
                                        :                (     0,   1, $h6-2 )
                                        )
                                      : ( ( $h6 <= 4 ) ? (     0, 4-$h6,     1 )
                                        : ( $h6 <= 5 ) ? ( $h6-4,     0,     1 )
                                        :                (     1,     0, 6-$h6 )
                                        );
}


sub hsb2rgb
{
    my ( $h, $s, $br ) = $_[0] && ref( $_[0]) eq 'ARRAY' ? @{$_[0]} : @_;
    make_zero2one( $h, $s, $br );
    my @rgb = map { $_ * $br } hs2rgb( $h, $s );
    wantarray() ? @rgb : \@rgb;
}


sub hsy2rgb
{
    my ( $h, $s, $y ) = $_[0] && ref( $_[0] ) eq 'ARRAY' ? @{$_[0]} : @_;
    make_zero2one( $h, $s, $y );

    my ( $r, $g, $b ) = hs2rgb( $h, $s );
    my $luma = rgb2gray0( $r, $g, $b );     #  How bright is the color?

    my @rgb;
    if ( $luma <  $y )          #  Too dim without decreasing saturation
    {
        my $s = ( 1 - $y ) / ( 1 - $luma ) ;
        my $m = 1 - $s;
        @rgb = map { $_ * $s + $m } ( $r, $g, $b );
    }
    else                        #  Too bright
    {
        my $k = $y / $luma;
        @rgb = map { $_ * $k } ( $r, $g, $b );
    }

    wantarray() ? @rgb : \@rgb;
}


#
#  Internal function for finding HSB from validated RGB
#

sub rgb2hsb0
{
    my ( $min, undef, $br ) = sort { $a <=> $b } @_;
    if ( $br == $min ) { return wantarray ? ( 0, 0, $br ) : [ 0, 0, $br ] }

    my ( $r, $g, $b ) = @_;
    my $s_br = $br - $min;    #  $s * $br
    my $h6;

    if    ( $r == $br ) { $h6 = ( $b > $g ) ? 6 - ($b-$g)/$s_br : ($g-$b)/$s_br }
    elsif ( $g == $br ) { $h6 = 2 + ($b-$r)/$s_br }
    else                { $h6 = 4 + ($r-$g)/$s_br }

    ( $h6/6, $s_br/$br, $br );
}


sub rgb2hsb
{
    my ( $r, $g, $b ) = $_[0] && ref( $_[0]) eq 'ARRAY' ? @{$_[0]} : @_;
    make_zero2one( $r, $g, $b );

    my @hsb = rgb2hsb0( $r, $g, $b );

    wantarray ? @hsb : \@hsb;
}


sub rgb2hsy
{
    my ( $r, $g, $b ) = $_[0] && ref( $_[0]) eq 'ARRAY' ? @{$_[0]} : @_;
    make_zero2one( $r, $g, $b );

    my @hsy = rgb2hsb0( $r, $g, $b );
    $hsy[2] = rgb2gray0( $r, $g, $b );

    wantarray() ? @hsy : \@hsy;
}


#===============================================================================
#  Conversions between linear RGB and CMY and CMYK values.
#
#     @rgb  = cmy2rgb( @cmy )
#     @rgb  = cmyk2rgb( @cmyk )
#     @cmy  = rgb2cym( @rgb )
#     @cmyk = rgb2cymk( @rgb )
#
#-------------------------------------------------------------------------------

sub cmy2rgb
{
    my ( $c, $m, $y ) = $_[0] && ref( $_[0]) eq 'ARRAY' ? @{$_[0]} : @_;
    make_zero2one( $c, $m, $y );
    my @rgb = ( 1 - $c, 1 - $m, 1 - $y );
    wantarray() ? @rgb : \@rgb;
}


sub cmyk2rgb
{
    my ( $c, $m, $y, $k ) = $_[0] && ref( $_[0]) eq 'ARRAY' ? @{$_[0]} : @_;
    make_zero2one( $c, $m, $y, $k );
    my $br = 1 - $k;
    my @rgb = ( ( 1 - $c ) * $br, ( 1 - $m ) * $br, ( 1 - $y ) * $br );
    wantarray() ? @rgb : \@rgb;
}


sub rgb2cmy
{
    my ( $r, $g, $b ) = $_[0] && ref( $_[0]) eq 'ARRAY' ? @{$_[0]} : @_;
    make_zero2one( $r, $g, $b );
    my @cmy = ( 1 - $r, 1 - $g, 1 - $b );
    wantarray() ? @cmy : \@cmy;
}


sub rgb2cmyk
{
    my ( $r, $g, $b ) = $_[0] && ref( $_[0]) eq 'ARRAY' ? @{$_[0]} : @_;
    make_zero2one( $r, $g, $b );
    my $br = $r > $g ? ( $r > $b ? $r : $b ) : ( $g > $b ? $g : $b );
    my @cmyk = $br > 0 ? ( 1 - $r/$br, 1 - $g/$br, 1 - $b/$br, 1 - $br )
                       : ( 0, 0, 0, 1 );
    wantarray() ? @cmyk : \@cmyk;
}


#===============================================================================
#  Conversions between linear RGB and gray (luma) values.
#-------------------------------------------------------------------------------
#
#    Various ITU recommendations and the new standard:
#        0.299 * $r + 0.587 * $g + 0.114 * $b;  #  Rec 601-1 (NTSC)
#        0.213 * $r + 0.715 * $g + 0.072 * $b;  #  Rec 709   (sRGB)
#        0.222 * $r + 0.707 * $g + 0.071 * $b;  #  ITU std (D65 white point)
#        0.330 * $r + 0.590 * $g + 0.080 * $b;  #  GJO
#
#-------------------------------------------------------------------------------
#
#   Internal function for gray (=luma) from validated RGB:
#
sub rgb2gray0
{
    my ( $r, $g, $b ) = @_;

    0.213 * $r + 0.715 * $g + 0.072 * $b;  #  Rec 709   (sRGB)
}


#-------------------------------------------------------------------------------
#  Produce the gray equivalent in brightness to an RGB value:
#
#     $gray = rgb2gray(  @rgb )
#     @rgb  = gray2rgb( $gray )
#
#-------------------------------------------------------------------------------

sub rgb2gray
{
    my ( $r, $g, $b ) = $_[0] && ref($_[0]) eq 'ARRAY' ? @{$_[0]} : @_;
    make_zero2one( $r, $g, $b );
    rgb2gray0( $r, $g, $b );
}


sub gray2rgb
{
    my $gray = zero2one( $_[0] );
    my @rgb = ( $gray, $gray, $gray );
    wantarray ? @rgb : \@rgb;
}


#===============================================================================
#  Conversions amongst RGB values, HTML hex code values, and color names.
#-------------------------------------------------------------------------------
#
#  Named colors in HTML and SVG specifications.
#  Beware that this copy has mixed case names.
#
my %name2html =
      ( AliceBlue            => '#F0F8FF',  # SVG 1.0
        AntiqueWhite         => '#FAEBD7',  # SVG 1.0
        Aqua                 => '#00FFFF',  # CSS 3.0
        Aquamarine           => '#7FFFD4',  # SVG 1.0
        Azure                => '#F0FFFF',  # SVG 1.0
        Beige                => '#F5F5DC',  # SVG 1.0
        Bisque               => '#FFE4C4',  # SVG 1.0
        Black                => '#000000',  # CSS 3.0
        BlanchedAlmond       => '#FFEBCD',  # SVG 1.0
        Blue                 => '#0000FF',  # CSS 3.0
        BlueViolet           => '#8A2BE2',  # SVG 1.0
        Brown                => '#A52A2A',  # SVG 1.0
        Burlywood            => '#DEB887',  # SVG 1.0
        CadetBlue            => '#5F9EA0',  # SVG 1.0
        Chartreuse           => '#7FFF00',  # SVG 1.0
        Chocolate            => '#D2691E',  # SVG 1.0
        Coral                => '#FF7F50',  # SVG 1.0
        CornflowerBlue       => '#6495ED',  # SVG 1.0
        Cornsilk             => '#FFF8DC',  # SVG 1.0
        Crimson              => '#DC143C',  # SVG 1.0
        Cyan                 => '#00FFFF',  # SVG 1.0 = Teal
        DarkBlue             => '#00008B',  # SVG 1.0
        DarkCyan             => '#008B8B',  # SVG 1.0
        DarkGoldenrod        => '#B8860B',  # SVG 1.0
        DarkGray             => '#A9A9A9',  # SVG 1.0
        DarkGreen            => '#006400',  # SVG 1.0
        DarkKhaki            => '#BDB76B',  # SVG 1.0
        DarkMagenta          => '#8B008B',  # SVG 1.0
        DarkOliveGreen       => '#556B2F',  # SVG 1.0
        Darkorange           => '#FF8C00',  # SVG 1.0
        DarkOrchid           => '#9932CC',  # SVG 1.0
        DarkRed              => '#8B0000',  # SVG 1.0
        DarkSalmon           => '#E9967A',  # SVG 1.0
        DarkSeaGreen         => '#8FBC8F',  # SVG 1.0
        DarkSlateBlue        => '#483D8B',  # SVG 1.0
        DarkSlateGray        => '#2F4F4F',  # SVG 1.0
        DarkTurquoise        => '#00CED1',  # SVG 1.0
        DarkViolet           => '#9400D3',  # SVG 1.0
        DeepPink             => '#FF1493',  # SVG 1.0
        DeepSkyBlue          => '#00BFFF',  # SVG 1.0
        DimGray              => '#696969',  # SVG 1.0
        Dimgrey              => '#696969',  # SVG 1.0
        DodgerBlue           => '#1E90FF',  # SVG 1.0
        FireBrick            => '#B22222',  # SVG 1.0
        FloralWhite          => '#FFFAF0',  # SVG 1.0
        ForestGreen          => '#228B22',  # SVG 1.0
        Fuchsia              => '#FF00FF',  # CSS 3.0 = Magenta
        Gainsboro            => '#DCDCDC',  # SVG 1.0
        GhostWhite           => '#F8F8FF',  # SVG 1.0
        Gold                 => '#FFD700',  # SVG 1.0
        Goldenrod            => '#DAA520',  # SVG 1.0
        Gray                 => '#808080',  # CSS 3.0
        Green                => '#008000',  # CSS 3.0
        GreenYellow          => '#ADFF2F',  # SVG 1.0
        Grey                 => '#808080',  # SVG 1.0
        Honeydew             => '#F0FFF0',  # SVG 1.0
        HotPink              => '#FF69B4',  # SVG 1.0
        IndianRed            => '#CD5C5C',  # SVG 1.0
        Indigo               => '#4B0082',  # SVG 1.0
        Ivory                => '#FFFFF0',  # SVG 1.0
        Khaki                => '#F0E68C',  # SVG 1.0
        Lavender             => '#E6E6FA',  # SVG 1.0
        LavenderBlush        => '#FFF0F5',  # SVG 1.0
        LawnGreen            => '#7CFC00',  # SVG 1.0
        LemonChiffon         => '#FFFACD',  # SVG 1.0
        LightBlue            => '#ADD8E6',  # SVG 1.0
        LightCoral           => '#F08080',  # SVG 1.0
        LightCyan            => '#E0FFFF',  # SVG 1.0
        LightGoldenrodYellow => '#FAFAD2',  # SVG 1.0
        LightGray            => '#D3D3D3',  # SVG 1.0
        LightGreen           => '#90EE90',  # SVG 1.0
        LightPink            => '#FFB6C1',  # SVG 1.0
        LightSalmon          => '#FFA07A',  # SVG 1.0
        LightSeaGreen        => '#20B2AA',  # SVG 1.0
        LightSkyBlue         => '#87CEFA',  # SVG 1.0
        LightSlateGray       => '#778899',  # SVG 1.0
        LightSteelBlue       => '#B0C4DE',  # SVG 1.0
        LightYellow          => '#FFFFE0',  # SVG 1.0
        Lime                 => '#00FF00',  # CSS 3.0
        LimeGreen            => '#32CD32',  # SVG 1.0
        Linen                => '#FAF0E6',  # SVG 1.0
        Magenta              => '#FF00FF',  # SVG 1.0 = Fuchsia
        Maroon               => '#800000',  # CSS 3.0
        MediumAquamarine     => '#66CDAA',  # SVG 1.0
        MediumBlue           => '#0000CD',  # SVG 1.0
        MediumOrchid         => '#BA55D3',  # SVG 1.0
        MediumPurple         => '#9370DB',  # SVG 1.0
        MediumSeaGreen       => '#3CB371',  # SVG 1.0
        MediumSlateBlue      => '#7B68EE',  # SVG 1.0
        MediumSpringGreen    => '#00FA9A',  # SVG 1.0
        MediumTurquoise      => '#48D1CC',  # SVG 1.0
        MediumVioletRed      => '#C71585',  # SVG 1.0
        MidnightBlue         => '#191970',  # SVG 1.0
        MintCream            => '#F5FFFA',  # SVG 1.0
        MistyRose            => '#FFE4E1',  # SVG 1.0
        Moccasin             => '#FFE4B5',  # SVG 1.0
        NavajoWhite          => '#FFDEAD',  # SVG 1.0
        Navy                 => '#000080',  # CSS 3.0
        OldLace              => '#FDF5E6',  # SVG 1.0
        Olive                => '#808000',  # CSS 3.0
        OliveDrab            => '#6B8E23',  # SVG 1.0
        Orange               => '#FFA500',  # SVG 1.0
        OrangeRed            => '#FF4500',  # SVG 1.0
        Orchid               => '#DA70D6',  # SVG 1.0
        PaleGoldenrod        => '#EEE8AA',  # SVG 1.0
        PaleGreen            => '#98FB98',  # SVG 1.0
        PaleTurquoise        => '#AFEEEE',  # SVG 1.0
        PaleVioletRed        => '#DB7093',  # SVG 1.0
        PapayaWhip           => '#FFEFD5',  # SVG 1.0
        PeachPuff            => '#FFDAB9',  # SVG 1.0
        Peru                 => '#CD853F',  # SVG 1.0
        Pink                 => '#FFC0CB',  # SVG 1.0
        Plum                 => '#DDA0DD',  # SVG 1.0
        PowderBlue           => '#B0E0E6',  # SVG 1.0
        Purple               => '#800080',  # CSS 3.0
        Red                  => '#FF0000',  # CSS 3.0
        RosyBrown            => '#BC8F8F',  # SVG 1.0
        RoyalBlue            => '#4169E1',  # SVG 1.0
        SaddleBrown          => '#8B4513',  # SVG 1.0
        Salmon               => '#FA8072',  # SVG 1.0
        SandyBrown           => '#F4A460',  # SVG 1.0
        SeaGreen             => '#2E8B57',  # SVG 1.0
        Seashell             => '#FFF5EE',  # SVG 1.0
        Sienna               => '#A0522D',  # SVG 1.0
        Silver               => '#C0C0C0',  # CSS 3.0
        SkyBlue              => '#87CEEB',  # SVG 1.0
        SlateBlue            => '#6A5ACD',  # SVG 1.0
        SlateGray            => '#708090',  # SVG 1.0
        Snow                 => '#FFFAFA',  # SVG 1.0
        SpringGreen          => '#00FF7F',  # SVG 1.0
        SteelBlue            => '#4682B4',  # SVG 1.0
        Tan                  => '#D2B48C',  # SVG 1.0
        Teal                 => '#008080',  # CSS 3.0 = Cyan
        Teal                 => '#008080',  # SVG 1.0
        Thistle              => '#D8BFD8',  # SVG 1.0
        Tomato               => '#FF6347',  # SVG 1.0
        Turquoise            => '#40E0D0',  # SVG 1.0
        Violet               => '#EE82EE',  # SVG 1.0
        Wheat                => '#F5DEB3',  # SVG 1.0
        White                => '#FFFFFF',  # CSS 3.0
        WhiteSmoke           => '#F5F5F5',  # SVG 1.0
        Yellow               => '#FFFF00',  # CSS 3.0
        YellowGreen          => '#9ACD32',  # SVG 1.0
    );

#
#  Additional color names for grays and anything else we might want to
#  recognize.
#

my %name2html_too =
      ( Gray0                => '#000000',
        Gray1                => '#030303',
        Gray2                => '#050505',
        Gray3                => '#080808',
        Gray4                => '#0a0a0a',
        Gray5                => '#0d0d0d',
        Gray6                => '#0f0f0f',
        Gray7                => '#121212',
        Gray8                => '#141414',
        Gray9                => '#171717',
        Gray10               => '#1a1a1a',
        Gray11               => '#1c1c1c',
        Gray12               => '#1f1f1f',
        Gray13               => '#212121',
        Gray14               => '#242424',
        Gray15               => '#262626',
        Gray16               => '#292929',
        Gray17               => '#2b2b2b',
        Gray18               => '#2e2e2e',
        Gray19               => '#303030',
        Gray20               => '#333333',
        Gray21               => '#363636',
        Gray22               => '#383838',
        Gray23               => '#3b3b3b',
        Gray24               => '#3d3d3d',
        Gray25               => '#404040',
        Gray26               => '#424242',
        Gray27               => '#454545',
        Gray28               => '#474747',
        Gray29               => '#4a4a4a',
        Gray30               => '#4d4d4d',
        Gray31               => '#4f4f4f',
        Gray32               => '#525252',
        Gray33               => '#545454',
        Gray34               => '#575757',
        Gray35               => '#595959',
        Gray36               => '#5c5c5c',
        Gray37               => '#5e5e5e',
        Gray38               => '#616161',
        Gray39               => '#636363',
        Gray40               => '#666666',
        Gray41               => '#696969',
        Gray42               => '#6b6b6b',
        Gray43               => '#6e6e6e',
        Gray44               => '#707070',
        Gray45               => '#737373',
        Gray46               => '#757575',
        Gray47               => '#787878',
        Gray48               => '#7a7a7a',
        Gray49               => '#7d7d7d',
        Gray50               => '#7f7f7f',
        Gray51               => '#828282',
        Gray52               => '#858585',
        Gray53               => '#878787',
        Gray54               => '#8a8a8a',
        Gray55               => '#8c8c8c',
        Gray56               => '#8f8f8f',
        Gray57               => '#919191',
        Gray58               => '#949494',
        Gray59               => '#969696',
        Gray60               => '#999999',
        Gray61               => '#9c9c9c',
        Gray62               => '#9e9e9e',
        Gray63               => '#a1a1a1',
        Gray64               => '#a3a3a3',
        Gray65               => '#a6a6a6',
        Gray66               => '#a8a8a8',
        Gray67               => '#ababab',
        Gray68               => '#adadad',
        Gray69               => '#b0b0b0',
        Gray70               => '#b3b3b3',
        Gray71               => '#b5b5b5',
        Gray72               => '#b8b8b8',
        Gray73               => '#bababa',
        Gray74               => '#bdbdbd',
        Gray75               => '#bfbfbf',
        Gray76               => '#c2c2c2',
        Gray77               => '#c4c4c4',
        Gray78               => '#c7c7c7',
        Gray79               => '#c9c9c9',
        Gray80               => '#cccccc',
        Gray81               => '#cfcfcf',
        Gray82               => '#d1d1d1',
        Gray83               => '#d4d4d4',
        Gray84               => '#d6d6d6',
        Gray85               => '#d9d9d9',
        Gray86               => '#dbdbdb',
        Gray87               => '#dedede',
        Gray88               => '#e0e0e0',
        Gray89               => '#e3e3e3',
        Gray90               => '#e5e5e5',
        Gray91               => '#e8e8e8',
        Gray92               => '#ebebeb',
        Gray93               => '#ededed',
        Gray94               => '#f0f0f0',
        Gray95               => '#f2f2f2',
        Gray96               => '#f5f5f5',
        Gray97               => '#f7f7f7',
        Gray98               => '#fafafa',
        Gray99               => '#fcfcfc',
        Gray100              => '#ffffff',

        #  And just for fun:

        IllinoisBlue         => '#003C7D',
        IllinoisFormalBlue   => '#6EBBBF',
        IllinoisFormalOrange => '#EF8A1C',
        IllinoisOrange       => '#F47F24',
        UIBlue               => '#003C7D',
        UIFormalBlue         => '#6EBBBF',
        UIFormalOrange       => '#EF8A1C',
        UIOrange             => '#F47F24',
    );

#  Version with pretty names, used for RGB to name:

my @name2rgb = map { [ $_ => htmlhex2rgb( $name2html{ $_ } ) ] } keys %name2html;

#  Version with lowercase names, used for name to RGB (includes nonstandard
#  names):

my %lc_name2rgb = ( ( map { lc( $_->[0] ) => $_->[1] } @name2rgb ),
                    ( map { lc( $_ ) => htmlhex2rgb( $name2html_too{ $_ } ) } keys %name2html_too )
                  );

#-------------------------------------------------------------------------------
#  Convert html to RGB and sRGB:
#
#    @rgb  = html2rgb( $htmlcolor )
#
#  where $htmlcolor is:
#
#      namedcolor or #xxxxxx or #xxx or xxxxxx or xxx
#
#-------------------------------------------------------------------------------
#
#  Internal function for converting html sRGB hex strings to RGB values.
#  It will take 3 or 6 hexadecimal digits, with or without a leading #.
#
sub htmlhex2rgb
{
    local $_ = $_[0] || '';
    my @hex_rgb = m/^#?([\da-f][\da-f])([\da-f][\da-f])([\da-f][\da-f])$/i ? ( $1, $2, $3 )
                : m/^#?([\da-f])([\da-f])([\da-f])$/i                      ? map { "$_$_" } ( $1, $2, $3 )
                :                                                            ( '00', '00', '00' );

    [ map { $index2rgb[ hex( $_ ) ] } @hex_rgb ];
}


sub html2rgb
{
    my $html = lc ( shift || '' );     #  Only lower case is indexed
    $html =~ s/\s+//g;                 #  No spaces
    $html =~ s/grey/gray/g;            #  Only USA spelling is indexed

    my $rgb = $lc_name2rgb{ $html } || htmlhex2rgb( $html );

    wantarray ? @$rgb : $rgb;
}


#-------------------------------------------------------------------------------
#  Convert an RGB and gray values to sRGB HTML strings
#
#     $html = rgb2html( @rgb )          #  linear RGB to sRGB html
#     $html = rgb2html_g10( @rgb )      #  gamma = 1.0, or @rgb is sRGB
#     $html = rgb2html_g18( @rgb )      #  gamma = 1.8
#     $html = rgb2html_g22( @rgb )      #  gamma = 2.2
#     $html = gray2html( $gray )        #  linear gray to sRGB html
#     $html = gray2html_g10( $gray )    #  gamma = 1.0
#     $html = gray2html_g18( $gray )    #  gamma = 1.8
#     $html = gray2html_g22( $gray )    #  gamma = 2.2
#
#-------------------------------------------------------------------------------
#
#  Internal routine to format string:
#
sub srgb2html0
{
    sprintf( '#%02x%02x%02x', map { int( 255.999 * $_ ) } @_ );
}


#-------------------------------------------------------------------------------
#  Convert an RGB value to an HTML string:
#-------------------------------------------------------------------------------
sub rgb2html
{
    my ( $r, $g, $b ) = defined($_[0]) && ref($_[0]) eq 'ARRAY' ? @{$_[0]} : @_;
    make_zero2one( $r, $b, $b );
    srgb2html0( linear2srgb( $r, $g, $b ) );
}

#
#  RGB to HTML, with gamma adjustment.  This should normally be done using
#  the sRGB correction invoked by rgb2html().  g10 (gamma 1.0) writes
#  unadjusted RGB values, so it would be appropriate for:
#
#     $html = rgb2html_g10( @sRGB ).
#
#  gamma = 1.0
#
sub rgb2html_g10
{
    my ( $r, $g, $b ) = defined($_[0]) && ref($_[0]) eq 'ARRAY' ? @{$_[0]} : @_;
    srgb2html0( map { zero2one($_) } ( $r, $g, $b ) );
}

#
#  gamma = 1.8
#
sub rgb2html_g18
{
    my ( $r, $g, $b ) = defined($_[0]) && ref($_[0]) eq 'ARRAY' ? @{$_[0]} : @_;
    srgb2html0( map { zero2one($_) ** (1/1.8) } ( $r, $g, $b ) );
}

#
#  gamma = 2.2
#
sub rgb2html_g22
{
    my ( $r, $g, $b ) = defined($_[0]) && ref($_[0]) eq 'ARRAY' ? @{$_[0]} : @_;
    srgb2html0( map { zero2one($_) ** (1/2.2) } ( $r, $g, $b ) );
}


#-------------------------------------------------------------------------------
#  Convert an gray value to an HTML string.
#-------------------------------------------------------------------------------
sub gray2html
{
    my $gray = linear2srgb( zero2one( $_[0] ) );
    srgb2html0( $gray, $gray, $gray );
}


#
#  gray to HTML, with gamma adjustment.  This should normally be done using
#  the sRGB correction invoked by gray2html().
#
#  gamma = 1.0
#
sub gray2html_g10
{
    my $gray = zero2one( $_[0] );
    srgb2html0( $gray, $gray, $gray );
}

#
#  gamma = 1.8
#
sub gray2html_g18
{
    my $gray = zero2one( $_[0] ) ** 1.8;
    srgb2html0( $gray, $gray, $gray );
}

#
#  gamma = 2.2
#
sub gray2html_g22
{
    my $gray = zero2one( $_[0] ) ** 2.2;
    srgb2html0( $gray, $gray, $gray );
}


#-------------------------------------------------------------------------------
#  Find the closest named color from CSS 3.0 and SVG 1.0:
#
#    $name = rgb2name( @rgb )
#
#-------------------------------------------------------------------------------
sub rgb2name
{
    my ( $r, $g, $b ) = ref( $_[0] ) eq 'ARRAY' ? @{ $_[0] } : @_;
    make_zero2one( $r, $g, $b );

    my $rgb = [ $r, $g, $b ];
    my ( $name, $dmin ) = ( '', 3 );
    foreach ( @name2rgb )
    {
        my $d = rgb_distance( $rgb, $_->[1] );
        if ( $d < $dmin ) { $name = $_->[0]; $dmin = $d; }
    }

    $name;
}


#
#  Internal function defining a distance between 2 rgb colors
#
#    $distance = rgb_distance( $rgb1, $rgb2 )
#
sub rgb_distance
{
    my ( $r1, $g1, $b1 ) = $_[0] && ( ref( $_[0] ) eq 'ARRAY' ) ? @{$_[0]} : (0,0,0);
    my ( $r2, $g2, $b2 ) = $_[1] && ( ref( $_[1] ) eq 'ARRAY' ) ? @{$_[1]} : (0,0,0);

    make_zero2one( $r1, $g1, $b1 );
    make_zero2one( $r2, $g2, $b2 );

    abs( $r1 - $r2 ) + abs( $g1 - $g2 ) + abs( $b1 - $b2 );
}


#===============================================================================
#  A few other functions.
#-------------------------------------------------------------------------------
#  Blend 2 or more RGB colors (actually takes any triplets):
#
#    \@color = blend_rgb_colors( $color1, $color2, ... )
#     @color = blend_rgb_colors( $color1, $color2, ... )
#
#-------------------------------------------------------------------------------
sub blend_rgb_colors
{
    my @rgb = ( 0, 0, 0 );
    foreach ( @_ )
    {
        next if ! ( $_ && ( ref($_) eq 'ARRAY' ) && ( @$_ >= 3 ) );  #  Skip bad colors
        my @clr = map { zero2one( $_ ) } @$_;
        foreach ( @rgb ) { $_ += shift @clr }
    }
    if ( @_ ) { foreach ( @rgb ) { $_ /= @_ } }

    wantarray() ? @rgb : \@rgb;
}


#-------------------------------------------------------------------------------
#  Blend 2 or more HTML colors (in linear RGB space) 
#
#     $html = blend_html_colors( $color1, $color2, ... )
#
#-------------------------------------------------------------------------------
sub blend_html_colors
{
    @_ or return '#000000';

    my @rgb = ( 0, 0, 0 );
    foreach ( @_ )
    {
        my @clr = html2rgb( $_ );
        foreach ( @rgb ) { $_ += shift @clr }
    }
    foreach ( @rgb ) { $_ /= @_ }

    rgb2html( @rgb );
}


1;
