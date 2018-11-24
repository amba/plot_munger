#!/usr/bin/env perl
use 5.010;
use warnings;
use strict;

use PDL;
use MooseX::Params::Validate 'validated_list';
use Carp;
use Data::Dumper;
use Getopt::Long qw/:config gnu_getopt/;;
use PDL::Graphics::Gnuplot ();
my $print_help;
my $x_col;
my $y_col;
my $z_col;
GetOptions("x=i" => \$x_col,
           "y=i", => \$y_col,
           "z=i", => \$z_col,
           "help|h" => \$print_help,
    )
    or die "GetOptions";

say "y_col: $y_col";

sub state_usage {
    say "Usage: plot_munger.pl OPTIONS input-ascii-file.dat";
}

sub state_help {
    state_usage ();
        say "
 Options:
  -h, --help                  give this help screen.
";
    
}

if ($print_help) {
    state_help();
    exit(0);
}

if (not defined $x_col or not defined $y_col
    or not defined $z_col) {
    die "need x, y, and z options";
}
my @col_indices =($x_col, $y_col, $z_col);
say "@col_indices";
my %unique_test = map {$_ => 1} @col_indices;
say (keys %unique_test);
if (keys %unique_test != 3) {
    die "x,y, and z need to be different";
}

my $filename = $ARGV[0];
if (not defined $filename) {
    warn "Error: no input file given\nTry 'plot_munger.pl --help'\n";
    exit(1);
}


my $blocks = read_gnuplot_format(file => $filename);
# dims of $blocks PDL: cols, line_in_block, block
say "datafile shape: ", $blocks->shape;

# 3D gnuplot data file (two x values, three y value):
# x11 y11 z11
# x12 y12 z12
# x13 y13 z13
#
# x21 y21 z21
# x22 y22 z22
# x23 y23 z23
#
# x-cordinate changes with each block: x11, x12 and x13 will be equal
# in most cases (exception: sweep of B-field or temperature where
# they will be almost equal.
#
# Parse into three 2x3 piddles for x, y and z data
# (first piddle dim (x) goes to the right):

# x11 x21
# x12 x22
# x13 x23

# y11 y21
# y12 y22
# y13 y23

# and same for z

my @maps = dog($blocks->xchg( 0, 2 ));
say $_->shape for @maps;
# first dim: block, second dim: line
my $x_block = $maps[$x_col-1];
my $y_block = $maps[$y_col-1];
my $z_block = $maps[$z_col-1];

my $terminal => 'qt';
my %terminal_options = ();
my %plot_options = (pm3d => 'implicit map corners2color c1', surface => 0, clut => 'sepia', grid => 1);
my $plot1 = PDL::Graphics::Gnuplot->new('qt', %terminal_options, \%plot_options);
$plot1->splot($x_block, $y_block, $z_block);

my $diff_y = diff_along_second_dim($y_block);
my $diff_z = diff_along_second_dim($z_block);
my $reduced_x_block = $x_block->slice(':,1:');
my $reduced_y_block = $y_block->slice(':,1:');
my $plot2 = PDL::Graphics::Gnuplot->new('qt', %terminal_options, \%plot_options);
$plot2->splot($reduced_x_block, $reduced_y_block, $diff_z);
$plot2->splot($reduced_x_block, $reduced_y_block, log(abs($diff_z / $diff_y)));

sub diff_along_second_dim {
    my $pdl  = shift;
    return ($pdl->slice(':,1:') - $pdl->slice(':,:-2'));
}








# produce 2D PDL for each block. Cat them into a 3d PDL
sub get_blocks {
    my ( $fh ) = validated_list(
        \@_,
        fh          => { isa => 'FileHandle'},
        );

    my @blocks;
    my @rows;
    my $num_columns;
    while ( my $line = <$fh> ) {
        if ( $line =~ /^#/ ) {
            next;
        }
        if ( $line =~ /^\s*$/ ) {

            # Finish block. Need check for number of rows if we have
            # multiple subsequent blank lines
            if ( @rows > 0 ) {

                # Give \@rows, not @rows to get a 2D piddle if we
                # only have a single row.
                push @blocks, pdl( \@rows );
                @rows = ();
            }
            next;
        }

        # awk splitting behaviour
        my @nums = split( ' ', $line );
        if (not defined $num_columns) {
            $num_columns = @nums;
        }
        else {
            if (@nums != $num_columns) {
                die "num cols not $num_columns";
            }
        }
        push @rows, [@nums];
    }
    if ( @rows > 0 ) {
        push @blocks, pdl( \@rows );
    }

    # bring blocks to same number of rows: reshape and add NaNs.
    my $max_rows = List::Util::max( map { ( $_->dims )[1] } @blocks );

    for my $block (@blocks) {
        my $rows = ( $block->dims() )[1];
        if ( $rows < $max_rows ) {
            $block->reshape( $num_columns, $max_rows );
            $block->slice(":,${rows}:-1") .= "NaN";
        }
    }

    return PDL::cat(@blocks);
}

sub read_gnuplot_format {
    my ( $file ) = validated_list(
        \@_,
        file        => { isa => 'Str'},
        );
    open my $fh, '<', $file
            or croak "cannot open file $file: $!";

    # return value is 3D PDL with following dims
    # 0st dim: column
    # 1st dim: row (in block)
    # 2nd dim: block
    return get_blocks( fh => $fh);
}
