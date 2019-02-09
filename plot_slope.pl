#!/usr/bin/env perl

#
# Work in progress
#

use 5.010;
use warnings;
use strict;

use PDL;
use Carp;
use Data::Dumper;
use Getopt::Long qw/:config gnu_getopt/;;
use PDL::Graphics::Gnuplot ();
use Scalar::Util 'looks_like_number';

my $print_help;
my @commands;
my $output_filename;
my $force_output_file;
my $palette = "bbwr";
my $xrange;
my $yrange;
my $logscale;
my $cbrange;

GetOptions(
      #     "cmd|c=s@", => \@commands,
    #"output|o=s" => \$output_filename,
    # "xrange|x=s" => \$xrange,
    # "yrange|y=s" => \$yrange,
    #"force|f" => \$force_output_file,
    # "logscale" => \$logscale,
    # "palette=s" => \$palette,
    # "cbrange=s" => \$cbrange,
    "help|h" => \$print_help,
    
    )
    or die "GetOptions";

sub state_usage {
    say "Usage: plot_slope.pl OPTIONS input-ascii-file.dat xcol:zcol";
}

sub state_help {
    state_usage ();
        say '
 Options:
  -h, --help                  give this help screen.

';
    
}

if ($print_help) {
    state_help();
    exit(0);
}

my $filename = $ARGV[0];
if (not defined $filename) {
    warn "Error: no input file given\nTry 'plot_munger.pl --help'\n";
    exit(1);
}

my $xy_arg = $ARGV[1];
if (not defined $xy_arg) {
    die "need xcol:ycol:zcol argument";
}

if ($xy_arg !~ /^([0-9]+):([0-9]+)/) {
    die "$xy_arg needs format xcol:ycol";
}
my ($x_col, $y_col) = ($1, $2);


my @col_indices =($x_col, $y_col);
my %unique_test = map {$_ => 1} @col_indices;
if (keys %unique_test != 2) {
    die "x and y need to be different";
}


my $blocks = read_gnuplot_format($filename);
my $blocks_xy = $blocks->dice([$x_col-1, $y_col-1],'X','X');

# dims of $blocks PDL: cols, line_in_block, block

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


my @dims = $blocks_xy->dims;
say "dims: @dims";

my $first_x = $blocks_xy->at(0,0,0);
my $last_x = $blocks_xy->at(0,0,-1);

my $first_y = $blocks_xy->at(1,0,0);
my $last_y = $blocks_xy->at(1,0,-1);

say "first x value: $first_x, last x value: $last_x";
say "first y value: $first_y, last y value: $last_y";

my $slope = abs(($first_x - $last_x) / ($first_y - $last_y));
say "slope: $slope";
say "slope * 3600: ", $slope * 3600;
# # assume regular data: x11 = x12 = x13, ...; y11 = y21 = y31, ...
# my $plot;
# my $output_block;
# if ($trace_index eq 'x') {
#     my $xvals = $blocks_xyz->slice('(0),(0),:')->unpdl;
#     my $best_index = find_closest_index($trace_value, $xvals);

#     my $yvals = $blocks_xyz->slice("(1),:,($best_index)");
#     my $zvals = $blocks_xyz->slice("(2),:,($best_index)");
    
#     $plot = plot($yvals, $zvals);
#     # output 2D PDL with dims (col, trace_point)
#     $output_block = $blocks->slice(":,:,($best_index)");
# }
# elsif ($trace_index eq 'y') {
#     my $yvals = $blocks_xyz->slice('(1),:,(0)')->unpdl;
#     my $best_index = find_closest_index($trace_value, $yvals);
#     my $xvals = $blocks_xyz->slice("(0),($best_index),:");
#     my $zvals = $blocks_xyz->slice("(2),($best_index),:");
#     $plot = plot($xvals, $zvals);
#     $output_block = $blocks->slice(":,($best_index),:");
# }


# if (defined $output_filename) {
#     $output_filename .= "_trace_$trace_index=$trace_value.dat";
#     if (-e $output_filename and not $force_output_file) {
#         die "file $output_filename already exists";
#     }
#     warn "writing output to $output_filename\n";
#     open my $fh, '>', $output_filename
#         or die "cannot open file $output_filename: $!";

#     write_output_datafile_block($fh, $output_block);
    
# }


# # keep process running. required for interactive features of plots
# sleep(100000);

# sub find_closest_index {
#     my $target = shift;
#     my $array_ref = shift;
#     my @array = @{$array_ref};
#     my $best_index;
#     my $min_deviation;
#     for my $index (0..$#array) {
#         my $deviation = abs($target - $array[$index]);
#         if (not defined $min_deviation or ($deviation < $min_deviation)) {
#             $min_deviation = $deviation;
#             $best_index = $index;
#         }
#     }
#     say "using x value ", $array[$best_index];
#     return $best_index;
# }






# # sub write_output_datafile {
# #     my ($filename, $output_block) = @_;
# #     open my $fh, '>', $filename
# #         or die "cannot open file $filename: $!";

# #     print {$fh} "# x\ty\tz\n";
# #     my $blocks = cat($x_block, $y_block, $z_block)->xchg(0,2);
# #     my @blocks = dog $blocks;
# #     for my $block (@blocks) {
# #         write_output_datafile_block($fh, $block);
# #     }
# # }

# sub write_output_datafile_block {
#     my ($fh, $block) = @_;
#     my @lines = dog($block);
#     LINE: for my $line (@lines) {
#         my @cols = @{unpdl $line};
#         my $output = "";
#         while (my ($idx, $col) = each (@cols) ){
#             if ($col eq 'NaN') {
#                 warn "NaN detected, skipping line '$line'";
#                 next LINE;
#             }
#             $output .= sprintf("%.10g", $col);
#             if ($idx != $#cols) {
#                     $output .= "\t";
#             }
#         }
#         $output .= "\n";
#         print {$fh} $output;
#     }
#     # finish block
#     print {$fh} "\n";
# }

# sub plot {
#     my ($y_block, $z_block) = @_;
#     my $terminal = 'qt';
#     my %terminal_options = ();
#     my %plot_options = (
#         title => "$trace_index=$trace_value",
#         grid => 1);

#     # if ($xrange) {
#     #     if ($xrange !~ /(.+):(.+)/) {
#     #         die "xrange arge needs format 'x1:x2'";
#     #     }
#     #     my ($x1, $x2) = ($1, $2);
#     #     if (not looks_like_number($x1) or not looks_like_number($x2)) {
#     #         die "xrange limits must be numbers";
#     #     }
#     #     %plot_options = (%plot_options, xrange => [$x1, $x2]);
#     # }
#     # if ($yrange) {
#     #     if ($yrange !~ /(.+):(.+)/) {
#     #         die "yrange arge needs format 'x1:x2'";
#     #     }
#     #     my ($y1, $y2) = ($1, $2);
#     #     if (not looks_like_number($y1) or not looks_like_number($y2)) {
#     #         die "yrange limits must be numbers";
#     #     }
#     #     %plot_options = (%plot_options, yrange => [$y1, $y2]);
#     # }
#     my $plot = PDL::Graphics::Gnuplot->new($terminal, %terminal_options, \%plot_options);

#     $plot->plot($y_block, $z_block);
#     return $plot;
# }


# produce 2D PDL for each block. Cat them into a 3d PDL
sub get_blocks {
    my $fh = shift;

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
    my $file = shift;
    open my $fh, '<', $file
            or die "cannot open file $file: $!";

    # return value is 3D PDL with following dims
    # 0st dim: column
    # 1st dim: row (in block)
    # 2nd dim: block
    return get_blocks($fh);
}
