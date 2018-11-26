#!/usr/bin/env perl
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
my $x_col;
my $y_col;
my $z_col;
my @commands;
my $output_filename;
my $force_output_file;

GetOptions("x=i" => \$x_col,
           "y=i", => \$y_col,
           "z=i", => \$z_col,
           "cmd|c=s@", => \@commands,
           "output|o=s" => \$output_filename,
           "force|f" => \$force_output_file,
           "help|h" => \$print_help,
    )
    or die "GetOptions";

sub state_usage {
    say "Usage: plot_munger.pl OPTIONS input-ascii-file.dat";
}

sub state_help {
    state_usage ();
        say '
 Options:
  -h, --help                  give this help screen.
  -x                          x column index (starts with 1, not 0)
  -y                          y column index (starts with 1, not 0)
  -z                          z column index (starts with 1, not 0)
  -c, --cmd                   transformation command, see below; can be provided
                              multiple times
  -o, --output                basename for output data file
  -f, --force                 overwrite existing files

 Known commands:
  - dzdy      : partial differentiate
  - abs        : absolute value
  - log        : natural logarithm
  - log10      : base 10 logarithm
  - min=$min   : lower cutoff value 
  - max=$max   : upper cutoff value
  - add=$value : add value

 Example: Calculate ln(|dzdy|):
 plot_munger.pl -x 1 -y 2 -z 3 --cmd dzdy --cmd abs --cmd log data.dat
';
    
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
my %unique_test = map {$_ => 1} @col_indices;
if (keys %unique_test != 3) {
    die "x,y, and z need to be different";
}

my $filename = $ARGV[0];
if (not defined $filename) {
    warn "Error: no input file given\nTry 'plot_munger.pl --help'\n";
    exit(1);
}


my $blocks = read_gnuplot_format($filename);
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

# x11 x21
# x12 x22
# x13 x23

# y11 y21
# y12 y22
# y13 y23

# and same for z

my @maps = dog($blocks->xchg( 0, 2 ));
# first dim: block, second dim: line

my $x_block = $maps[$x_col-1];
my $y_block = $maps[$y_col-1];
my $z_block = $maps[$z_col-1];

splot($x_block, $y_block, $z_block);

($x_block, $y_block, $z_block) = apply_commands(
    [@commands], $x_block, $y_block, $z_block
    );

splot($x_block, $y_block, $z_block);

if (defined $output_filename) {
    $output_filename .= '_' . join('_', @commands) . '.dat';
    if (-e $output_filename and not $force_output_file) {
        die "file $output_filename already exists";
    }
    warn "writing output to $output_filename\n";
    write_output_datafile($output_filename, $x_block, $y_block, $z_block);
    
}

sub write_output_datafile {
    my ($filename, $x_block, $y_block, $z_block) = @_;
    open my $fh, '>', $filename
        or die "cannot open file $filename: $!";

    print {$fh} "# x\ty\tz\n";
    my $blocks = cat($x_block, $y_block, $z_block)->xchg(0,2);
    my @blocks = dog $blocks;
    for my $block (@blocks) {
        write_output_datafile_block($fh, $block);
    }
}

sub write_output_datafile_block {
    my ($fh, $block) = @_;
    my @lines = dog($block);
    LINE: for my $line (@lines) {
        my @cols = @{unpdl $line};
        my $output = "";
        while (my ($idx, $col) = each (@cols) ){
            if ($col eq 'NaN') {
                warn "NaN detected, skipping line '$line'";
                next LINE;
            }
            $output .= sprintf("%.10g", $col);
            if ($idx != $#cols) {
                    $output .= "\t";
            }
        }
        $output .= "\n";
        print {$fh} $output;
    }
    # finish block
    print {$fh} "\n";
}


sub apply_commands {
    my ($cmds, @blocks) = @_;
    my @commands = @{$cmds};
    for my $cmd (@commands) {
        @blocks = apply_command($cmd, @blocks);
    }
    return @blocks;
}

sub apply_command {
    my ($cmd, $x_block, $y_block, $z_block) = @_;
    
    if ($cmd eq 'abs' or $cmd eq 'log' or $cmd eq 'log10') {
        $z_block = $z_block->$cmd;
    }
    elsif ($cmd eq 'dzdy') {
        my $diff_y = diff_along_second_dim($y_block);
        my $diff_z = diff_along_second_dim($z_block);
        $x_block = $x_block->slice(':,1:');
        $y_block = $y_block->slice(':,1:');
        $z_block = $diff_z / $diff_y;
    }
    elsif ($cmd =~ /^(min|max|add)=(.*)/) {
        # value command
        $cmd = $1;
        my $value = $2;
        if (not looks_like_number($value)) {
            die "$value not a number";
        }
        if ($cmd eq 'max') {
            my $mask = $z_block > $value;
            $z_block = $z_block * (1 - $mask) + $mask * $value;
        }
        elsif ($cmd eq 'min') {
            my $mask = $z_block > $value;
            $z_block = $z_block * $mask + (1 - $mask) * $value;
        }
        elsif ($cmd eq 'add') {
            $z_block = $z_block + $value;
        }
    }
    else {
        die "unknown command $cmd";
    }
    return ($x_block, $y_block, $z_block);
}

sub diff_along_second_dim {
    my $pdl  = shift;
    return ($pdl->slice(':,1:') - $pdl->slice(':,:-2'));
}



sub splot {
    my ($x_block, $y_block, $z_block) = @_;
    my $terminal = 'qt';
    my %terminal_options = ();
    my %plot_options = (
        pm3d => 'implicit map corners2color c1',
        surface => 0,
        clut => 'sepia',
        cbtics => {format => "%g"},
        grid => 1);
    my $plot = PDL::Graphics::Gnuplot->new($terminal, %terminal_options, \%plot_options);
    $plot->splot($x_block, $y_block, $z_block);
}


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
