#!/usr/bin/env perl
use 5.020;
use warnings;
use strict;

my $file = $ARGV[0] or die "need filename argument";
open my $fh, '>', $file or die "cannot open file";
say {$fh} "# x\t\ty\t\tval1\t\tval2";
my $N = 100;

for my $i (1..$N) {
    for my $j (1..$N) {
        my $value = sin($i/10 - $j/10);
        say {$fh} "$i $j $value $value";
    }
    print {$fh} "\n";
}
