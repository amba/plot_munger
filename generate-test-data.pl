#!/usr/bin/env perl
use 5.020;
use warnings;
use strict;

my $file = $ARGV[0] or die "need filename argument";
open my $fh, '>', $file or die "cannot open file";
my $N = 100;

for my $i (1..$N) {
    for my $j (1..$N) {
        my $value;
        if ($i % 20 == 0 || $j % 20 == 0) {
            $value = 1;
        }
        else {
            $value = 0;
        }
        say {$fh} "$i $j $value $value";
    }
    print {$fh} "\n";
}
