#!/usr/bin/perl
use strict;
use warnings;
my $fh;
my @files;
open ($fh,"<wald-difex-files.txt") || die "c";
while (<$fh>){
    s/[\r\n]//g;
    push(@files,$_);
}

my %hash;
open ($fh,"<degreport-p0.05-cov-genes-table.upC.downC.txt") || die "c";
<$fh>;
while (<$fh>){
    s/[\r\n]//g;
    s/"//g;
    my $line = $_;
    my @lines = split ("\t",$line);
    $hash{$lines[0]} = $lines[1];
}
close $fh;

foreach my $f (@files){
    my $ofh;
    my $outf = $f;
    $outf =~ s/.txt/.Aging.txt/g;
    open ($ofh,">$outf") || die "C";
    open ($fh,"<$f") || die "C";
    my $head = <$fh>;
    $head =~ s/[\r\n]//g;
    print $ofh "$head\tAgingCluster\n";

    while (<$fh>){
        s/[\r\n]//g;
        s/"//g;
        my $line = $_;
        my @lines = split ("\t",$line);
        if (defined $hash{$lines[1]}){
            for (my $i = 1; $i < scalar @lines; $i++){
                print $ofh "$lines[$i]\t";
            }
            print $ofh $hash{$lines[1]},"\n";
        }
        if (!defined $hash{$lines[1]}){
            for (my $i = 1; $i < scalar @lines; $i++){
                print $ofh "$lines[$i]\t";
            }
            print $ofh "0","\n";
        }
    }
    close $fh;
    close $ofh;
}
