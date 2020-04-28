#!/usr/bin/perl
use strict;
use warnings;
my $fh;

my @wanttissues = ("Lung","Colon","Kidney","Liver","Small Intestine");
my %tissuehash;
foreach my $t (@wanttissues){
    $tissuehash{$t} = 1;
}

my %clinhash;
open ($fh,"<sample-info.txt") || die "no";
<$fh>;
while (<$fh>){
    s/[\r\n]//g;
    s/"//g;
    my $line = $_;
    my @lines = split("\t",$line);
    my $tissue = $lines[1];

    if (defined $tissuehash{$tissue}){
        $clinhash{$lines[0]} = 1;
        my $sid = $lines[0];
        $sid =~ s/-/./g; # change hyphens to dots to be consistent with R column names later
        print $sid;
        for (my $i = 1; $i < scalar @lines; $i++){
            print "\t$lines[$i]";
        }
        print "\n";
    }
}
close $fh;