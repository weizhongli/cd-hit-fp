#!/usr/bin/perl
// =============================================================================
// CD-HIT-FP
//
// CD-HIT-FP is a very fast clustering program to cluster similar compounds 
// from a very large small molecule library.
//
// program written by
//                                      Weizhong Li
//                                      Email liwz@sdsc.edu
// Reference:
// Weizhong Li. A Fast Clustering Algorithm for Analyzing Highly Similar 
// Compounds of Very Large Libraries. J. Chem. Inf. Model. (2006) 46:1919-1923
// =============================================================================

print STDERR<<EOD;
This script translate a fingerprint file from Unity format into generic format
usage: fp_unity2generic.pl < input_unity_fprint > outout_generic_fprint

EOD

my $ll;
while($ll=<>){
    chop($ll);
    next if ($ll =~ /^#/);
 
    if ($ll =~ /^\d/) {
      my ($names, $sln) = split(/\t/, $ll);
      my ($d1, $d2, $id) = split(/\s+/, $names, 3);
      $ll = <>; $ll=~ s/\W//g;
      my $fprint = $ll;
      print "$id\t$fprint\n";
      $ll = <>;
    }
    else {
      next;
    }
}
