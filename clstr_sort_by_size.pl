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

#usage ./clstr_sort_by_size.pl < sth.clstr > sorted_sth.clstr

my $sort_by_what = "no";

my @clstr = ();
my $no = 0;
my $len = 0;
my $clstr = "";

while($ll = <>) {
  if ($ll =~ /^>/) {
    if ($clstr) {
      push(@clstr, [$len, $no, $clstr]);
    }
    $no = 0;
    $len = 0;
    $clstr = "";
  }
  else {
    $clstr .= $ll;
    $no++;
    chop($ll);
    my @lls = split(/\t/,$ll);
    my $this_len = $lls[2];
    if ($this_len > $len) {$len = $this_len;}
  }
}
    if ($clstr) {
      push(@clstr, [$len, $no, $clstr]);
    }

if ($sort_by_what eq "no") {
  @clstr = sort {$b->[1] <=> $a->[1]} @clstr;
}
elsif ($sort_by_what eq "len") {
  @clstr = sort {$b->[0] <=> $a->[0]} @clstr;
}

my $clstr_no = 0;
foreach $c (@clstr) {
  print ">Cluster $clstr_no $c->[1]\n";
  print $c->[-1];
  $clstr_no++;
}

