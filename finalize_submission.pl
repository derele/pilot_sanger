#!/usr/bin/perl

use warnings;      # avoid d'oh! bugs
use strict;        # avoid d'oh! bugs
$|=1;
use Data::Dumper;  # easy printing, sometimes only used during coding


my $similar = "";
my $primer = "M13F(GTAAAACGACGGCCAGT)";

while (<>) {
  if ($_ =~ m/(PUT_ID:)(.*) (Score = .*)/g){
    print("$1\n");
    $similar = $2;
  }
  elsif ($_ =~ m/(PUT_ID:)\n/g){
    print("$1\n");
    $similar = "";
  }
  elsif ($_ =~ m/(COMMENT:)/g) {
    print("$1$similar\n");
  }
  elsif ($_ =~ m/CLONE: \w+_*+r_\d+\w\d+/) {
    $primer = "M13R(GGCAGGAAACAGCTATGACC)";
    print $_ ;
  }
  elsif ($_ =~ m/CLONE: \w+_*+f_\d+\w\d+/) {
    $primer = "M13F(GTAAAACGACGGCCAGT)";
    print $_ ;
  }
  elsif ($_ =~ m/CLONE: \w+_[A-Z|1-9]+_\d+\w\d+/) {
    $primer = "M13F(GTAAAACGACGGCCAGT)";
    print $_ ;
  }
  elsif ($_ =~ m/(SEQ_PRIMER:) M13/) {
    print ("$1 $primer\n") ;
  }
  
  else {
    print $_ ;
  } 
}


