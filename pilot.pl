#!/usr/bin/perl

use warnings;      # avoid d'oh! bugs
use strict;        # avoid d'oh! bugs
$|=1;
use Data::Dumper;  # easy printing, sometimes only used during coding
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlast;

# In one master hash 3 hashes for the different categorizations of
# sequences hashes are stored containing an array for
# each species, and for each library this array contains the seq objects

my %seq;
my @files = </home/ele/Data/pilot/trace2seq/trimmed_fasta/*>;

for my $file (@files) {
  my $in  = Bio::SeqIO->new(-file => "$file" , '-format' => 'Fasta');
  while ( my $seq = $in->next_seq ) {
    my $sequence= $seq->seq;
    my $id= $seq -> display_id; 
    if ($id =~ m/(?<lib>(?<spec>^\w\w)_(\w|\d){2,5})_.*/g) { 
      my $spec=$+{spec};                                 
      my $lib=$+{lib};
      $lib =~ s/(f|r|e|_)$//;
      push (@{$seq{"$spec"."_total"}{'total'}},$seq); 
      push (@{$seq{$lib}{'total'}},$seq);
      my $lengthfac = length($sequence)/70;     # get the lenthfraction
      $lengthfac = sprintf("%.0f",  $lengthfac);# and round
      if (length($sequence)<150) {
        push (@{$seq{"$spec"."_total"}{'short'}},$seq);
        push (@{$seq{$lib}{'short'}},$seq);
      }    
      elsif($sequence =~ m/(.*A{5,}|T{5,}|G{5,}|C{5,}.*){$lengthfac,}/ig) {
        push (@{$seq{"$spec"."_total"}{'poly'}},$seq);
        push (@{$seq{$lib}{'poly'}},$seq);
      }
      else {
        push (@{$seq{"$spec"."_total"}{'good1'}},$seq); 
        push (@{$seq{$lib}{'good1'}},$seq); 
      }
    }
  }
}
#print Dumper (\%seq);

# write the sequences still regardet good to files

my $out  = Bio::SeqIO->new(-file => ">Ac_good1.fasta" , '-format' => 'Fasta');
for my $seq(@{$seq{Ac_total}{good1}}) {
  $out -> write_seq($seq);
}

$out  = Bio::SeqIO->new(-file => ">Aj_good1.fasta" , '-format' => 'Fasta');
for my $seq(@{$seq{Aj_total}{good1}}) {
  $out -> write_seq($seq);
}


__END__

=head1 NAME

pilot.pl

=head1 SYNOPSIS

Not intended for users called by pilot.Rnw

=head1 DESCRIPTION

A custom perl script designed to pipe into a R-script.
The purpose of the script is explained in a pdf.

=head1 ARGUMENTS

 --help      print Options and Arguments 
 --man       print complete man page

=head1 AUTHOR

Emanuel Heitlinger, emanuelheitlinger@gmail.com

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 by Emanuel Heitlinger

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
