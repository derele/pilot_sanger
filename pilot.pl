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
my @files = <trace2seq/trimmed_fasta/*>;

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

## blastsearches... this is were reproducability gets lost as long as you 
##don't have the same databases I am using. Databases can passibly be 
##included in the whole thing, for the meantime I just comment out the 
##calls that produce the comparison-tables of the blast-searches and use this files.

system("./cobl.pl -i Ac_good1.fasta -d /home/ele/db/blastdb/nempep4/nempep4 -p blastx -d /home/ele/db/blastdb/fishpep1/fishpep1 -p blastx -d /home/ele/db/blastdb/fishrRNA/fishrRNA -p blastn -d /home/ele/db/blastdb/nemrRNA/nemrRNA -p blastn -o Ac_cobl_table.csv -e 0.6");

system("./cobl.pl -i Aj_good1.fasta -d /home/ele/db/blastdb/fishrRNA/fishrRNA -p blastn -d /home/ele/db/blastdb/fishpep1/fishpep1 -p blastx -o Aj_cobl_table.csv -e 0.6");

#read the cobl.pl (blast-comparison) tables
open(ACCOBL, "Ac_cobl_table.csv") or die "can`t open file: $!";
open(AJCOBL, "Aj_cobl_table.csv") or die "can`t open file: $!";;
my @cobl = (<ACCOBL>, <AJCOBL>);
close (ACCOBL);
close (AJCOBL);

my %gc_of;

COBL:for (@cobl) {
  my @temp = split (/\t/, $_ );
  my $read = $temp[0];
  my $db = $temp[1];
  my $annot = $temp[4];
 SEQOB:for my $spec("Ac_total", "Aj_total"){
    for (@{$seq{$spec}{'good1'}}){
      my $seq_obj = $_;
      my $id =  $seq_obj -> primary_id;
      my $gc = gc_percent($seq_obj -> seq);
      if ($id =~ m/($read)/) {
        my $lib=$+;
        $lib =~ s/(f|r)*_(\d){2,3}_*\w\d\d$//;
        if ($db =~ m/fishpep1/ and $spec =~ m/Ac/) {
          $seq_obj -> desc("similar to "."$annot");
          push (@{$seq{$spec}{'fishpep'}}, $seq_obj);
          push (@{$seq{$lib}{'fishpep'}}, $seq_obj);
          $gc_of{$spec}{'fishpep'}{$id} = $gc;
          next COBL;
        }
        if ($db =~ m/rRNA/) {
          $seq_obj -> desc("similar to "."$annot");
          push (@{$seq{$spec}{'rRNA'}}, $seq_obj);
          push (@{$seq{$lib}{'rRNA'}}, $seq_obj);
          next COBL;
        }
        push (@{$seq{$spec}{'good2'}}, $seq_obj);
        push (@{$seq{$lib}{'good2'}}, $seq_obj);
        $gc_of{$spec}{'good2'}{$id} = $gc;
        next COBL;
      }
    }
  }
}

#print out the numbers per library 
while (my($spec, $outer) = each(%seq)) {
  while (my ($diagnosis, $middle) = each(%$outer)){
    my $number = scalar @$middle;
    print "$spec,$diagnosis,$number,count\n";
  }
}

#print out a table of GC contents
while (my ($spec, $read_hash) = each (%gc_of)) {
  while (my ($diagnosis, $inner_hash)= each (%$read_hash)) {
    while (my ($read, $gc)= each (%$inner_hash)) {
      print "$spec,$diagnosis,$read,$gc\n",
    }
  }
}

sub gc_percent{
  my $sequence=$_[0];
  my $gc_count;
  my @DNA=split('',$sequence);
  for my $base(@DNA){
    if ($base =~ /^(G|C)/i){
      ++$gc_count;
    }
  }
  my $gc_percent = ($gc_count/length($sequence))*100;
  $gc_percent = sprintf("%.2f", $gc_percent);
  return($gc_percent);
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
