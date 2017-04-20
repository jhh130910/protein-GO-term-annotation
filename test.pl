#!/usr/bin/perl -w

for (1..10){
   print &add($_);
}

print &add(1,2,3,4);

#print &help()

sub add {
        my @arr = @_;
        %h = map { $_,$_*2;}@arr;
        foreach ( keys %h ){
                print "Raw_data: $_\t Multiply: $h{$_}\n";
        }
}

if ($$$){
        $a->[1] > $a->[2] ? &help() : 0 ;

        return;
}

=head1 Description 
  test file 
=cut
