#!/usr/bin/perl

while($ll=<>){
	chomp($ll);
	if($ll=~/^>/){
		$ll=~/^>(\d+)/;$id=$1;
		$l=0;
		for($i=0;$i<6;$i++){
			$ll=<>;chomp($ll);
			if($ll=~/Possible/){next;}
			$fingerprint.=$ll;$l++;
		}
		if($l==5){
			$ll=<>;chomp($ll);
			$fingerprint.=$ll;
		}
		$fingerprint=~s/ //g;
		$len=length($fingerprint);if($len !=256){print "#error: $id\t$fingerprint\n";}
		print "$id\t$fingerprint\n";
		$fingerprint="";
	}
}
