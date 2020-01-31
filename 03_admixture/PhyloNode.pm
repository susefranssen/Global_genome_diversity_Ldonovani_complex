package PhyloNode;
use strict;
use Class::Struct;
struct PhyloNode => { #Definition of a PhyloNode
	label => '$',
	leaf => '$',
	sibling => 'PhyloNode',
	child => 'PhyloNode',
	anc => 'PhyloNode',	
	branchlen => '$', #this never seems to get used.. Not sure why! - use edgelength
	leafnum => '$',
	weight => '$',
	degree => '$',
	marked => 0,
	cluster => '@',
	edgelength => 0
};

sub AddWeight {
	my $referent = shift;
	$referent->{weight} += shift;
}

sub IncrementDegree {
	my $self = shift;
	$self->{degreee}++;
}

1;


	


