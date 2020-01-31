package PhyloTree;
use PhyloNode;
use strict;

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = { 
		name => "anonymous",
		rooted => 0,
		Nleaves => 0,
		Ninternals => 0,
		Nnodes => 0,
		rootnode => 0,
		internallabels => 0,
		edgelengths => 0,
		tstring => 0,
		nodelist => [],	
		leaflist => [],
		internalslist => [],
		@_
	};
	bless ($self, $class);
	return $self;
}


sub getMRCA {
	my $self = shift;
	my $node1 = shift;
	my $node2 = shift;
	$self->UnmarkAllNodes($self->{rootnode});
	my $nodeptr = $node1;
	while ($nodeptr) {
		$nodeptr->{marked} = 1;
		$nodeptr = $nodeptr->{anc};
	}
	my $nodeptr = $node2;
	while ($nodeptr && $nodeptr->{marked} == 0 ) {
		$nodeptr = $nodeptr->{anc};
	}
	return $nodeptr;
}

#sub clearclusters_recursor {
#	my $self = shift;
#	my $node = shift;
#	if ($node) {
#		$node->{cluster} = [];
#		clearclusters_recursor($node->{child});
#		clearclusters_recursor($node->{sibling});
#	}
#}

sub remakeclusters_recursor {
	my $self = shift;
	my $node = shift;
	if ($node) {
		$self->remakeclusters_recursor($node->{child});
		$self->remakeclusters_recursor($node->{sibling});
		#print "HERE, for ".$node."!\n";
		if ($node->{leaf} != 1) {
			#my @newcluster;
			@{$node->{cluster}} = ();
			my $nn = $node->{child};
			while ($nn) {
				foreach my $kidscluster ( @{$nn->{cluster}} ) {
					push(@{$node->{cluster}},$kidscluster);
		#			print "PUT '".$kidscluster."' AT ".$node."\n";
				}
				$nn = $nn->{sibling};
			}
			@{$node->{cluster}} = sort(@{$node->{cluster}});
#			$node->{cluster} = \@sortedcluster;
		}
	}
}


sub makeclusters_recursor {
	my $self = shift;
	my $node = shift;
	my $translationhash = shift;
	if ($node) {
		$self->makeclusters_recursor($node->{child},$translationhash);
		$self->makeclusters_recursor($node->{sibling},$translationhash);
		#print "HERE, for ".$node."!\n";
		if ($node->{leaf} == 1) {
			my $leafnum = $$translationhash{$node->{label}};
		#	print "GOT LN OF '".$leafnum."' FOR ".$node->{label}."\n";
			@{$node->{cluster}} = ();
			push(@{$node->{cluster}},$leafnum);
		} else {
			#my @newcluster;
			@{$node->{cluster}} = ();
			my $nn = $node->{child};
			while ($nn) {
				foreach my $kidscluster ( @{$nn->{cluster}} ) {
					push(@{$node->{cluster}},$kidscluster);
		#			print "PUT '".$kidscluster."' AT ".$node."\n";
				}
				$nn = $nn->{sibling};
			}
			@{$node->{cluster}} = sort(@{$node->{cluster}});
#			$node->{cluster} = \@sortedcluster;
		}
	}
}

sub makeclusters {
	my $self = shift;
	my $labelstolabelnums = shift;
#	clearclusters_recursor ( $self->{rootnode} );
	$self->makeclusters_recursor ( $self->{rootnode},$labelstolabelnums );
}

sub clone {
	my $model = shift;
	my $self = $model->new(%$model, @_);
	#ALSO NEEDS TO CLONE ALL OF THE NODES
	
	return $self;
}
	
# Parser FSA states
my ($stGETNAME, $stGETINTERNODE, $stNEXTMOVE,
		$stFINISHCHILDREN, $stQUIT, $stACCEPTED) = (1, 2, 3, 4, 5, 6); 
# Parser Token types
my ($tkSTRING, $tkNUMBER, $tkOTHER, $tkENDOFSTRING, $tkBAD, $tkSPACE, $tkLPAR, $tkRPAR,
	$tkCOMMA, $tkSEMICOLON, $tkCOLON, $tkTAB, $tkNEWLINE) = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13);
# Parser Error types
my ($errSYNTAX, $errENDOFSTRING, $errMISSINGLPAR, $errUNBALANCED, $errSTACKNOTEMPTY, $errSEMICOLON) = 
	(1,2,3,4,5,6);
	
my %tkTypeLookup = (
	1 , "tkString",
	2 , "tkNumber",
	3 , "tkOther",
	4 , "tkEndOfString",
	5 , "tkBAD",
	6 , "tkSPACE",
	7 , "tkLPAR",
	8 , "tkRPAR",
	9 , "tkCOMMA",
	10 , "tkSEMICOLON",
	11 , "tkCOLON",
	12 , "tkTAB",
	13 , "tkNEWLINE"
);

sub UnmarkAllNodes {
	my $self = shift;
	my $node = shift;
	if ($node) {
		$self->UnmarkAllNodes( $node->{child} );
		$self->UnmarkAllNodes( $node->{sibling} );
		$node->{marked} = 0;
	}
}


sub makenodelistRecursor {
	my $self = shift;
	my $node = shift;
	if ($node) {
		push( @{$self->{nodelist}}, $node );
		$self->makenodelistRecursor($node->{child});
		$self->makenodelistRecursor($node->{sibling});
	}
}

sub makeleaflistRecursor {
	my $self = shift;
	my $node = shift;
	if ($node) {
		if ( $node->{leaf} == 1) { push( @{$self->{leaflist}}, $node ); }
		$self->makeleaflistRecursor($node->{child});
		$self->makeleaflistRecursor($node->{sibling});
	}
}

sub makeinternalslistRecursor {
	my $self = shift;
	my $node = shift;
	if ($node) {
		$self->makeinternalslistRecursor($node->{child});
		$self->makeinternalslistRecursor($node->{sibling});
		if ( $node->{leaf} == 0) { push( @{$self->{internalslist}}, $node ); }
	}
}

sub MakeNodeList {
	my $self = shift;
	$self->{nodelist} = [];
	$self->makenodelistRecursor($self->{rootnode});
}

sub MakeInternalsList {
	my $self = shift;
	$self->{internalslist} = [];
	$self->makeinternalslistRecursor($self->{rootnode});
}

sub MakeLeafList {
	my $self = shift;
	$self->{leaflist} = [];
	$self->makeleaflistRecursor($self->{rootnode});
}

sub MarkDownFrom {
	my $self = shift;
	my $node = shift;
	if ($node) {
		$node->{marked} = 1;
		$self->MarkDownFrom($node->{anc});
	}
}

sub MarkUpFromRecursor {
	my $self = shift;
	my $node = shift;
	if ($node) {
		$node->{marked} = 1;
		$self->MarkUpFromRecursor($node->{child});
		$self->MarkUpFromRecursor($node->{sibling});
	}
}

sub MarkUpFrom {
	my $self = shift;
	my $node = shift;
	if ($node) {
		$node->{marked} = 1;
		$self->MarkUpFromRecursor($node->{child});
	}
}

sub GetPathLengthDists {
	my $self = shift;
	print scalar ( @{$self->{leaflist}} )."\n";
	for (my $i = 0; $i != scalar ( @{$self->{leaflist}} ); $i++)
	{
		$a = ( @{$self->{leaflist}} )[$i];
		#print "A=".$a."->".$a->{label}."\n";
		print $a->{label};
		foreach (my $k = 0; $k != $i; $k++) { print "\t"; }
		foreach (my $j = ($i+1); $j != scalar ( @{$self->{leaflist}} ); $j++)
		{
			$b = ( @{$self->{leaflist}} )[$j];
			#print "B=".$b."->".$b->{label}."\n";
			$self->UnmarkAllNodes($self->{rootnode});
			$self->MarkDownFrom($b);
			my $totallen = 0;
			my $thisnode = $a;
			while( ($thisnode->{marked} == 0) && $thisnode != $self->{rootnode} ) {
				$totallen += $thisnode->{edgelength};
				$thisnode = $thisnode->{anc};
			}
		#	print "GOT ".$totallen." TO CROSSING POINT \n";
			if ($thisnode->{marked} == 0 ) { die("FATAL ERROR - NO MARKS FOUND IN GetPathLengthDists\n"); }
		#	print "BEGINNING WALK UP\n";
			while ($thisnode->{leaf} == 0 ) {
				my $nextnode = $thisnode->{child};
		#		print "NN=".$nextnode."\t".$nextnode->{marked}."\n";
			#	die();
				while ($nextnode->{marked} == 0) { $nextnode = $nextnode->{sibling}; } #print "NNW=".$nextnode."\t".$nextnode->{mark}."\n"; }
		#		print "GOT MARKED CHILD\n";
				$totallen += $nextnode->{edgelength};
				$thisnode = $nextnode;
			}
			print "\t".$totallen;
		}
		print "\n";
	}

}

sub WriteNHtraversal {
	my $self = shift;
	my $node = shift;
	my $showinternals = shift;
	if ($node) {
#		print "LOOKING AT NODE ".$node.", WITH CHILD=".$node->{child}." AND SIB=".$node->{sibling}." AND ANC=".$node->{anc}."\n";
		if ($node->{leaf} ) {
			print $node->{label};
			if ($self->{edgelengths}) {
				print ":".$node->{edgelength};
			}
		} else { print "("; }
		$self->WriteNHtraversal($node->{child},$showinternals);
		if ($node->{sibling}) { print ","; }
		elsif ($node != $self->{rootnode} ) {	
			print ")";
			if ($node->{anc}->{label} && $showinternals) {
				if ($node->{anc}->{label} ne "" ) {
					print "'".$node->{anc}->{label}."'";
				}
			}
			if ($self->{edgelengths} && ($node->{anc} != $self->{rootnode}) )
			{
				print ":".$node->{anc}->{edgelength};
			}
		}
		$self->WriteNHtraversal($node->{sibling},$showinternals);
	}
}

#Write a Tree In New Hampshire Format
sub writeNH {
	my $self = shift;
	my $showinternals = shift;
#	print "HAVE ".$self->{Nleaves}." AND ".$self->{Ninternals}." AND ".$self->{Nnodes}."\n";
	$self->WriteNHtraversal($self->{rootnode},$showinternals);
	print ";";
}

sub WriteNHstringtraversal {
	my $self = shift;
	my $node = shift;
	my $showinternals = shift;
	if ($node) {
#		print "LOOKING AT NODE ".$node.", WITH CHILD=".$node->{child}." AND SIB=".$node->{sibling}." AND ANC=".$node->{anc}."\n";
		if ($node->{leaf} ) {
			$self->{tstring} .= $node->{label};
			if ($self->{edgelengths}) {
				$self->{tstring} .= ":".$node->{edgelength};
			}
		} else { $self->{tstring} .= "("; }
		$self->WriteNHstringtraversal($node->{child},$showinternals);
		if ($node->{sibling}) {$self->{tstring} .= ","; }
		elsif ($node != $self->{rootnode} ) {	
			$self->{tstring} .= ")";
			if ($node->{anc}->{label} && $showinternals) {
				if ($node->{anc}->{label} ne "" ) {
					$self->{tstring} .= "'".$node->{anc}->{label}."'";
				}
			}
			if ($self->{edgelengths} && ($node->{anc} != $self->{rootnode}) )
			{
				$self->{tstring} .= ":".$node->{anc}->{edgelength};
			}
		}
		$self->WriteNHstringtraversal($node->{sibling},$showinternals);
	}
}

sub writeNHstring {
	my $self = shift;
	my $showinternals = shift;
	$self->{tstring} = "";
	$self->WriteNHstringtraversal($self->{rootnode},$showinternals);
	$self->{tstring} .= ";";
}


#PARSES A NEW HAMPSHIRE FORMAT TREE STRING
#RETURNS 0,0 FOR SUCCESS, OR POSITION AND TYPE OF ERROR
sub parseNH {
	#print "$tkTypeLookup{12}\n";
	my @stk; #does the job of the stack
#Makes a kind of toy FSA parser, based on Rod's code
	my $self = shift;
	my $treestring = shift;
#	print $self->{name}." IS PARSING ".$treestring."\n";
	my $curNode = new PhyloNode;
#	print "REF ROOT IS ".$curNode."\n";
	$self->{rootnode} = $curNode;
	my $initiallength = length($treestring);
	my $state = $stGETNAME;
	my $Error = 0;
	#while ($nextType != 4) {
	my ($nextType,$nextToken,$treestring) = getNextToken($treestring);
		
		while (($state != $stQUIT) && ($state != $stACCEPTED))
		{
		#	print "STATE IS ".$state." TOKEN IS ".$nextType." (".$tkTypeLookup{$nextType}.") nextToken is '".$nextToken."'\n";
		#	print " STRING IS".$treestring."\n";
			if ($state == $stGETNAME) {
				SWITCH : {
					if ($nextType == $tkSPACE || $nextType == $tkTAB || $nextType == $tkNEWLINE) {
						($nextType,$nextToken,$treestring) = getNextToken($treestring);		
						last SWITCH;
					} 
					if ($nextType == $tkSTRING) {
						# to do: handle translation
						$self->{Nleaves}++;
						$curNode->{leaf} = 1;
						$curNode->{leafnum} = $self->{Nleaves};
						$curNode->{weight} =  1;
						$curNode->{label} =  $nextToken;
						$curNode->{degree} = 0;
						($nextType,$nextToken,$treestring) = getNextToken($treestring);
						$state = $stGETINTERNODE;
						last SWITCH;
					}
					if ($nextType == $tkNUMBER) {
						#to do: handle translation
						$self->{Nleaves}++;
						$curNode->{leaf}  = 1;
						$curNode->{leafnum} = $self->{Nleaves};
						$curNode->{weight} = 1;
						$curNode->{label} = $nextToken;
						$curNode->{degree} = 0;
						($nextType,$nextToken,$treestring) = getNextToken($treestring);;
						$state = $stGETINTERNODE;
						last SWITCH;
					}
					if ($nextType == $tkLPAR) {
						$state = $stNEXTMOVE;
						last SWITCH;
					}
					if ($nextType == $tkENDOFSTRING) {
						$Error = $errENDOFSTRING;
						$state = $stQUIT;
						last SWITCH;
					}
					$Error = $errSYNTAX;
					$state = $stQUIT;
				}
			} elsif ($state == $stGETINTERNODE) {
				SWITCH : {
					if  ($nextType == $tkSPACE || $nextType == $tkTAB || $nextType == $tkNEWLINE) {
						($nextType,$nextToken,$treestring) = getNextToken($treestring);		
						last SWITCH;
					} 
					if ($nextType == $tkCOLON || $nextType == $tkCOMMA || $nextType == $tkRPAR ) {
						$state = $stNEXTMOVE;
						last SWITCH;
					}
					if ($nextType == $tkENDOFSTRING) {
						$Error = $errENDOFSTRING;
						$state = $stQUIT;
						last SWITCH;
					}
					$Error = $errSYNTAX;
					$state = $stQUIT;
				}
			} elsif ($state == $stNEXTMOVE) {
			
			#	print "DOING NXT MOVE\n";
				SWITCH : {
				#	if ($nextType == $tkSPACE || $nextType == $tkTAB || $nextType == $tkNEWLINE) {
				#		($nextType,$nextToken,$treestring) = getNextToken($treestring);		
				#		last SWITCH;
				#	} 
					if ($nextType == $tkCOLON) {
						$nextToken = $tkSPACE;
					#	while ($nextToken == $tkSPACE) { 
							($nextType,$nextToken,$treestring) = getNextToken($treestring);
					#		print "NT=".$nextToken."\n";
					#	}
						#f = atof (p.GetTokenAsCstr());
						$curNode->{edgelength} = $nextToken;
						$self->{edgelengths} = 1;
						#print "GOT EDGELENGTH OF ".$nextToken."\n";
						($nextType,$nextToken,$treestring) = getNextToken($treestring);
						last SWITCH;
					}
					if  ($nextType == $tkSPACE || $nextType == $tkTAB || $nextType == $tkNEWLINE) {
						($nextType,$nextToken,$treestring) = getNextToken($treestring);		
						last SWITCH;
					} 
					# The next node encountered will be a sibling
					# of Curnode and a descendant of the node on
					# the top of the node stack.
					if ($nextType == $tkCOMMA) {
						my $q = new PhyloNode;
						$curNode->{sibling} = $q;
						if (scalar(@stk) == 0)
						{
							$Error = $errMISSINGLPAR;
							$state = $stQUIT;
						}
						else
						{
							my $stktop = $stk[scalar(@stk)-1];
							$q->{anc} = $stktop;
							$stktop->AddWeight ( $curNode->weight );
							$stktop->IncrementDegree;
							$curNode = $q;
							$state = $stGETNAME;
							($nextType,$nextToken,$treestring) = getNextToken($treestring);
						}
						last SWITCH;
					}
					
					#The next node will be a child of CurNode, hence
					# we create the node and push CurNode onto the
					# node stack.
					if ($nextType == $tkLPAR) {
				#		print "DOING LPAR \n";
						$self->{Ninternals}++;
						push(@stk,$curNode);
						my $q = new PhyloNode;
						$curNode->{child} =  $q;
						$q->{anc} = $curNode;
						$curNode->IncrementDegree();
						$curNode = $q;
						($nextType,$nextToken,$treestring) = getNextToken($treestring);
						$state = $stGETNAME;
						last SWITCH;
				#		print "DONE LPAR NM\n";
					}
					
					#We've finished ready the descendants of the node
					#at the top of the node stack so pop it off.
					if ($nextType == $tkRPAR) {
						if (scalar(@stk) == 0)
						{
							$Error = $errUNBALANCED;
							$state = $stQUIT;
						}
						else
						{
							my $q = $stk[scalar(@stk)-1];
							$q->AddWeight ($curNode->weight);
							$curNode = $q;
							pop(@stk);
							$state = $stFINISHCHILDREN;
							($nextType,$nextToken,$treestring) = getNextToken($treestring);
					
						}
						last SWITCH;
					}
					if ($nextType == $tkSEMICOLON) {
						if (scalar(@stk) == 0)
						{
							$state = $stACCEPTED;
						}
						else
						{
							$Error = $errSTACKNOTEMPTY;
							$state = $stQUIT;
						}
						last SWITCH;
					}
					if ($nextType == $tkENDOFSTRING) {
						$Error = $errENDOFSTRING;
						$state = $stQUIT;
						last SWITCH;
					}
				#	print "DO I END UP HERE\n";
					$Error = $errSYNTAX;
					$state = $stQUIT;
				}
			} elsif ( $state == $stFINISHCHILDREN ) {
				SWITCH : {
					if ($nextType == $tkSTRING || $nextType == $tkNUMBER) {
						# internal label
						$self->{Internallabels} = 1;
						$curNode->{label} = $nextToken;
						($nextType,$nextToken,$treestring) = getNextToken($treestring);
						last SWITCH;
					}
					if ($nextType == $tkCOLON ) {
						($nextType,$nextToken,$treestring) = getNextToken($treestring);
						$curNode->{edgelength} = $nextToken;
						$self->{edgelengths} = 1;
						($nextType,$nextToken,$treestring) = getNextToken($treestring);
						last SWITCH;
				    }
				    if ($nextType == $tkSPACE || $nextType == $tkTAB || $nextType == $tkNEWLINE ) {
						($nextType,$nextToken,$treestring) = getNextToken($treestring);
						last SWITCH;
					}
					# We've completed traversing the descendants of the
					# node at the top of the stack, so pop it off.
					if ($nextType == $tkRPAR ) {
						if (scalar(@stk) == 0)
						{
							$Error = $errUNBALANCED;
							$state = $stQUIT;
						}
						else
						{
							my $q = $stk[scalar(@stk)-1];;
							$q->AddWeight ($curNode->weight);
							$curNode = $q;
							pop(@stk);
							($nextType,$nextToken,$treestring) = getNextToken($treestring);
						}
						last SWITCH;
					}
					# The node at the top of the stack still has some
					# descendants.
					if ($nextType == $tkCOMMA ) {
						my $q = new PhyloNode;
						$curNode->{sibling} = $q;
						if (scalar(@stk) == 0)
						{
							$Error = $errMISSINGLPAR;
							$state = $stQUIT;
						}
						else
						{
							my $stktop = $stk[scalar(@stk)-1];;
							$q->{anc} = $stktop;
							$stktop->AddWeight ( $curNode->weight );
							$stktop->IncrementDegree;
							$curNode = $q;
							$state = $stGETNAME;
							($nextType,$nextToken,$treestring) = getNextToken($treestring);
						}
						last SWITCH;
					}
					if ($nextType == $tkSEMICOLON) {
						$state = $stNEXTMOVE;
						last SWITCH;
					}
					if (scalar(@stk) == 0) { $Error = $errSEMICOLON; }	
					else { $Error = $errSYNTAX; }
					$state = $stQUIT;
				}
			}
		#	print "WHILEAWAY state=".$state." token=".$nextType."\n";
		}
	#print "ENDWHILE\n";
	#}
	# Handle errors
	my $pos = $initiallength - length($treestring);
	$self->{Nnodes} = $self->{Ninternals} + $self->{Nleaves};
	if ($state == $stQUIT)
	{
		
		#Clean up to go here
		return ($Error,$pos);
	}
	else
	{
		my $reftoroot = $self->{rootnode};
	#	print "REF ROOT IS ".$reftoroot."\n";
		PhyloNode::weight($reftoroot,$self->{Nleaves});
		return (0);
	}

	die();
}

sub getNextToken {
	my  $text = shift;
	my $token = "";
	my $tokenType = 0;
	my $pos = 0;
	my $ch = "";
	#print length($text)." ".$text."\n";
	if (length($text) == 0)
	{
		$token = "";
		$tokenType = $tkENDOFSTRING;
	} else {
		$ch = substr($text,$pos,1);
	#	print "CHAR IS '".$ch."'\n";
		if ($ch eq "'")
		{
			# Single quoted token
			my $done = 0;
			$pos++;
			while ($done == 0)
			{
				$ch = substr($text,$pos,1);
				# Look ahead for double quote
				if ($ch eq "'")
				{
					$ch = substr($text,$pos,1);
					$done = ($ch eq "'");
				}
				if ( ($done == 0) && ($ch ne "\n") && ($ch ne "\r"))
				{
					$token .= $ch;
				}
				$pos++;
			}
		#	print STDERR "done token iterator\n";
			$tokenType = $tkSTRING; #classify the token
		}
		elsif ( $ch =~ m/[a-zA-Z]/ )
		{
		#	print "Char Alpha\n";
			while ( $ch =~ m/[a-zA-Z0-9]/ || ($ch eq '_') || ($ch eq '.') && ($pos < length($text)) )
			{
				if ($ch eq '_') {
					$token .= ' ';
				} else {
					$token .= $ch;
				}
				$pos++;
				$ch = substr($text,$pos,1);
			}
			$tokenType = $tkSTRING; #classify the token
		}
		else
		{
			if ( $ch =~ m/[0-9]/)
			{
				
				while ($ch =~  m/[0-9.eE\-]/ && $pos < length($text) )
				{
					$token .=$ch;
					$pos++;
					$ch = substr($text,$pos,1);
				}
				$tokenType = $tkNUMBER;
			}
			elsif ( $ch eq '-' )
			{
				$token .= $ch;
				$pos++;
				$ch = substr($text,$pos,1);
				#print "NEGNUM ".$ch."\n";
				if ( $ch =~ m/[0-9]/)
				{
				#	print "GOT DIGIT\n";
					while ($ch =~  m/[0-9.eE\-]/ && $pos < length($text) )
					{
						$token .=$ch;
						$pos++;
						$ch = substr($text,$pos,1);
					}
					$tokenType = $tkNUMBER;
				#	print "NEG NUM ".$token."\n";
				} else {
					$pos--;
					$tokenType = $tkBAD; #just a minus ain't allowed 
				}
			}
			else
			{
				$token .= $ch;
				$pos++;
				my %tokensmap = (
					"(",$tkLPAR,
					")",$tkRPAR,
					" ",$tkSPACE,
					",",$tkCOMMA,
					";",$tkSEMICOLON,
					":",$tkCOLON,
					"\t",$tkTAB,
					"\r",$tkNEWLINE,
					"\n",$tkNEWLINE
				);
				#print "testing ".$tokensmap{"("}."\n";
				if (exists($tokensmap{$ch}) ) {
					$tokenType = $tokensmap{$ch};
				} else {
					$tokenType = $tkBAD;
				}
			}
		}
	}
	#Cut the rest of $text
	$text = substr($text,$pos);
#	print "returning ".$tokenType." ".$token." '".$text."'\n";
	return ($tokenType,$token,$text);
}


#THIS BIT OF NASTINESS GENERATES METHODS FOR ACCESSING THE LISTED DATA MEMBERS,
#THINGS NOT ON THIS LIST ARE EFFECTIVELY PRIVATE-ISH
for my $field (qw(name rooted rootnode Ninternals)) {
	my $slot = __PACKAGE__ . "::$field";
	no strict "refs";
	*$field = sub {
		my $self = shift;
		$self->{$slot} = shift if @_;
		return $self->{$slot};
	};
}


1;
