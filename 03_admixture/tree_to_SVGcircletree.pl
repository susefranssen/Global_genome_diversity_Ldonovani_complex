#!/usr/bin/perl
#!/usr/bin/perl
use strict;
use warnings;
use PhyloTree;
use PhyloNode;
use Math::Trig;
 use List::Util qw(sum min max);
use POSIX qw/floor ceil/;
use Getopt::Long;

my $help = 0;

$| = 1;

my $HGItree = "(((amphimedon_queenslandica:0.5443040880954275,trichoplax_adhaerens:0.3953122846133697):1.00000050002909E-6,nematostella_vectensis:0.22398094747577915):0.0922564,((((((((trichobilharzia_regenti:0.10456343193508609,((((schistosoma_mattheei:0.015141019159580974,(schistosoma_haematobium:0.007090888917114553,schistosoma_curassoni:0.0061420916376505996):0.0023707660317076247):0.0036453357950648165,schistosoma_margrebowiei:0.013284666801122303):0.018672763515083474,(schistosoma_mansoni:0.011146368026663758,schistosoma_rodhaini:0.008896743528063718):0.019358934316114858):0.037636292950444704,schistosoma_japonicum:0.07676690262539353):0.03510259908118353):0.18434724378031422,((fasciola_hepatica:0.09043718513775853,echinostoma_caproni:0.07953746395799559):0.14719024368808126,clonorchis_sinensis:0.242757007642578):0.08041829285964996):0.14940169988353152,protopolystoma_xenopodis:0.47945618895051173):0.04325462946947227,((schistocephalus_solidus:0.0349559539268644,(diphyllobothrium_latum:0.059454177485766826,spirometra_erinaceieuropaei:0.056403927322079564):0.009570082746471682):0.18259126705059933,(mesocestoides_corti:0.11551188101780575,((hymenolepis_diminuta:0.04769489380172161,(hymenolepis_nana:0.03951330087434702,hymenolepis_microstoma:0.022550267684797183):0.033565247983616806):0.20296637462611863,((echinococcus_multilocularis:0.009740569974250652,echinococcus_granulosus:0.00504556581955149):0.03127917250748528,((taenia_asiatica:0.008631046010548602,taenia_solium:0.010331528289254443):0.028752479529428075,hydatigera_taeniaeformis:0.06297359454692343):0.011800959843500679):0.07077626458761359):0.06053824634942566):0.14922731223559302):0.22617109039059097):0.21773710085862757,schmidtea_mediterranea:0.5665926768662565):0.3358941392698669,(crassostrea_gigas:0.3910317213284673,capitella_teleta:0.31554257451327544):0.05679453390806605):0.02865926643438824,((((((panagrellus_redivivus:0.369231810357707,bursaphelenchus_xylophilus:0.46465362558437756):0.04999401512263331,(meloidogyne_hapla:0.2533919419107295,globodera_pallida:0.281132242466975):0.3173073719938062):0.07643940528222536,(rhabditophanes_kr3021:0.28113998125943296,(parastrongyloides_trichosuri:0.09735144511657838,((strongyloides_papillosus:0.016567024127332252,strongyloides_venezuelensis:0.0220529895830295):0.04937561722199901,(strongyloides_ratti:0.041673811685990464,strongyloides_stercoralis:0.0507905823151483):0.03401292055394807):0.055414542025801636):0.15780931368897227):0.3731145597872935):0.057885560446111824,(((enterobius_vermicularis:0.14032259417537848,syphacia_muris:0.18659398523647652):0.23757016283077853,((anisakis_simplex:0.17158521023505333,(((ascaris_suum:0.00447791910440041,ascaris_lumbricoides:0.004381606850447894):0.006639108813746893,parascaris_equorum:0.12773560060704994):0.07036478672854261,toxocara_canis:0.07540556953926114):0.02684328147648546):0.0982689139263542,((gongylonema_pulchrum:0.15126259456599558,((((elaeophora_elaphi:0.05950042294219867,(acanthocheilonema_viteae:0.04872463614740447,litomosoides_sigmodontis:0.08165942350688099):0.009573032482803567):0.014923058270534244,((wuchereria_bancrofti:0.02318475660462276,(brugia_pahangi:0.003689398765705496,(brugia_timori:0.01058684179133393,brugia_malayi:0.002774164989874138):0.00570281170241602):0.01758130020348728):0.043569136993611424,loa_loa:0.05552946606564138):0.009783860024658142):0.013887533669589834,(dirofilaria_immitis:0.06790222500473068,(onchocerca_flexuosa:0.038029745861596014,(onchocerca_ochengi:0.0031682368368410587,onchocerca_volvulus:0.004801574105576899):0.021209390203435607):0.03586748937046973):0.015307763131987152):0.07130057777650407,thelazia_callipaeda:0.17325099158568333):0.06485227395980092):0.17193239468849106,dracunculus_medinensis:0.36172239659452404):0.02912446693145849):0.04561951883872225):0.12530630068321127,(pristionchus_pacificus:0.4444537528785171,(caenorhabditis_elegans:0.3320351469419328,(((strongylus_vulgaris:0.11079053872330093,(cylicostephanus_goldi:0.0720317522724643,oesophagostomum_dentatum:0.041482155152890964):0.0172602833625667):0.03129534994020681,(((ancylostoma_caninum:0.013015152215173848,ancylostoma_duodenale:0.015610848446522788):0.0033129197346754746,ancylostoma_ceylanicum:0.04192876597572771):0.02901346280584023,necator_americanus:0.07444679171558675):0.007888864522624238):0.030757912405848103,(((nippostrongylus_brasiliensis:0.08152385549745127,heligmosomoides_bakeri:0.10790859998455278):0.02177081979177046,((haemonchus_contortus:0.007732427526416987,haemonchus_placei:0.026368859273294444):0.0380386183423617,teladorsagia_circumcincta:0.06797672134835334):0.03681714769102292):0.025669082246367946,(dictyocaulus_viviparus:0.15732215010794692,(angiostrongylus_cantonensis:0.0323232193304421,angiostrongylus_costaricensis:0.050166555716476054):0.06184845345632645):0.05286796233561285):0.01392303849566058):0.15625214450174735):0.09371748464442949):0.12505079620888182):0.027692649820317213):0.3816184126668093,(romanomermis_culicivorax:0.5575929074322425,(((trichinella_nativa:0.024927631849874926,trichinella_spiralis:0.023001304243239613):0.40549546727003655,(trichuris_muris:0.09649136648633516,(trichuris_trichiura:0.05153621223445387,trichuris_suis:0.012761234932885052):0.059548208143779494):0.2305585140570685):0.19994987136100173,soboliphyme_baturini:0.4779435592108415):0.06422441889530313):0.06298543976962022):0.1799350347804883,(ixodes_scapularis:0.4570979503001921,drosophila_melanogaster:0.5676504531061201):0.08948591629143347):0.022938957253102455):0.059609643116709035,(ciona_intestinalis:0.6145863204798903,(danio_rerio:0.1586765444769755,homo_sapiens:0.16812890467338645):0.2138315042157433):0.058289245736278326):0.0922564);";

my $pete_tree = "((schmidtea_mediterranea:1,macrostomum_lignano:1):1,((gyrodactylus_salaris:1,protopolystoma_xenopodis:1):1,((((((hydatigera_taeniaeformis:1,(taenia_saginata:1,taenia_asiatica:1,taenia_solium:1):1):1,(echinococcus_multilocularis:1,echinococcus_granulosus:1,echinococcus_canadensis:1):1):1,(hymenolepis_microstoma,hymenolepis_diminuta:1,hymenolepis_nana:1):1):1,mesocestoides_corti:1):1,((diphyllobothrium_latum:1,spirometra_erinaceieuropaei:1):1,schistocephalus_solidus:1):1):1,(((opisthorchis_viverrini:1,clonorchis_sinensis:1):1,(fasciola_hepatica:1,echinostoma_caproni:1):1):1,(trichobilharzia_regenti,(schistosoma_japonicum,((schistosoma_rodhaini:1,schistosoma_mansoni:1):1,(schistosoma_haematobium:1,schistosoma_margrebowiei:1,schistosoma_curassoni:1,schistosoma_mattheei:1):1):1):1):1):1):1):1):1;";

my $hpix = 800;
my $wpix = 800;

#start and end in radians, measured from x-axis.
my $start_angle = 0;
my $end_angle = 2 * pi;

#my $start_angle = 7.0 * pi / 6.0;
#my $end_angle = 11.0 * pi / 6.0;

my $multiply_by_pi = 0;

#there should be 4 'rings' - space_for_bars, space_for_wedges, space_for_labels, $diagram_base_diameter
#any spacing ("padding") should be built into these regions..
#space for text and first colors wedges
my $space_for_labels = 100;
#space for second and subsequent wedges
my $space_for_wedges = 20;
#space for bars araound wedges
my $space_for_bars = 150;
#gap between wedges
my $gap_between_wedges = 2;
#font
my $fontsize = 7;
#padding between edges of tree / wedge and text
my $text_padding = 2;
#gap between bars and edge of outermost wedge edge
my $bar_padding = 2;
my $highlight_bars = "";
my %highlight_bars;

my $bar_width = 10;
my $root_branch_length = 0;
my $no_neg_branches = 0;
#not yet implemented
my $small_gap_between_text = 2;
my $colors_to_clades = "";
my $numwedges = -1;
#does colors_to_clades have a title line to discard?
my $title_line_in_CCfile = 1;
my $colorsppbywedge = 0;
my $colorbarbywedge = 0;

my $textcolor = "#000000";
my $skipSVGheader =0;
my $hideSPPnames = 0;
my $drawlegend = 0;
my $legend_margin = 50;
my $legend_width = 0;
my $legend_height = 0;
#point to a file that has species -> colors for text labels
my $species_name_colors = "";
my $species_name_colors_column = 1;
my $also_color_dots = 0;

#----
my $hgi_fix = 0;

#calculate_if not set explicitly
#amount added/removed from edges of wedges to make pretty
my $wedge_angle_fiddle = 0;

my $helminthfambars = "";


#admixture stuff
#allow multiple qmatrix expressions
my @qmatrix;
my $pedfile = "";
my $qcolors = "";
my $fix_spp_names_Ldon=0;
my $annularsectors = 0;
my $annularsectorgap_p = 0.2;
my $min_drawable = 0.000011;
my $nomatchcolors = 0;
my $color_ref_ring = -1;
my $not_concentric = 0;
my $multicolors = 0;
my $learncolors =0;
my $labelwedges = 0;
my $wedgelabelnudge = 0.1;
my $wedgelabelsize = $fontsize;
#
my $ok = GetOptions("root_branch_length=f"=>\$root_branch_length,"no_negative_branches"=>\$no_neg_branches,"highlight_bars=s"=>\$highlight_bars,"wedgelabelsize=i"=>\$wedgelabelsize,"learncolors"=>\$learncolors,"label_wedges"=>\$labelwedges,"wedgelabelnudge=f"=>\$wedgelabelnudge,"also_color_dots"=>\$also_color_dots,"species_name_colors_column=i"=>\$species_name_colors_column,"species_name_colors=s"=>\$species_name_colors,"multicolors"=>\$multicolors,
"not_concentric"=>\$not_concentric,"color_ref_ring=i"=>\$color_ref_ring,"nomatchcolors"=>\$nomatchcolors,"min_drawable=f"=>\$min_drawable,"annularsectors"=>\$annularsectors,"annularsectorgap_p=f"=>\$annularsectorgap_p,"ldonglobal_fix_spp_names"=>\$fix_spp_names_Ldon,"qmatrix:s"=>\@qmatrix,"pedfile=s"=>\$pedfile,"qcolors=s"=>\$qcolors,"start_angle=f"=>\$start_angle,"end_angle=f"=>\$end_angle,"multiply_by_pi"=>\$multiply_by_pi,"legendwidth=i"=>\$legend_width,
"legendheight=i"=>\$legend_height,"drawlegend"=>\$drawlegend,"bar_width=i"=>\$bar_width,"space_for_labels=i"=>\$space_for_labels, "space_for_bars=i"=>\$space_for_bars, 
"space_for_wedges=i"=>\$space_for_wedges, "hideSPPnames"=>\$hideSPPnames,"skipSVGheader"=>\$skipSVGheader,"hpix=i"=>\$hpix,"wpix=i"=>\$wpix,
"50hgi_fix"=>\$hgi_fix,"50helminthsbars=s"=>\$helminthfambars, "help"=>\$help,"specieswedge=s"=>\$colors_to_clades,"numwedges=i"=>\$numwedges,"specieswedgetitle!"=>\$title_line_in_CCfile,"color_spp_by_wedge=i"=>\$colorsppbywedge,"color_bar_by_wedge=i"=>\$colorbarbywedge,"fontsize=f"=>\$fontsize);

if (!$ok || scalar @ARGV != 2 || $help) { 
	print "usage is perl $0 treefile output.svg\n";
	print "reads a treefile n newick format, parses and outputs a circular svg of it\n";
	print "can use __50HGI__ in place of treefile to use built-in 50 helminths output tree, or __PETE__ for Pete Olson't tapeworm\n";
	print "can use - in place of output.svg for stdout\n";
	print "options:\n";
	print "--help\tshow this\n";
	print "--hpix\theight of drawing in 'pixels' (default 800)\n";
	print "--wpix\twidth of drawing in 'pixels' (defaut 800)\n";
	print "--hideSPPnames\tdon't show labels on leaves\n";
	print "--space_for_labels=i\tspace for leaf labels and first wedges (default=100)\n";
	print "--start_angle=f\tstart angle in radians, measured from x axis (default 0)\n";
	print "--end_angle=f\tend angle in radians, measure from x axis (default 2*pi)\n";
	print "--multiply_by_pi\tconvenience function: whatever you specified for the above two options, multiply it by pi\n";
	print "--no_negative_branches\tset negative branches to 0\n";
	print "--root_branch_length\toverride root branch length in file, and set it to this value\n";
	print "\n\n";
	print "___OPTIONS RELATED TO COLORS AND WEDGES ETC___\n";
	print "--specieswedge=s\tturn on '50 helminths' style wedges on species names\n";
	print "\t\tfilename of 'colors_to_clades file with tab-delim lines, first column is species, next is a group factor, next is color\n";
	print "--numwedges=i\thow many pairs of 'factor/color' columns to read for specieswedge\n";
	print "--nospecieswedgetitle\tspecieswedgefile doesn't have a title line\n";
	print "--label_wedges=i\tlabel wedges with columns\n";
	print "--wedgelabelnudge=f\tfraction of wedge width to push labels out by\n";
	print "--wedgelabelsize=f\tsize of wedge labels (default=fontsize)\n";
	print "--fontsize=f\tsize of font\b";
	print "--color_spp_by_wedge=i\tuse color from wedge i to color species names instead of drawing wedges\n";
	print "--space_for_wedges=i\tspace for all subsequent wedges combined (default=20 if more than 1 wedge; otherwise 0)\n";
	print "--species_name_colors=s\tuse colors from file s to color species name text (tab-delim, two cols, first col is species names, second col is color hexadecmial RGB code, with hash)\n";
	print "--species_name_colors_column=i\twhich column of the species_name_colors file actually contains colors (default 1, 0 indexed)!\n";
	print "--also_color_dots\talso color dotted line the same as species label\n";
	print "\n\n";
	print "__OPTIONS RELATED TO BAR CHARTS__\n";
	print "--50helminthsbars=s\tspecify string with '50 helminths' style family membership string\n";
	print "--color_bar_by_wedge=i\tuse color from wedge i to color bars\n";
	print "--space_for_bars=i\tspace for bars(default=100 if 50helmnithbars specified, otherwise 0)\n";
	print "--drawlegend\tadd a legend for bar heights\n";
	print "--legendwidth=i\thow big a legend (default - make bars same size as plot, autosize legend to fit)\n";
	print "--legendheight=i\thow big a legend (default - make bars same size as plot, autosize legend to fit)\n";
	print "--highlight_bars=s\tcommand-delim list of species to outline bars for\n";
	print "\n\n";
	print "__OPTIONS RELATED TO ADMIXTURE PLOT__\n";
	print "--qmatrix=s\tspecify admixture .Q file for structure plot\n";
	print "--pedfile=s\tspecify pedfile input to admixture to get names\n";
	print "--qcolors=s\ttextfile with list of hex colors, one per line, same number of lines as cols in Q matrix (i.e. admixture's K parameter). If not specified, the first K of Kelly's 22 divergent colors are used!\n";
	print "--ldonglobal_fix_spp_names\tstrip trailing country name from tree labels\n";
	print "--annularsectors\tdraw annular sectors, not just bars\n";
	print "--annularsectorgap_proportion\tgap between adjacent sectors, as a proportion (0-1, default 0.2)\n";
	print "--min_drawable\tsmallest proportion to bother drawing (default 0.000011, as min. in L.don data is 0.00001)\n";
	print "--nomatchcolors\tdon't bother trying to match colors between concentric rings\n";
	print "--color_ref_ring=i\tuse ring i (0..n) as base for matching colors - only really makes sense to specify highest K. Default is to use the outer ring, and be concentric\n";
	print "--not_concentric\tif you don't specify color_ref_ring, the default is to use each ring in turn as reference forthe next.. unless you set this, then if uses *just* the outer ring as a reference\n";
	print "--multicolors\t do not make sure different clusters don't have the same color\n";
	print "--learncolors\t try to learn colors from species name colors\n";
	print "\n\n";
	print "__SPECIAL OPTIONS UNLIKELY TO MAKE SENSE TO NON-JAMES__\n";
	print "--50hgi_fix\treorder annoying node at base of flatworms, and suppress drawing outgroup!\n";
	print "--fix_spp_names_Ldon\tfudge to wrangle species names for James's L.donovani data\n";
	print "--skipSVGheader\tdon't include xml headers, SVG open or close tag: used to wrap this drawing within other SVG code..caution! you must provide these, and know something about drawing size (hpix and wpix)\n"; 
		exit();
}

my @KellyColorRGB = ("#fdfdfd", "#1d1d1d", "#ebce2b", "#702c8c", "#db6917", "#96cde6", "#ba1c30", "#c0bd7f", "#7f7e80", "#5fa641", "#d485b2", "#4277b6", "#df8461", "#463397", "#e1a11a", "#91218c", "#e8e948", "#7e1510", "#92ae31", "#6f340d", "#d32b1e", "#2b3514");

if ($labelwedges) { print STDERR "WARNING! - labelwedges uses textPath elements that don't render properly in some cases: e.g. cairoSVG. Chrome seems to handle them well!\n"; } 
my $outfile = $ARGV[1];

my $underscore_spaces = 0;

if ($learncolors && $species_name_colors eq "") { 
	die("\ndoesn't make sense to specify learncolors and not species_name_colors: i've nothing to learn from!\n");
}

if ($colorsppbywedge && $species_name_colors ne "") { 
	die("\ndoesn't make sense to specify both colorsppbywedge and species_name_colors\n");
}

if ($numwedges < 2) { $space_for_wedges = 0; }
if ($numwedges > -1 && $colors_to_clades eq "")  {
	die("\ndoesn't make sense to specify numwedges wihtout pointing toa specieswedge file\n");
}

if ($colorsppbywedge  && $colors_to_clades eq "")  {
	die("\ndoesn't make sense to specify colorsppbywedge wihtout pointing toa specieswedge file\n");
}

if ($colorbarbywedge  && $colors_to_clades eq "")  {
	die("\ndoesn't make sense to specify colorbarbywedge wihtout pointing toa specieswedge file\n");
}

if ($colorsppbywedge > 0 && $colorsppbywedge > $numwedges)  {
	die("\ndoesn't make sense to specify colorsppbywedge for a wedge bigger than numwedges!\n");
}

if ($colorbarbywedge > 0 && $colorbarbywedge > $numwedges)  {
	die("\ndoesn't make sense to specify colorbarbywedge for a wedge bigger than numwedges!\n");
}


if ($colors_to_clades ne "" && $numwedges == -1) { 
	print STDERR "specieswedge specified but numwedges not set. Assuming 1\n";
	$numwedges = 1;
}

if (scalar(@qmatrix > 0)  && $helminthfambars ne "") { die("can only draw one set of bars\n"); }
if (scalar(@qmatrix > 0) ne "" &&  $pedfile eq ""  ) { die("must specify Qmatrix and pedfile \n"); }

if ($helminthfambars eq "" && $drawlegend ) { die("can't drawlegend when there are no bars\n"); } 

if ($highlight_bars ) { 
	my @bar_bits = split(/\,/,$highlight_bars);
	print STDERR "read ".(scalar @bar_bits)." labels to highlight bars\n";
	foreach my $b (@bar_bits) { 
		(my $spp = $b) =~ s/_/ /g;
		$highlight_bars{$spp}++;
	}	
}


if ($multiply_by_pi) { 
	$start_angle = $start_angle * pi;
	$end_angle = $end_angle * pi;
	
}

#print STDERR "ENTRIES OF Q\n";
#foreach my $q (@qmatrix) { print STDERR $q."\n"; }
#array of hashes. each hash is species -> group for that level of group;
my @groups;
#array of hashes. each hash is group name -> color for that level of group.
#will store final color in file for that group name.
my @group_colors;

if ($colors_to_clades ne "" ) { 
	open(CC,"<",$colors_to_clades) or die("cannot open ".$colors_to_clades);
	while (<CC>) { 
		my $wanted = 1;
		if ($title_line_in_CCfile && $. == 1) { $wanted = 0; }
		
		if ($wanted) { 
			my @bits = split;
			if (scalar @bits < (1 + (2*$numwedges)) ) { 
				die("\nI don't seem to have enough tab-delim field on line ".$.." of file ".$colors_to_clades." for ".$numwedges." wedges\n");
			}
			my $sppname = shift @bits;
			$sppname =~ s/\_/ /g;
			foreach (my $i = 0; $i != $numwedges; $i++) { 
				my $group = shift @bits;
				my $colour = shift @bits;
				$groups[$i]{$sppname} = $group;
				$group_colors[$i]{$group} = $colour;
			}
		}
	} 
	close CC;
}

my %bar_data;
if ($helminthfambars ne "") { 
	my @bits = split(/\,/,$helminthfambars);
	foreach my $b (@bits) { 
	#	print STDERR $b."\n"; 
		$b =~ m/([a-zA-Z0-9\_]+)\(([0-9]+)\)/;
	#	print STDERR $1."->".$2."\n";
		my $count = $2;
		(my $spp = $1) =~ s/_/ /g;
	#	print STDERR $spp."->".$count."\n";
		$bar_data{$spp} = $count;
	}
	
	#die();
} elsif (scalar(@qmatrix == 0))  { 
	$space_for_bars = 0;
}


my $tstr = "";
my $pt = new PhyloTree;
if ($ARGV[0] eq "__50HGI__" ) { 
  $tstr = $HGItree; 
} elsif ($ARGV[0] eq "__PETE__") { 
	$tstr = $pete_tree;
} 
else  { 
	open(INF,"<$ARGV[0]") or die("cannot read ".$ARGV[0]);
	while(<INF>) { 
		chomp;
		$tstr .= $_;
	}
	close INF;
}


#print STDERR $tstr."\n";
$pt->parseNH($tstr);
#print STDERR "parsed!\n";
$pt->MakeLeafList();
$pt->MakeNodeList();
$pt->makeclusters();

if ($no_neg_branches) { 
	foreach my $n (@{$pt->{nodelist}}) { 
		if ($n->{edgelength} < 0 ) { $n->{edgelength} = 0; } 
	}
}
#foreach my $l (@{$pt->{leaflist}}) { 
#	print STDERR $l." ".$l->{leaf}." ".$l->{label}."\n";
#}
#die();

if ($hgi_fix) { 
	my $smed = -1;
	my $sman = -1;
	foreach my $n (@{$pt->{leaflist}}) { 
		if ($n->{label} eq "schmidtea mediterranea") { $smed = $n; }
		if ($n->{label} eq "schistosoma mansoni") {$sman = $n; }
	}
	
	my $root_of_flatties = $pt->getMRCA($smed,$sman);
	my $n1 = $root_of_flatties->{child};
	my $n2 = $root_of_flatties->{child}->{sibling};

	
	my $n2_currsib = $n2->{sibling};
	$root_of_flatties->{child} = $n2;
	$root_of_flatties->{child}->{sibling} = $n1;
	$n2->{sibling} = $n1;
	$n1->{sibling} = $n2_currsib;

	my $n3 = $n1->{child}->{child};
	my $n4 = $n1->{child}->{child}->{sibling};
	
	my $n4_currsib = $n4->{sibling};

	$n1->{child}->{child} = $n4;
	$n1->{child}->{child}->{sibling} = $n3;
	$n4->{sibling} = $n3;
	$n3->{sibling} = $n4_currsib;
	$root_branch_length = 0.05;
	
	$pt->MakeLeafList();
	$pt->MakeNodeList();
	$pt->makeclusters();

}

if (! $pt->{rootnode}->{edgelength} ) { $pt->{rootnode}->{edgelength} = $root_branch_length; } 

my %leaf_lookup;

my %node_angles;

my $number_of_leaves = $pt->{Nleaves};


#my $angle_per_leaf = (2 * pi) / $number_of_leaves;
my $angle_per_leaf = ($end_angle - $start_angle) / $number_of_leaves;

if (! $wedge_angle_fiddle ) { 
	$wedge_angle_fiddle = $angle_per_leaf / 2.5;
}

print STDERR "each leaf has ".$angle_per_leaf." radians\n";

my $curr_angle = $start_angle;
#find leaf angles
foreach my $n (@{$pt->{nodelist}}) {
    if ($n->{leaf} ) {
		$node_angles{$n} = $curr_angle;
		$curr_angle += $angle_per_leaf;
	}
}

#print STDERR "".(scalar keys %node_angles)." angles stored:\n";
#foreach my $n (keys %node_angles) { print STDERR $n."->".$node_angles{$n}."\n"; }
#find internal node angles - mean of descendents
visit_postorder($pt->{rootnode},\%node_angles);
#print STDERR "".(scalar keys %node_angles)." angles stored:\n";
#foreach my $n (keys %node_angles) { print STDERR $n."->".$node_angles{$n}."\n"; }
my $diagram_base_diameter = min($wpix,$hpix) - (2.0 *($space_for_labels+$space_for_bars+$space_for_wedges));
print STDERR "regions are: bars=".$space_for_bars.",wedges=".$space_for_wedges.",labels=".$space_for_labels.",tree_diam=".$diagram_base_diameter."\n";
my $diagram_base_radius = $diagram_base_diameter / 2.0;

#get all root-tip lengths
my $max_depth = -1;
foreach my $n (@{$pt->{nodelist}}) {
    if ($n->{leaf} ) {
    	my $p = $n;
		my $pathlen = 0;
    	while ($p) { 
    		if (exists $p->{edgelength} ) { 
    			$pathlen += $p->{edgelength};
    		}
    		$p = $p->{anc};
    	}
 #   	print STDERR $n->{label}."->".$pathlen."\n";
		if ($pathlen > $max_depth) { $max_depth = $pathlen; }
	}
	
}

my %species_name_colors;
if ($species_name_colors ne "") { 
	open(SNC,"<",$species_name_colors) or die("cannot open file ".$species_name_colors." to read species name <-> color map\n");
	while(<SNC>) { 
		chomp;
		my @bits = split(/\s+/);
		print "'".$bits[0]."'->'".$bits[$species_name_colors_column]."'\n";
		$species_name_colors{$bits[0]} = $bits[$species_name_colors_column];
	}
	close SNC;
	print STDERR "read colors for ".(scalar keys %species_name_colors)." unique species names\n";
}

my $radius_per_brlen = $diagram_base_radius / $max_depth;
print STDERR "got ".$radius_per_brlen." mm per unit br. length\n";

my %branch_start_radius ;
my %branch_end_radius;
my %arc_start;
my %arc_end;
node_placement_preorder($pt->{rootnode},\%branch_start_radius,\%branch_end_radius,\%arc_start,\%arc_end,$radius_per_brlen,\%node_angles);

my $svgfileglob;

if ($outfile eq "-") { 
	$svgfileglob = *STDOUT;  
} else {
open($svgfileglob ,">",$outfile) or die("cannot write to SVG file ".$outfile);
}
if (! $skipSVGheader) { 
	print $svgfileglob  "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
	print $svgfileglob  '<!-- Generator: JAMES\'s tree_to_SVGcircletree.pl  -->'."\n";
	print $svgfileglob  '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">'."\n";
	print $svgfileglob  '<svg version="1.1" id="Layer_1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0mm" y="0mm" width="'.$wpix.'" height="'.$hpix.'" viewBox="0 0 '.$wpix.' '.$hpix.'" enable-background="new 0 0 '.$wpix.' '.$hpix.'" xml:space="preserve">'."\n";
}



#---------------------------------------------------------
#print some wedges!
if ($colors_to_clades ne "") { 
	#HOW BIG - first wedge covers species names, second and third are thin
	#array of $numwedges entries with diameters of each wedge
	my @wedge_widths;
	my @wedge_starts;
	my @wedge_label_radius;
	
	$wedge_starts[0] = $diagram_base_radius;
	$wedge_widths[0] = $space_for_labels;
	
	if ($numwedges > 1 ) { 
		my $space_per_wedge = $space_for_wedges / ($numwedges-1);
		my $wedge_start_here = $diagram_base_radius + $space_for_labels +$gap_between_wedges;
	
		foreach (my $i = 1; $i != $numwedges; $i++) { 
			$wedge_starts[$i] = $wedge_start_here;
			$wedge_widths[$i] = $space_per_wedge-$gap_between_wedges;
			#place labels, if any!
			$wedge_label_radius[$i] = $wedge_start_here + ($wedgelabelnudge * ($space_per_wedge-$gap_between_wedges));
			$wedge_start_here += $space_per_wedge;	
		}
	}
	#find angles of wedges for each piece
	my $labelhash = return_label_hash($pt);
	
	foreach (my $i = 0; $i != $numwedges; $i++ ) { 
		my %groups_to_members;
		my %nodes_to_group;
		
		foreach my $spp ( keys %{$groups[$i]} ){
			my $node = ${$labelhash}{$spp};
			my $group = ${$groups[$i]}{$spp};
			$nodes_to_group{$node} = $group;
			if (exists $groups_to_members{ $group } ) { 
				push( @{$groups_to_members{$group}},$node);
			} else { $groups_to_members{$group} = [$node]; }  
		}
		
		print STDERR "for wedge ".$i." got ".(scalar keys %groups_to_members)." groups\n";
		
		#new style is merging adjacent blocks
		#old style alternative is now commented out below..
		#do I need this?
		my %angles_to_nodes;
		foreach my $k ( keys %node_angles) { 
			$angles_to_nodes{$node_angles{$k}} = $k;
		}
		
		#each wedge identified by the minimal leaf pointer
		
#		my %min_angle_wedge;
#		my %max_angle_wedge;
#	 	my %wedge_group;
	 	my %node_labels;
	 	
		foreach my $k (@{$pt->{leaflist}}) { 
		#	print STDERR $k." ".$k->{leaf}." ".$k->{label}."\n";
		
#			$min_angle_wedge{$k} = $node_angles{$k};
#			$max_angle_wedge{$k} = $node_angles{$k};
		#	if (exists $nodes_to_group{$k} ) { 
#				$wedge_group{$k} = $nodes_to_group{$k};
		#	} else { die("can't find group for leaf ".$k." (".$k->{label}.")"); }
			$node_labels{$k} = $k->{label};
		}
		#die();
		my @leaf_nodes_sorted_by_angle = sort {  $node_angles{$a} cmp $node_angles{$b} } @{$pt->{leaflist}};
		
		#print STDERR join(",",@leaf_nodes_sorted_by_angle)."\n";
		my @output;
		#foreach my $n () { 
		foreach my $n (@leaf_nodes_sorted_by_angle) { 
			my $lab = $node_labels{$n};
			my $angle = $node_angles{$n};
			my $group = $nodes_to_group{$n};
			
			push(@output,$lab.",".$n.",".$angle.",".$group); 
		}
	#	print STDERR join("; ",@output)."\n";
		my $start_i = 0; 
		my $last_group = $nodes_to_group{$leaf_nodes_sorted_by_angle[$start_i]};
		#first_walk_backwards!
		my $walk_i = (scalar @leaf_nodes_sorted_by_angle) - 1;
		my $group_here = $nodes_to_group{$leaf_nodes_sorted_by_angle[$walk_i]};
		while($group_here eq $last_group) { 
		#	print STDERR $walk_i." ".$group_here." ".$last_group." ".$start_i."\n";
			$start_i = $walk_i;
			$walk_i--;
			$group_here = $nodes_to_group{$leaf_nodes_sorted_by_angle[$walk_i]};
			
		}
		if ($start_i != 0 ) { 
			print STDERR "first group starts before circle origin at ".$start_i."\n";
		} else { 
			print STDERR "first group starts at circle origin at ".$start_i."\n";
		}
	#	print STDERR "now walking forwards!\n";
		my $first_group_start = $start_i;
		$walk_i = 0;
		
		my $group_marker = 0;
	
		
		my @group_starts;
		my @group_ends;
		my @group_labels;
		
		$group_starts[$group_marker] = $first_group_start;
		
		my $group_max = $first_group_start;
		#if first group starts AT origin, need to iterate all the way round
		if ($group_max == 0 ) { $group_max = scalar @leaf_nodes_sorted_by_angle; } 	
		while($walk_i < $group_max) { 
			my $group_here = $nodes_to_group{$leaf_nodes_sorted_by_angle[$walk_i]};
	#		print STDERR "at ".$walk_i." in ".$group_here.", previously. ".$last_group."\n";
			while($group_here eq $last_group && $walk_i < $group_max) { 
		#		print STDERR $walk_i." ".$group_here." ".$last_group." ".$start_i."\n";
				$start_i = $walk_i;
				$walk_i++;
				$group_here = $nodes_to_group{$leaf_nodes_sorted_by_angle[$walk_i]};
			}
			#save_group_details;
			$group_ends[$group_marker] = $start_i;
			$group_labels[$group_marker] = $last_group;
			if ($walk_i < $group_max ) {
				$group_starts[$group_marker+1] = $start_i+1;
				
			} 	
	#		print STDERR "found group end at ".$start_i." for group ".$group_marker."\n";
	#		print STDERR "next group starts at ".($start_i+1)." for group ".($group_marker+1)."\n";
			$group_marker++;
			
			$start_i = $walk_i;
			$walk_i++;
			$last_group = $group_here;
		}
	#	print STDERR "got ".(scalar @group_starts)." groups\n";
		foreach (my $gri = 0; $gri != scalar @group_starts; $gri++) { 
#			print STDERR "group ".$gri." (".$group_labels[$gri].") is from ".$group_starts[$gri]."->".$group_ends[$gri]."\n";
			my $node_name = $node_labels{$leaf_nodes_sorted_by_angle[$group_starts[$gri]]};
			#wedge_angles max and min
			my $angle1 = $node_angles{$leaf_nodes_sorted_by_angle[$group_starts[$gri]]};
			my $angle2 = $node_angles{$leaf_nodes_sorted_by_angle[$group_ends[$gri]]};
			my $color = $group_colors[$i]{$groups[$i]{$node_name}};
			my $min = min($angle1,$angle2);
			my $max = max($angle1,$angle2);
			my $group = $groups[$i]{$node_name};
	#		print STDERR "group ".$gri." is called ".$group."\n";
			
			if ($group ne "Outgroup" || !$hgi_fix) {
	#			print STDERR "set of ".(scalar @{$groups_to_members{$k}})."\n";
	#
	#			my $mrca = findMRCAofmanyleaves($groups_to_members{$k},$pt);
	#			print STDERR "GOT MRCA=".$mrca."\n";
	#			my ($min, $max) = returnMinMaxAnglesFromMRCA($mrca,\%node_angles);
	#			print STDERR "min_angle = ".$min." max=".$max."\n";

				my $large_arc = $max - $min <= pi ? 0 : 1;
				my $sweep_swap = 0;
			
			
				if ($group_starts[$gri] > $group_ends[$gri]) { 
					$large_arc = abs(1-$large_arc);
					$sweep_swap = 1;	
			#		print STDERR "SWEEP SWAPPED FOR ".$gri."\n";	
					$min += $wedge_angle_fiddle;
					$max -= $wedge_angle_fiddle;
				}	else { 
					$min -= $wedge_angle_fiddle;
					$max += $wedge_angle_fiddle;
				}
			
	#			print STDERR "my outer radius is ".($wedge_starts[$i] +$wedge_widths[$i] )."\n";
				my $arc1_startx = $wedge_starts[$i] * cos($min);
				my $arc1_starty =  $wedge_starts[$i] * sin($min);
				my $arc1_endx  =  $wedge_starts[$i] * cos($max);
				my $arc1_endy =  $wedge_starts[$i] * sin($max);
				my $arc2_startx = ($wedge_starts[$i] + $wedge_widths[$i] ) * cos($min);
				my $arc2_starty =  ($wedge_starts[$i] + $wedge_widths[$i])* sin($min);
				my $arc2_endx  =  ($wedge_starts[$i] +$wedge_widths[$i] )* cos($max);
				my $arc2_endy =  ($wedge_starts[$i] + $wedge_widths[$i]) * sin($max);
				print $svgfileglob  "<path d=\"M".($diagram_base_radius+$space_for_wedges + $space_for_bars+ $space_for_labels + $arc1_startx).",".($diagram_base_radius +$space_for_wedges+ $space_for_bars+ $space_for_labels + $arc1_starty);
				if ($sweep_swap) { 
					print $svgfileglob  " A".$wedge_starts[$i].",".$wedge_starts[$i]." 0 ".$large_arc.",0 ";
				} else{ 
					print $svgfileglob  " A".$wedge_starts[$i].",".$wedge_starts[$i]." 0 ".$large_arc.",1 ";
				}
				print $svgfileglob  " ".($diagram_base_radius +$space_for_wedges + $space_for_bars + $space_for_labels + $arc1_endx).",".($diagram_base_radius +$space_for_wedges+ $space_for_bars+ $space_for_labels + $arc1_endy);
				print $svgfileglob  " L".($diagram_base_radius +$space_for_wedges+ $space_for_bars+ $space_for_labels + $arc2_endx)." ".($diagram_base_radius +$space_for_wedges+ $space_for_bars+ $space_for_labels + $arc2_endy)." ";
				if ($sweep_swap) { 
					print $svgfileglob  " A".($wedge_starts[$i] + $wedge_widths[$i]) .",".($wedge_starts[$i] + $wedge_widths[$i]) ." 0 ".$large_arc.",1 ";
				} else{ 
					print $svgfileglob  " A".($wedge_starts[$i] + $wedge_widths[$i]) .",".($wedge_starts[$i] + $wedge_widths[$i]) ." 0 ".$large_arc.",0 ";
				}
				print $svgfileglob  " ".($diagram_base_radius +$space_for_wedges+ $space_for_bars + $space_for_labels + $arc2_startx).",".($diagram_base_radius +$space_for_wedges+ $space_for_bars+ $space_for_labels + $arc2_starty);
				print $svgfileglob  " L".($diagram_base_radius +$space_for_wedges + $space_for_bars+ $space_for_labels + $arc1_startx)." ".($diagram_base_radius +$space_for_wedges+ $space_for_bars+ $space_for_labels + $arc1_starty)." ";
				print $svgfileglob  "\" style=\"stroke:none; fill:".$color.";\" fill-opacity=\"0.5\"/>\n";
				
				if ($i > 0 && $labelwedges) { 
					#write labels for each wedge
					#first make arc at wedge_label_radius
					my $label_arc_startx = $wedge_label_radius[$i] * cos($min);
					my $label_arc_starty =  $wedge_label_radius[$i] * sin($min);
					my $label_arc_endx  =  $wedge_label_radius[$i] * cos($max);
					my $label_arc_endy =  $wedge_label_radius[$i] * sin($max);
					print $svgfileglob  "<path id=\"curve_".$i."_".$gri."\" d=\"M".($diagram_base_radius+$space_for_wedges + $space_for_bars+ $space_for_labels + $label_arc_startx).",".($diagram_base_radius +$space_for_wedges+ $space_for_bars+ $space_for_labels + $label_arc_starty);
					if ($sweep_swap) { 
						print $svgfileglob  " A".$wedge_label_radius[$i].",".$wedge_label_radius[$i]." 0 ".$large_arc.",0 ";
					} else{ 
						print $svgfileglob  " A".$wedge_label_radius[$i].",".$wedge_label_radius[$i]." 0 ".$large_arc.",1 ";
					}
					print $svgfileglob  " ".($diagram_base_radius +$space_for_wedges + $space_for_bars + $space_for_labels + $label_arc_endx).",".($diagram_base_radius +$space_for_wedges+ $space_for_bars+ $space_for_labels + $label_arc_endy);
					print $svgfileglob "\" style=\"stroke:transparent; fill:transparent;\" />\n";
					my $fontstyle = "normal";
					my $arclen = 2 * pi * (( $max - $min)  / (2 * pi)) * $wedge_label_radius[$i];
					#the 9 / 16 voodoo here is guess from http://www.upsdell.com/BrowserNews/res_fontmetrics.htm#Arial_400 !
					my $textlen = (9/ 16) * $wedgelabelsize * length($group_labels[$gri]);
					my $start_fraction = (0.5 - (($textlen ) / (2  * $arclen))) * 100;
					print STDERR "arclen=".$arclen." textlen=".$textlen."\n";
					print STDERR  " starting label at offset ".$start_fraction."\n";
					print $svgfileglob "<text alignment-baseline=\"baseline\"  font-family=\"arial\" font-style=\"".$fontstyle."\" font-size=\"".$wedgelabelsize."\" fill=\"".$textcolor."\"><textPath xlink:href=\"#curve_".$i."_".$gri."\" startOffset=\"".$start_fraction."%\">".$group_labels[$gri]."</textPath></text>\n";
				}
			
			}
		
		}
		
	  
	}	
}


#-------------------------
# THIS BIT DRAWS THE ACTUAL TREE!!!
#-------------------------
#---------------------------------------------------------
#position and print tree branches and labels
#
my $max_radius = -1;
#
#
foreach my $n (@{$pt->{nodelist}}) {
#	print STDERR "node angle ".$n->{label}."->".$node_angles{$n}."\n";
	my $angle = $node_angles{$n};
	my $cosa = cos($angle);
	my $sina = sin($angle);
	#my $height = $genome_size{$spp} * $mm_per_bp;
	
	#draw text labels
	#draw dotted lines to leaves
	if (! $hideSPPnames) { 
		if ($n->{leaf}) { 
			my $label = $n->{label};
			if ($hgi_fix) { 
		#	print STDERR $label."\n";
				$label = ucfirst($label);
				$label =~ s/kr3021/sp./;
		#		print STDERR $label."\n";
			}
			my $transx = ($diagram_base_radius + $text_padding) * $cosa;
			my $transy = ($diagram_base_radius + $text_padding) * $sina;

			print $svgfileglob  "<g transform=\"matrix(".$cosa.",".$sina.",".(-1.0 * $sina).",".$cosa.",".($diagram_base_radius + $space_for_labels +$space_for_wedges+ $space_for_bars + $transx).",".($diagram_base_radius + $space_for_labels +$space_for_wedges+ $space_for_bars+ $transy).")\">\n";
			my $vertical_flip  = 0;
			if ($node_angles{$n} > 0.5 * pi && $node_angles{$n} < 1.5 * pi)  {
				$vertical_flip = 1;
				print STDERR "flipping ".$n->{label}."\n";
				print $svgfileglob "<g transform=\"rotate (180)\">\n";
			}
		
			if ($colorsppbywedge)  {
				my $group = $groups[$colorsppbywedge-1]{$n->{label}};
				my $col = $group_colors[$colorsppbywedge-1]{$group};
			#	print STDERR "'".$n->{label}."' is in '".$group."', and so will be '".$col."'\n";
				$textcolor = $col;
			} elsif (scalar keys %species_name_colors > 0) {
				my $new_label = $label;
				if ($fix_spp_names_Ldon) { 
					#my @bits = split(" ",$label);
                                	#pop @bits;
                                	#$new_label = join(" ",@bits);
                                	($new_label = $label) =~ s/\s/_/g;
				}
				if (exists $species_name_colors{$new_label} ) { 
					$textcolor = $species_name_colors{$new_label};
				} else  {
					print STDERR "warning - no colors found for species ".$label." ('".$new_label."')\n";
				}
				
				
			}
		#	if (0) { 	
			my $fontstyle = "normal";
			
			if ($hgi_fix) { $fontstyle = "italic"; }
			if (length($label) > (1.8*(($space_for_labels - $text_padding) / $fontsize ))) { 
				#squish text to fit space
				if ($vertical_flip) { 
					print $svgfileglob  "<text alignment-baseline=\"middle\" text-anchor=\"end\" x=\"0\" y=\"0\"  font-family=\"arial\" font-style=\"".$fontstyle."\" font-size=\"".$fontsize."\" fill=\"".$textcolor."\" textLength=\"".($space_for_labels - $text_padding)."\" lengthAdjust=\"spacingAndGlyphs\">".$label."</text>\n";
				} else { 
					print $svgfileglob  "<text alignment-baseline=\"middle\" text-anchor=\"start\" x=\"0\" y=\"0\"  font-family=\"arial\" font-style=\"".$fontstyle."\" font-size=\"".$fontsize."\" fill=\"".$textcolor."\" textLength=\"".($space_for_labels - $text_padding)."\" lengthAdjust=\"spacingAndGlyphs\">".$label."</text>\n";
				}
			} else {  
				if ($vertical_flip) { 
					print $svgfileglob  "<text alignment-baseline=\"middle\" text-anchor=\"end\" x=\"0\" y=\"0\"  font-family=\"arial\" font-style=\"".$fontstyle."\" font-size=\"".$fontsize."\" fill=\"".$textcolor."\">".$label."</text>\n";
				} else {
					print $svgfileglob  "<text alignment-baseline=\"middle\" text-anchor=\"start\" x=\"0\" y=\"0\"  font-family=\"arial\" font-style=\"".$fontstyle."\" font-size=\"".$fontsize."\" fill=\"".$textcolor."\">".$label."</text>\n";
				}
			}
			if ($vertical_flip) { print $svgfileglob "</g>\n"; }
			print $svgfileglob  "</g>\n";

		}
	}
	#draw straight branches
	
	my $start_radius = $branch_start_radius{$n};
	my $end_radius = $branch_end_radius{$n};
	if ($end_radius > $max_radius) { $max_radius = $end_radius; }
#	print STDERR "Node ".$n." (";
#	if ($n->{leaf}) { print STDERR $n->{label}; } else { print STDERR "INTERNAL"; }
#	print STDERR ") angle=".$angle." ".$start_radius."-".$end_radius;
	my $arc_start_angle  = -1;
	my $arc_end_angle  = -1;
	if (!$n->{leaf}) { 
		$arc_start_angle = $arc_start{$n};
		$arc_end_angle = $arc_end{$n};
#		print STDERR " ".$arc_start_angle."-".$arc_end_angle."\n"; 
	}
	my $startx = $start_radius* $cosa;
	my $starty = $start_radius * $sina;
	my $endx = $end_radius * $cosa;
	my $endy = $end_radius * $sina;
	print $svgfileglob  "<line x1=\"".($diagram_base_radius + $space_for_labels +$space_for_wedges+ $space_for_bars+ $startx)."\" y1=\"".($diagram_base_radius+$space_for_wedges + $space_for_bars+ $space_for_labels + $starty)."\" x2=\"".($diagram_base_radius +$space_for_wedges+ $space_for_bars + $space_for_labels + $endx)."\" y2=\"".($diagram_base_radius +$space_for_wedges+ $space_for_bars + $space_for_labels +$endy)."\" style=\"stroke:rgb(0,0,0);stroke-width:1\" />\n";
	if ($n->{leaf} ) { 
		my $even_more_endx = $diagram_base_radius * $cosa;
		my $even_more_endy = $diagram_base_radius * $sina;	
		my $line_stroke_color = "rgb(128,128,128)";
		if ($also_color_dots ) { $line_stroke_color = $textcolor; }
		print $svgfileglob  "<line x1=\"".($diagram_base_radius + $space_for_labels +$space_for_wedges+ $space_for_bars + $endx)."\" y1=\"".($diagram_base_radius+$space_for_wedges+ $space_for_bars + $space_for_labels +$endy)."\" x2=\"".($diagram_base_radius+$space_for_wedges + $space_for_bars + $space_for_labels + $even_more_endx)."\" y2=\"".($diagram_base_radius +$space_for_wedges+ $space_for_bars + $space_for_labels +$even_more_endy)."\" stroke-dasharray=\"3, 3\"  style=\"stroke:".$line_stroke_color.";stroke-width:0.75\" />\n";
		
	}
	#draw arc
	if (!$n->{leaf} ) { 
		my $large_arc = $arc_end_angle - $arc_start_angle <= pi ? 0 : 1;
		my $arc_startx = $end_radius * cos($arc_start_angle);
		my $arc_starty = $end_radius * sin($arc_start_angle);
		my $arc_endx  = $end_radius * cos($arc_end_angle);
		my $arc_endy = $end_radius * sin($arc_end_angle);
		print $svgfileglob  "<path d=\"M".($diagram_base_radius+$space_for_wedges + $space_for_bars+ $space_for_labels + $arc_startx).",".($diagram_base_radius +$space_for_wedges+ $space_for_bars+ $space_for_labels + $arc_starty);
		print $svgfileglob  " A".$end_radius.",".$end_radius." 0 ".$large_arc.",1 "; 
		print $svgfileglob  " ".($diagram_base_radius+ $space_for_bars +$space_for_wedges+ $space_for_labels + $arc_endx).",".($diagram_base_radius +$space_for_wedges+ $space_for_bars+ $space_for_labels + $arc_endy)."\" style=\"stroke: #000000; stroke-width:1; fill:none;\"/>\n";
	}
	
}
#
#
print STDERR "max radius is ".$max_radius."\n";

#------------------------------------------------------------
#---------------------------------------------------------
#______________
#print some bars

my $grey_median = 1;

if ($helminthfambars ne "") { 
	my $max = max (values %bar_data);
	my $dist_per_unit = ($space_for_bars-$bar_padding) / $max; 
	print STDERR "longest bar is ".$max." so we have ".$dist_per_unit." dist/unit\n";
	foreach my $b (keys %bar_data){ 
		my $bar_height = ( $bar_data{$b} * $dist_per_unit ) / 1.0;
		my $labelhash = return_label_hash($pt);
		my $bar_angle = $node_angles{${$labelhash}{$b}};
	#	print STDERR "bar of ".$bar_height." at ".$bar_angle."\n";
		my $cosa = cos($bar_angle);
		my $sina = sin($bar_angle);
		#my $height = $genome_size{$spp} * $mm_per_bp;
	
		my $transx = ($diagram_base_radius + $space_for_wedges+ $space_for_labels + $bar_padding) * $cosa;
		my $transy =  ($diagram_base_radius + $space_for_wedges+ $space_for_labels + $bar_padding) * $sina;

		print $svgfileglob  "<g transform=\"matrix(".$cosa.",".$sina.",".(-1.0 * $sina).",".$cosa.",".($diagram_base_radius + $space_for_labels +$space_for_wedges+ $space_for_bars + $transx).",".($diagram_base_radius + $space_for_labels +$space_for_wedges+ $space_for_bars+ $transy).")\">\n";
		my $barcolor = "#000000";
		if ($colorbarbywedge)  {
			my $group = $groups[$colorbarbywedge-1]{$b};
			my $col = $group_colors[$colorbarbywedge-1]{$group};
		#	print STDERR "'".$n->{label}."' is in '".$group."', and so will be '".$col."'\n";
			$barcolor = $col;
		}
#		print STDERR "bar:".${$labelhash}{$b}." ".$b."\n";
		if (exists $highlight_bars{ $b  } ) { 
			#print STDERR "adding highlight to bar ".${$labelhash}{$b}."\n";
			print $svgfileglob   "<rect x=\"".(0)."\" y=\"".(-0.5 * $bar_width)."\" height=\"".$bar_width."\" width=\"".$bar_height."\"   stroke-dasharray=\"3, 3\" style=\"stroke:black;stroke-width:1;fill:".$barcolor."\"/>\n";
		} else { 
		print $svgfileglob   "<rect x=\"".(0)."\" y=\"".(-0.5 * $bar_width)."\" height=\"".$bar_width."\" width=\"".$bar_height."\"  style=\"fill:".$barcolor."\"/>\n";
		}
		print $svgfileglob  "</g>\n";
	}
	if ($drawlegend) { 
		#choose a quadrant;
		my $q1_start = pi / 8;
		my $q1_end = 3 * pi / 8;
		my $q2_start = 5 * pi / 8;
		my $q2_end = 7 * pi / 8;
		my $q3_start = 9 * pi / 8;
		my $q3_end = 11 * pi / 8;
		my $q4_start = 13 * pi / 8;
		my $q4_end = 15 * pi / 8;
		
		my @q1 ;
		my @q2 ;
		my @q3 ;
		my @q4 ;
		
		foreach my $n (@{$pt->{nodelist}}) {
		    if ($n->{leaf} ) {
		    	my $count = 0;
		    	if ($bar_data{$n->{label}}) { 
		    		$count = $bar_data{$n->{label}};
				}
				if ( $node_angles{$n} >= $q1_start && $node_angles{$n} <= $q1_end) { push(@q1,$count); } 
				if ( $node_angles{$n} >= $q2_start && $node_angles{$n} <= $q2_end) { push(@q2,$count); } 
				if ( $node_angles{$n} >= $q3_start && $node_angles{$n} <= $q3_end) { push(@q3,$count); } 
				if ( $node_angles{$n} >= $q4_start && $node_angles{$n} <= $q4_end) { push(@q4,$count); } 

			}
		}
		#print STDERR "Qs are (".join(",",@q1)."),(".join(",",@q2)."),(".join(",",@q3)."),(".join(",",@q4).")\n";
		my @q_maxs = (max(@q1),max(@q2),max(@q3),max(@q4));
		print STDERR "max per q=".join(",",@q_maxs)."\n";
		my $overall_min = min @q_maxs;
		my $quadrant_i = 0;
		while ($q_maxs[$quadrant_i] != $overall_min) {
			$quadrant_i++;
		}
		print STDERR "smallest quadrant is ".$quadrant_i."\n";
		
		#find max, min and median
		my @values_to_show;
		my @value_labels;
		my $max_bar = max(values %bar_data);
		my $min_bar = min(values %bar_data);
		my $median_bar = median(values %bar_data);
		my $which_max;
		my $which_min;
		my $which_median;
		my @bar_colors;
		
		foreach my $b (keys %bar_data) { 
			if ($bar_data{$b} == $max_bar) { $which_max = $b;}
			if ($bar_data{$b} == $min_bar) { $which_min = $b;}
			if ($bar_data{$b} == floor($median_bar) ) { $which_median = $b;}
			
		}
		if ($max_bar == $min_bar ) { 
			$values_to_show[0] = $max_bar;
			$value_labels[0] = "value";
			if ($colorbarbywedge) { 
				my $group = $groups[$colorbarbywedge-1]{$which_max};
				my $col = $group_colors[$colorbarbywedge-1]{$group};
				$bar_colors[0] = $col;
			} else { 
				$bar_colors[0] = "#9E9E9E";
			}
		} elsif ( $max_bar == $median_bar ) { 
			$values_to_show[0] = $max_bar;
			$values_to_show[1] = $min_bar;
			$value_labels[0] = "median/max";
			$value_labels[1] = "min";
			if ($colorbarbywedge) { 
				my $group = $groups[$colorbarbywedge-1]{$which_max};
				my $col = $group_colors[$colorbarbywedge-1]{$group};
				$bar_colors[0] = $col;
				$group = $groups[$colorbarbywedge-1]{$which_min};
				$col = $group_colors[$colorbarbywedge-1]{$group};
				$bar_colors[1] = $col;
			} else { 
				$bar_colors[0] = "#9E9E9E";
				$bar_colors[1] = "#9E9E9E";

			}
		} elsif ($min_bar == $median_bar  ) { 
			$values_to_show[0] = $max_bar;
			$values_to_show[1] = $min_bar;
			$value_labels[0] = "max";
			$value_labels[1] = "median/min";
			if ($colorbarbywedge) { 
				my $group = $groups[$colorbarbywedge-1]{$which_max};
				my $col = $group_colors[$colorbarbywedge-1]{$group};
				$bar_colors[0] = $col;
				$group = $groups[$colorbarbywedge-1]{$which_min};
				$col = $group_colors[$colorbarbywedge-1]{$group};
				$bar_colors[1] = $col;

			} else { 
				$bar_colors[0] = "#9E9E9E";
				$bar_colors[1] = "#9E9E9E";

			}
		} else {
			$values_to_show[0] = $max_bar;
			$values_to_show[1] = $median_bar;
			$values_to_show[2] = $min_bar;
			$value_labels[0] = "max";
			$value_labels[1] = "median";
			$value_labels[2] = "min";
			if ($colorbarbywedge) { 
				my $group = $groups[$colorbarbywedge-1]{$which_max};
				my $col = $group_colors[$colorbarbywedge-1]{$group};
				$bar_colors[0] = $col;
				$group = $groups[$colorbarbywedge-1]{$which_median};
				$col = $group_colors[$colorbarbywedge-1]{$group};
				$bar_colors[1] = $col;
				if ($grey_median) { $bar_colors[1] = "#000000"; }
				
				$group = $groups[$colorbarbywedge-1]{$which_min};
				$col = $group_colors[$colorbarbywedge-1]{$group};
				$bar_colors[2] = $col;

			} else { 
				$bar_colors[0] = "#9E9E9E";
				$bar_colors[1] = "#9E9E9E";
				$bar_colors[2] = "#9E9E9E";
			}
		}
		#generate SVG for plot
		print STDERR "colors are: ".join(",",@bar_colors)."\n";
		#my $code_for_barplot =		BarPlotSVG(\@values_to_show,\@value_labels,\@bar_colors,$legend_width,$legend_height);
		#print STDERR "bar dimensions are: ".$dist_per_unit."pixel/unit ".$bar_width."px bar widths\n";
		my $code_for_barplot = "HELLO!";
		if ($legend_width && $legend_height) { 
			$code_for_barplot =		BarPlotSVG(\@values_to_show,\@value_labels,\@bar_colors,$legend_width,$legend_height);
		} else { 
		 ($code_for_barplot,$legend_width,$legend_height) =		BarPlotSVG_flex(\@values_to_show,\@value_labels,\@bar_colors,$dist_per_unit,$bar_width);
		}
		print STDERR "legend will be ".$legend_width."x".$legend_height."\n";
		if ($quadrant_i == 0) { 
			print $svgfileglob "<g transform=\"translate(".($wpix-($legend_width + $legend_margin) ) .",".($hpix-($legend_height + $legend_margin) ).")\">\n";		
		} elsif ($quadrant_i == 1) { 
			print $svgfileglob "<g transform=\"translate(".$legend_margin.",".($hpix-($legend_height + $legend_margin) ).")\">\n";		

		} elsif ($quadrant_i == 2) { 
			print $svgfileglob "<g transform=\"translate(".$legend_margin.",".$legend_margin.")\">\n";
		} elsif ($quadrant_i == 3) { 
			print $svgfileglob "<g transform=\"translate(".($wpix-($legend_width + $legend_margin) ) .",".$legend_margin.")\">\n";		

		} else { die("BROKEN QUADRANTS ".$quadrant_i);}
		print $svgfileglob $code_for_barplot;
		print $svgfileglob "</g>\n";
	}

}

#---------------------------------------------------------
#______________
#print ADMIXTURE qfile data

my %sample_name_to_qrow;
if (scalar(@qmatrix > 0) ) { 
open(PED,"<",$pedfile) or die("cannot open PED file ".$pedfile);
while (<PED>) { 
	chomp;
	my @bits = split(/\s+/);
	$sample_name_to_qrow{$.} = $bits[1];
}
print STDERR "read ".(keys %sample_name_to_qrow)." rows from PED\n";
close PED;
}
#array of hashes
my @all_qdata;
#just an array
my @all_matrix_widths;
#array of arrays
my @allcolors;

for(my $qmi = 0; $qmi != scalar @qmatrix; $qmi++) { 
	my $qmatrix = $qmatrix[$qmi];

	my $matrix_width = -1;
	my %qdata;
	open(QDATA,"<",$qmatrix) or die("cannot open qfile '".$qmatrix."'\n"); 
	while(<QDATA>) { 
		chomp;
		my @bits = split(/\s+/);
		if ($matrix_width > -1 ) { 
			if ($matrix_width != scalar @bits) { die("problem - Q matrix width not consistent"); }
		} else { 
			$matrix_width = scalar @bits;
		}
		my $name = $sample_name_to_qrow{$.};
		$qdata{$name} = \@bits;
	}
	close QDATA;
	push(@all_matrix_widths,$matrix_width);
	push(@all_qdata,\%qdata);
	my @qcols;
	if ($qcolors ne "") { 
		
	
		open(QCOLS,"<",$qcolors) or die("cannot open Q colors file ".$qcolors);
		while(<QCOLS>) { 
			chomp;
			push(@qcols,$_);
		}
		close QCOLS;
	} else {
		for (my $i = 1; $i != (1 + $matrix_width); $i++) { 
			push(@qcols,$KellyColorRGB[$i]);
		}
	}
	print STDERR "read ".(scalar @qcols)." q-matrix colors\n";
	push(@allcolors,\@qcols);
}

#differe 'learn colors' algorithm - 0 means iterate through q values, and find best color in order..
#but if 'small' q groups go first, colors can look funny.
#new algorithm finds maximum global weight, bias towards big, consistently colored groups
#assign colors to biggest weight groups first..
#inefficient (for q groups and n colors, can take q*q*n
my $flip_learn_colors = 1;

if ($learncolors && $species_name_colors ne "") { 
	my @backup_cols;
	if ($qcolors ne "") { 
		
	
		open(QCOLS,"<",$qcolors) or die("cannot open Q colors file ".$qcolors);
		while(<QCOLS>) { 
			chomp;
			push(@backup_cols,$_);
		}
		close QCOLS;
	} else {
		for (my $i = 2; $i != scalar @KellyColorRGB; $i++) { 
			push(@backup_cols,$KellyColorRGB[$i]);
		}
	}
	#have names for each species: try and use them to infer good colors for Q matrix
	#which ring should I color
	if ($color_ref_ring == -1 && ! $not_concentric ) { $color_ref_ring = (scalar @qmatrix) - 1; }
	print STDERR "assigning label colors onto ring ".$color_ref_ring."\n";
	# how many colors do I need?
	my $width_here = $all_matrix_widths[$color_ref_ring];
	print STDERR "need to find ".$width_here." colors for this ring\n";
	#array of colors of width_here
	my @qcols;
	#score total of weights per color
	my %cols_weight_per_Q;
	
	my @taxlist = keys %{$all_qdata[$color_ref_ring]};
	my %new_species_name_colors;
	if ($fix_spp_names_Ldon) { 
		#horrible - names for Ldon Qmatrix don't match Tree
		foreach my $k (keys %species_name_colors) { 
			my @bits = split("_",$k);
			pop @bits;
			my $new_label = join("_",@bits);
			#$new_label =~ s/\s/_/;
			print STDERR "replacing label '".$k."' with '".$new_label."'\n";
			$new_species_name_colors{$new_label} = $species_name_colors{$k};
		}
	  #  ($new_label = $b) =~ s/\s/_/g;
	} else { 
		foreach my $k (keys %species_name_colors) { 
			$new_species_name_colors{$k} = $species_name_colors{$k};
		}
	}
	foreach my $b (@taxlist) { 
		my $color_here = "-1";
		
		if (exists $new_species_name_colors{$b} ) { 
			$color_here = $new_species_name_colors{$b};
		} else  {
			print STDERR "(learncolors) warning - no colors found for species '".$b."'\n";
		}
		if ($color_here ne "-1" ) { 
			foreach (my $i = 0; $i != $width_here; $i++ ) { 
				my $data_point_here = ${$all_qdata[$color_ref_ring]}{$b}[$i];
				$cols_weight_per_Q{$i}{$color_here}	+= $data_point_here;
			}
		}
	}
	my %already_used_colors;
	my $max_index_in_altcolors_used = -1;
	
	if ($flip_learn_colors) { 
		my %unique_colors;
		foreach (my $i = 0; $i != $width_here; $i++ ) { 
				foreach my $k (keys %{$cols_weight_per_Q{$i}} ) {
					$unique_colors{$k}++;
				}
		}
		print STDERR "have ".(scalar keys %unique_colors)." colors to assign to ".$width_here." groups\n";
		my %already_filled_Q;
		
		while (scalar keys %already_filled_Q != $width_here && scalar keys %already_used_colors != scalar keys %unique_colors) { 
			my $global_max_weight = -1;
			my $global_max_col = "none";
			my $global_max_q = -1;
			print STDERR "filled ".(scalar keys %already_filled_Q)." out of ".$width_here." Qs\n";
			print STDERR "used up ".(scalar keys %already_used_colors)." out of ".(scalar keys %unique_colors)." colors\n";
			foreach (my $i = 0; $i != $width_here; $i++ ) { 
				foreach my $k (keys %{$cols_weight_per_Q{$i}} ) {
					if ( ! exists $already_used_colors{$k} && ! exists $already_filled_Q{$i} ) { 
						if ($cols_weight_per_Q{$i}{$k} > $global_max_weight ) {
							$global_max_weight = $cols_weight_per_Q{$i}{$k};
							$global_max_col = $k;
							$global_max_q = $i;
						}
					}
				}		
			}	
			print STDERR "found max weight of ".$global_max_weight." for color ".$global_max_col." at q=".$global_max_q."\n";
			$already_used_colors{$global_max_col}++;
			$already_filled_Q{$global_max_q}++;
			$qcols[$global_max_q] = $global_max_col;
		}
		if (scalar keys %already_filled_Q < $width_here) { 
		#need to use some backup colors
			foreach (my $i = 0; $i != $width_here; $i++ ) { 
				if (! exists $already_filled_Q{$i} ) { 
					
					$qcols[$i] = $backup_cols[++$max_index_in_altcolors_used];
					print STDERR "using backup color (".$backup_cols[$max_index_in_altcolors_used].") for q=".$i."\n";
					$already_used_colors{$backup_cols[$max_index_in_altcolors_used]}++;
				}
			}
		}
		
	} else { 
		foreach (my $i = 0; $i != $width_here; $i++ ) { 
			if (exists $cols_weight_per_Q{$i} ) { 
				my $max_weight = -1;
				my $max_col = "none";
				print STDERR "color weight matrix for q=".$i." is:\n";
				foreach my $k (keys %{$cols_weight_per_Q{$i}} )  {
					print STDERR "\t".$k."\t".$cols_weight_per_Q{$i}{$k}."\n";
				}
				foreach my $k (keys %{$cols_weight_per_Q{$i}} ) { 
					if ($cols_weight_per_Q{$i}{$k} > $max_weight ) { 
						$max_weight = $cols_weight_per_Q{$i}{$k};
						$max_col = $k;
					}
				}
				if ($max_col ne "none" && !(exists $already_used_colors{$max_col}) ) { 
					$qcols[$i] = $max_col;
					print STDERR "assigning ".$max_col." for q=".$i."\n";
					$already_used_colors{$max_col}++;
				} else {
					#assign back_up color
					print STDERR "using backup color for q=".$i."\n";
					$qcols[$i] = $backup_cols[++$max_index_in_altcolors_used];
					$already_used_colors{$backup_cols[$max_index_in_altcolors_used]}++;
				}
			}
		}
	}
	$allcolors[$color_ref_ring] = \@qcols;
} 
	
if ( ! $nomatchcolors) {
	if ($color_ref_ring == -1 && ! $not_concentric) {
		#assign colors in concentric rings
		
		if (scalar @qmatrix > 1 ) { 
			print STDERR "trying to find matches between ".(scalar @qmatrix)." concentric rings\n";
			my @taxlist = keys %{$all_qdata[0]};
			#iterate through rings != $color_ref_ring;
			
			
			
			
			for(my $qmi = scalar @qmatrix - 2; $qmi >= 0; $qmi--) { 
			my %sum_across_all;
				my $my_ref_ring = $qmi + 1;
				my $matrix_width_ref = $all_matrix_widths[$my_ref_ring];
				#hash of hashes for this ring, index by this ring_index, ref_ring_index
		#		print STDERR "comparing ring ".$qmi." with ".$my_ref_ring."\n";
				my $matrix_width = $all_matrix_widths[$qmi];	
				foreach my $b (@taxlist) { 
					foreach (my $i = 0; $i != $matrix_width; $i++ ) { 
						my $data_point_here = ${$all_qdata[$qmi]}{$b}[$i];
						foreach (my $j = 0; $j != $matrix_width_ref; $j++ ) { 
							
							my $data_point_there = ${$all_qdata[$my_ref_ring]}{$b}[$j];
						#	if ($data_point_here > $min_drawable ) { 
							#	if ($i == 2) { print STDERR $b." ".$i." ".$j." ".$data_point_here." ".$data_point_there."\n"; }
								$sum_across_all{$i}{$j} += abs($data_point_here - $data_point_there);
						#	}
						}
					}
				}
				my @map_j;
				my %already_used_colors;
			#	print STDERR "matrix of matches:\n";
			#	foreach (my $i = 0; $i != $matrix_width; $i++) { 
			#		my @dists;
			#		foreach (my $j = 0; $j != $matrix_width_ref; $j++ ) { 
			#			push(@dists,$sum_across_all{$i}{$j});
			#		}
			#		print STDERR join("\t",@dists)."\n";
			#	}
			#	print STDERR "\n";
				foreach (my $i = 0; $i != $matrix_width; $i++ ) { 
					my $max_value = -1;
					my $max_j = -1;
					my $global_max = -1;
					my $global_max_j = -1;
					foreach (my $j = 0; $j != $matrix_width_ref; $j++ ) { 

						if ($sum_across_all{$i}{$j}  < $max_value || $max_value == -1) { 
							if (! exists $already_used_colors{$j}|| $multicolors) { 
								$max_value = $sum_across_all{$i}{$j} ;
								$max_j = $j;
							} else {
					#			print STDERR "i=".$i." color ".$j." already used.. skipping\n";
							}
							
						}
						if ($sum_across_all{$i}{$j}  < $global_max || $global_max == -1) { 
							$global_max = $sum_across_all{$i}{$j} ;
							$global_max_j = $j;
						}
						
						#print $i.",".$j."=".$sum_across_all{$i}{$j}."\n";
					}	
					if ($max_j == -1) { die("can't find a value!!\n"); }
				#	print STDERR "setting i=".$i." to ".$max_j."\n";
					if ($max_j != $global_max_j) { print STDERR "warning: assigning color ".$max_j." for ".$i." when should be ".$global_max_j."\n"; }
					push(@map_j,$max_j);
					$already_used_colors{$max_j}++;
					$allcolors[$qmi][$i] = $allcolors[$my_ref_ring][$max_j];
					
				}
	#			print STDERR " map is : ".join(" ",@map_j)."\n";
	#			print STDERR " ref ring is ".join(" ",@{$allcolors[$my_ref_ring]})."\n";
	#			print STDERR " this ring is ".join(" ",@{$allcolors[$qmi]})."\n";
		
			}
		}
	} else {
		#assign colors against reference ring 
		$color_ref_ring = (scalar @qmatrix) - 1; 
		if (scalar @qmatrix > 1 ) { 
			#trying to find matches between rings;
			print STDERR "trying to find matches between rings\n";
			#matching from outside in: will usually be biggest ring.
			my @taxlist = keys %{$all_qdata[0]};
			#iterate through rings != $color_ref_ring;
			my $matrix_width_ref = $all_matrix_widths[$color_ref_ring];
			for(my $qmi = 0; $qmi != scalar @qmatrix; $qmi++) { 
				my %sum_across_all;
			
				if ($qmi != $color_ref_ring) { 
					#hash of hashes for this ring, index by this ring_index, ref_ring_index
					print STDERR "comparing ring ".$qmi." with ".$color_ref_ring."\n";
					my $matrix_width = $all_matrix_widths[$qmi];	
					foreach my $b (@taxlist) { 
						foreach (my $i = 0; $i != $matrix_width; $i++ ) { 
							my $data_point_here = ${$all_qdata[$qmi]}{$b}[$i];
							foreach (my $j = 0; $j != $matrix_width_ref; $j++ ) { 
								my $data_point_there = ${$all_qdata[$color_ref_ring]}{$b}[$j];
								$sum_across_all{$i}{$j} += abs($data_point_here - $data_point_there);
							}
							
						}
					}
					my @map_j;
					my %already_used_colors;
					foreach (my $i = 0; $i != $matrix_width; $i++ ) { 
						my $max_value = -1;
						my $max_j = -1;
						my $global_max = -1;
						my $global_max_j = -1;
					
						foreach (my $j = 0; $j != $matrix_width_ref; $j++ ) { 
							if ($sum_across_all{$i}{$j}  < $max_value || $max_value == -1) { 
								if (! exists $already_used_colors{$j} || $multicolors) { 
									$max_value = $sum_across_all{$i}{$j} ;
									$max_j = $j;
								} else { 
									print STDERR "color ".$j." already used.. skipping\n";
								}
							}
							#print $i.",".$j."=".$sum_across_all{$i}{$j}."\n";
							if ($sum_across_all{$i}{$j}  < $global_max || $global_max == -1) { 
								$global_max = $sum_across_all{$i}{$j} ;
								$global_max_j = $j;
							}	
						
						}
						
						if ($max_j == -1) { die("can't find a value!!\n"); }
						if ($max_j != $global_max_j) { print STDERR "warning: assigning color ".$max_j." for ".$i." when should be ".$global_max_j."\n"; }
				
						push(@map_j,$max_j);
						$already_used_colors{$max_j}++;
						$allcolors[$qmi][$i] = $allcolors[$color_ref_ring][$max_j];
					}
				}
			}
		}
	}
 
}

#try and keep same order of colors for consecutive rings in @matrix
my $reorder_colors = 1;

for(my $qmi = 0; $qmi != scalar @qmatrix; $qmi++) { 
	#my $qmatrix = $qmatrix[$qmi];
	my $dist_per_unit = (( $space_for_bars /  scalar @qmatrix) - $bar_padding); 
	
	my %qdata = %{$all_qdata[$qmi]};
	my $matrix_width = $all_matrix_widths[$qmi];
	
	
	
	#allcolors is an array of arrays
	my @qcols = @{$allcolors[$qmi]};
	
	
	
# make list of leaf nodes in sensible order
	my @leaf_nodes_sorted_by_angle = sort {  $node_angles{$a} cmp $node_angles{$b} } @{$pt->{leaflist}};
	my @ordered_angles;
#fast look-up in that list
	my %node_to_sortedarray_pos;
	for ( my $i = 0; $i != scalar @leaf_nodes_sorted_by_angle; $i++) { 
		$ordered_angles[$i] = $node_angles{$leaf_nodes_sorted_by_angle[$i]};
	#	print STDERR $i." ".$leaf_nodes_sorted_by_angle[$i]->{label}." ".$node_angles{$leaf_nodes_sorted_by_angle[$i]}."\n";
		$node_to_sortedarray_pos{$leaf_nodes_sorted_by_angle[$i]->{label}} = $i; 
	}
	
#	my @stuff = keys %qdata;
#	my @test_set = (0,15,30,50);
#	my @test_qs;
#	foreach my $t (@test_set) { push(@test_qs,$stuff[$t]); }
#	my $b = $stuff[0];
#	foreach my $b (@test_qs){ 
	my $bar_height = 1.0 * $dist_per_unit;

	foreach my $b (keys %qdata){ 
		#print STDERR "data is ".join(" ",@{$qdata{$b}})."\n";
		
		
		my $labelhash = return_label_hash($pt);
		my $bar_angle = -1;
		
		#----- this rubbish is to deal with the fact that names in q matrix differ a bit from the tree ---
		#-- in a sane world, this wouldn't be needed. If $fix_spp_names_Ldon is 0, then this doesn't do anything --
		# labelhash with "Ldon style" (qmatrix) style names, not tree names
		my %new_labelhash;
		#translator - give it the qmatirx ($b) name, get the tree name.
		my %name_translator;
		
		if (! $fix_spp_names_Ldon) {
			#do nothing
			foreach my $names (keys %{$labelhash}) { 
				$new_labelhash{$names} = ${$labelhash}{$names};
				$name_translator{$names} = $names;
			}
	
		} else {
			#fix species names
			#
			foreach my $names (keys %{$labelhash}) { 
				my @bits = split(/[\s+_]/,$names);
				pop @bits;
				my $new_label = join(" ",@bits);
				$new_label =~ s/\s/_/;
			#	print STDERR " renaming '".$names."'-->'".$new_label."'\n";
				$new_labelhash{$new_label} = ${$labelhash}{$names};	
				$name_translator{$new_label} = $names;
			}
	#		print STDERR "'".$b."'\n";	
		}
		#-----
		
		$bar_angle = $node_angles{$new_labelhash{$b}};
		my $start_x = 0;
		if ($annularsectors) { 
			#draw things as proper slices of a circle, rather than as rectangles
			#
			
			
		#	print STDERR "sectors of total width ".$bar_height." centered at ".$bar_angle."\n";
			
			#angles shared for particular individual
			my $central_angle_i = $node_to_sortedarray_pos{$name_translator{$b}};
			my $central_angle = $ordered_angles[$central_angle_i];
		#	print STDERR "label is at i=".$central_angle_i."(out of ".(scalar @leaf_nodes_sorted_by_angle).") angle=".$central_angle."\n";
			my $angle1 = 999999;
			my $angle2 = 999999;
			if ($central_angle_i > 0) { 
				$angle1 = $ordered_angles[$central_angle_i -1];
			} else { $angle1 =  $ordered_angles[scalar @leaf_nodes_sorted_by_angle - 1]; $angle1 -= 2 * pi; }
			if ($central_angle_i < (scalar @leaf_nodes_sorted_by_angle - 1) ) {
				$angle2 = $ordered_angles[$central_angle_i + 1];
			} else { $angle2 = $ordered_angles[0]; $angle2 += 2 * pi}

		#	print STDERR "angle1 = ".$angle1." angle2 = ".$angle2."\n";
		
		
			my $start_angle = $central_angle - (($central_angle - $angle1) * ((1.0 - $annularsectorgap_p) / 2.0));
			my $end_angle = $central_angle + (($angle2 - $central_angle) * ((1.0 - $annularsectorgap_p) / 2.0));
			if ($start_angle > 2*pi ) { $start_angle -= 2 *pi; }
			if ($end_angle > 2*pi ) { $end_angle -= 2 *pi; }
		#	print STDERR "start_angle = ".$start_angle." end_angle = ".$end_angle."\n";
			#die();
			#this is code for EACH sector!
			print $svgfileglob  "<g>\n";
			#loop across Ks
			my @wedge_starts;
			my @wedge_widths;
			
		#	foreach (my $i = 0; $i != 1; $i++ ) { 
			my @sorted_i = (0..$matrix_width );
			
			if ($reorder_colors) { 
				my @indices = (0..$matrix_width-1);
				@sorted_i = sort  { $qcols[$::a] cmp $qcols[$::b] } @indices;
	#			@sorted_i = sort  { $::a <=> $::b } @indices;

			} 
			
				foreach my $i (@sorted_i) { 
			
			
					my $this_height = $bar_height * ${$qdata{$b}}[$i];
				
					#skip drawing arcs that are too small to see!
					if (${$qdata{$b}}[$i] > $min_drawable) { 
				
						my $barcolor = $qcols[$i];
						#not sure if i need this stuff..
						my $min = min($start_angle,$end_angle);
						my $max = max($start_angle,$end_angle);
		
						my $large_arc = $max - $min <= pi ? 0 : 1;
						my $sweep_swap = 0;
						#my $this_height = $bar_height * 1.0;
				
						#print STDERR "height is ".$this_height."\n";
		
						if ($start_angle > $end_angle) { 
							$large_arc = abs(1-$large_arc);
							$sweep_swap = 1;	
							print STDERR "SWEEP SWAPPED FOR ".$b."\n";	
							print STDERR "label is at i=".$central_angle_i."(out of ".(scalar @leaf_nodes_sorted_by_angle).") angle=".$central_angle."\n";
							print STDERR "angle1 = ".$angle1." angle2 = ".$angle2."\n";
							print STDERR "start_angle = ".$start_angle." end_angle = ".$end_angle."\n";
			
						}
				#			$min += $wedge_angle_fiddle;
				#			$max -= $wedge_angle_fiddle;
				#		}	else { 
				#			$min -= $wedge_angle_fiddle;
				#			$max += $wedge_angle_fiddle;
				#		}
			
			#			print STDERR "my outer radius is ".($wedge_starts[$i] +$wedge_widths[$i] )."\n";
						my $real_start_x = $diagram_base_radius + $space_for_wedges+ $space_for_labels + $start_x + $bar_padding + (($dist_per_unit + $bar_padding) * $qmi) ;
						my $arc1_startx = $real_start_x * cos($min);
						my $arc1_starty =  $real_start_x * sin($min);
						my $arc1_endx  =  $real_start_x * cos($max);
						my $arc1_endy =  $real_start_x * sin($max);
						my $arc2_startx = ($real_start_x + $this_height ) * cos($min);
						my $arc2_starty =  ($real_start_x + $this_height)* sin($min);
						my $arc2_endx  =  ($real_start_x +$this_height )* cos($max);
						my $arc2_endy =  ($real_start_x + $this_height) * sin($max);
						print $svgfileglob  "<path d=\"M".($diagram_base_radius + $space_for_wedges+ $space_for_labels + $space_for_bars + $arc1_startx).",".($diagram_base_radius + $space_for_labels +$space_for_wedges+ $space_for_bars +$arc1_starty);
						if ($sweep_swap) { 
							print $svgfileglob  " A".$real_start_x.",".$real_start_x." 0 ".$large_arc.",0 ";
						} else{ 
							print $svgfileglob  " A".$real_start_x.",".$real_start_x." 0 ".$large_arc.",1 ";
						}
						print $svgfileglob  " ".($diagram_base_radius + $space_for_wedges+ $space_for_labels + $space_for_bars + $arc1_endx).",".($diagram_base_radius + $space_for_wedges+ $space_for_labels + $space_for_bars + $arc1_endy);
						print $svgfileglob  " L".($diagram_base_radius + $space_for_wedges+ $space_for_labels + $space_for_bars + $arc2_endx)." ".($diagram_base_radius + $space_for_wedges+ $space_for_labels + $space_for_bars + $arc2_endy)." ";
						if ($sweep_swap) { 
							print $svgfileglob  " A".($real_start_x + $this_height) .",".($real_start_x  + $this_height) ." 0 ".$large_arc.",1 ";
						} else{ 
							print $svgfileglob  " A".($real_start_x + $this_height) .",".($real_start_x  + $this_height) ." 0 ".$large_arc.",0 ";
						}
						print $svgfileglob  " ".($diagram_base_radius + $space_for_wedges+ $space_for_labels + $space_for_bars + $arc2_startx).",".($diagram_base_radius + $space_for_wedges+ $space_for_labels + $space_for_bars + $arc2_starty);
						print $svgfileglob  " L".($diagram_base_radius + $space_for_wedges+ $space_for_labels + $space_for_bars + $arc1_startx)." ".($diagram_base_radius + $space_for_wedges+ $space_for_labels + $space_for_bars + $arc1_starty)." ";
						print $svgfileglob  "\" style=\"stroke:none; fill:".$barcolor.";\"/>\n"; #fill-opacity=\"0.5\"/>\n";
				
					}
					$start_x += $this_height;
				
			}
			print $svgfileglob  "</g>\n";
		} else { 
			#MUCH EASIER _ JUST DRAW BARS!!
			print STDERR "bar of ".$bar_height." at ".$bar_angle."\n";
			my $cosa = cos($bar_angle);
			my $sina = sin($bar_angle);
			#my $height = $genome_size{$spp} * $mm_per_bp;
	
			my $transx = ($diagram_base_radius + $space_for_wedges+ $space_for_labels + $bar_padding + (($dist_per_unit + $bar_padding) * $qmi) )  * $cosa;
			my $transy =  ($diagram_base_radius + $space_for_wedges+ $space_for_labels + $bar_padding + (($dist_per_unit + $bar_padding) * $qmi) ) * $sina;

			print $svgfileglob  "<g transform=\"matrix(".$cosa.",".$sina.",".(-1.0 * $sina).",".$cosa.",".($diagram_base_radius + $space_for_labels +$space_for_wedges+ $space_for_bars + $transx).",".($diagram_base_radius + $space_for_labels +$space_for_wedges+ $space_for_bars+ $transy).")\">\n";
	#	my $barcolor = "#000000";
	#	if ($colorbarbywedge)  {
	#		my $group = $groups[$colorbarbywedge-1]{$b};
	#		my $col = $group_colors[$colorbarbywedge-1]{$group};
		#	print STDERR "'".$n->{label}."' is in '".$group."', and so will be '".$col."'\n";
	#		$barcolor = $col;
	
			foreach (my $i = 0; $i != $matrix_width; $i++ ) { 
				my $this_height = $bar_height * ${$qdata{$b}}[$i];
				if (${$qdata{$b}}[$i] > $min_drawable) { 
					my $barcolor = $qcols[$i];
					print $svgfileglob   "<rect x=\"".$start_x."\" y=\"".(-0.5 * $bar_width)."\" height=\"".$bar_width."\" width=\"".$this_height."\"  style=\"fill:".$barcolor."\"/>\n";
				}
				$start_x += $this_height;
			}
			print $svgfileglob  "</g>\n";
		}
	}
}

if (! $skipSVGheader )  {
	print $svgfileglob  "</svg>\n";
}
if ($outfile ne "-") { 
	close $svgfileglob ;
}
	


#___________________________
#traversal functions
#---------------	

#this version re-sizes to the barwidth and height_per_unit
sub BarPlotSVG_flex {
	my $values_array = shift;
	my $labels_array = shift;
	my $colors_array = shift;
	
	my $pixel_per_unit = shift;
	my $pixel_per_bar = shift;
	
	my $text_vertical = 1;
	
	my $SVGstring = "";
	my $xmargin = 30;
	my $ymargin = 10;
	my $ticklen = 2;
	##proportion og barwidth used as a gap
	my $bar_gap = 0.2;
	
	my $maxval = max( @{$values_array} );
		my $order_of_mag = 10**floor(log($maxval) / log(10));
	my $close_supremum = ceil( $maxval / $order_of_mag) * $order_of_mag;	
	print STDERR "drawing axis up to ".$close_supremum." for ".$maxval."\n";
	my @axis_ticks = (0..ceil( $maxval / $order_of_mag));
	foreach(@axis_ticks){ $_ *= $order_of_mag; };
	print STDERR "ticks at ".join(" ",@axis_ticks)."\n";
	my $axis_max = $axis_ticks[scalar @axis_ticks -1];
	
	my $height = ($axis_max * $pixel_per_unit) + $ymargin;
	
	#my $pixel_per_unit = ($height - $ymargin) / $axis_max;
#	my $pixel_per_bar = ( ( $width-$ticklen)  - $xmargin) / ( scalar @{$values_array});
	
	my $width = ($pixel_per_bar * ( scalar @{$values_array})) + $xmargin + $ticklen;
	for (my $i = 0; $i != scalar @{$values_array}; $i++ ) { 
		my $col = ${$colors_array}[$i];
		my $value = ${$values_array}[$i];
		$SVGstring .= "<rect x=\"".(($i * $pixel_per_bar) + $xmargin + ($pixel_per_bar * $bar_gap)) ."\" y=\"".($height - ($value * $pixel_per_unit))."\" height=\"".($value * $pixel_per_unit)."\" width=\"".($pixel_per_bar * ( 1.0 - $bar_gap))."\" fill=\"".$col."\"/>\n";

		#IF YOU MESS WITH POSITIONING HERE, MAKE SURE YOU ALSO CHANGE THE CENTER OF ROTATION!s
		if ($text_vertical == 0) { 
		$SVGstring .= "<text transform=\"rotate(45 ".(($i * $pixel_per_bar) + $xmargin + ($pixel_per_bar * $bar_gap) + (0.3 * ($pixel_per_bar * ( 1.0 - $bar_gap))))." ".($height).")\" text-anchor=\"start\" alignment-baseline=\"hanging\" x=\"".(($i * $pixel_per_bar) + $xmargin + ($pixel_per_bar * $bar_gap) + (0.3 * ($pixel_per_bar * ( 1.0 - $bar_gap))))."\" y=\"".($height)."\" font-family=\"arial\"  font-size=\"".($pixel_per_bar/ 1.0)."\"  >".(${$labels_array}[$i])."</text>\n";
		} else {
		$SVGstring .= "<text transform=\"rotate(90 ".(($i * $pixel_per_bar) + $xmargin + ($pixel_per_bar * $bar_gap) + (0.3 * ($pixel_per_bar * ( 1.0 - $bar_gap))))." ".($height + 2).")\" text-anchor=\"start\" alignment-baseline=\"middle\" x=\"".(($i * $pixel_per_bar) + $xmargin + ($pixel_per_bar * $bar_gap) + (0.3 * ($pixel_per_bar * ( 1.0 - $bar_gap))))."\" y=\"".($height + 2)."\" font-family=\"arial\"  font-size=\"".($pixel_per_bar/ 0.7)."\"  >".(${$labels_array}[$i])."</text>\n";
		
		}
	}	
	#x axis
	$SVGstring .= "<line x1=\"".$xmargin."\" y1=\"".$height."\" x2=\"".$width."\" y2=\"".$height."\"  style=\"stroke:rgb(0,0,0);stroke-width:0.5\"  />\n";
	#y axis
	$SVGstring .= "<line x1=\"".$xmargin."\" y1=\"".(0)."\" x2=\"".$xmargin."\" y2=\"".$height."\"  style=\"stroke:rgb(0,0,0);stroke-width:0.5\"  />\n";
	#ticks
	foreach my $tick (@axis_ticks) { 
		$SVGstring .= "<line x1=\"".$xmargin."\" y1=\"".($height - ($tick * $pixel_per_unit))."\" x2=\"".($xmargin - $ticklen)."\" y2=\"".($height - ($tick * $pixel_per_unit))."\"  style=\"stroke:rgb(0,0,0);stroke-width:0.5\" />\n";
		$SVGstring .= "<text x=\"".($xmargin - $ticklen)."\" alignment-baseline=\"middle\" text-anchor=\"end\" y=\"".($height - ($tick * $pixel_per_unit))."\"  font-family=\"arial\"  font-size=\"".(($xmargin -$ticklen)/ 2.0)."\" >".$tick."</text>\n";
	}
	
	return ($SVGstring, $width, $height);
}

#this version is a fixed size
sub BarPlotSVG { 
	my $values_array = shift;
	my $labels_array = shift;
	my $colors_array = shift;
	my $width = shift;
	my $height = shift;
	
	my $SVGstring = "";
	my $xmargin = 30;
	my $ymargin = 10;
	my $ticklen = 2;
	##proportion og barwidth used as a gap
	my $bar_gap = 0.2;
	
	my $maxval = max( @{$values_array} );
		my $order_of_mag = 10**floor(log($maxval) / log(10));
	my $close_supremum = ceil( $maxval / $order_of_mag) * $order_of_mag;	
	print STDERR "drawing axis up to ".$close_supremum." for ".$maxval."\n";
	my @axis_ticks = (0..ceil( $maxval / $order_of_mag));
	foreach(@axis_ticks){ $_ *= $order_of_mag; };
	print STDERR "ticks at ".join(" ",@axis_ticks)."\n";
	my $axis_max = $axis_ticks[scalar @axis_ticks -1];
	my $pixel_per_unit = ($height - $ymargin) / $axis_max;
	my $pixel_per_bar = ( ( $width-$ticklen)  - $xmargin) / ( scalar @{$values_array});
	
	for (my $i = 0; $i != scalar @{$values_array}; $i++ ) { 
		my $col = ${$colors_array}[$i];
		my $value = ${$values_array}[$i];
		$SVGstring .= "<rect x=\"".(($i * $pixel_per_bar) + $xmargin + ($pixel_per_bar * $bar_gap)) ."\" y=\"".($height - ($value * $pixel_per_unit))."\" height=\"".($value * $pixel_per_unit)."\" width=\"".($pixel_per_bar * ( 1.0 - $bar_gap))."\" fill=\"".$col."\"/>\n";
		$SVGstring .= "<text transform=\"rotate(45 ".(($i * $pixel_per_bar) + $xmargin + ($pixel_per_bar * $bar_gap) + (0.3 * ($pixel_per_bar * ( 1.0 - $bar_gap))))." ".($height).")\" text-anchor=\"start\" alignment-baseline=\"hanging\" x=\"".(($i * $pixel_per_bar) + $xmargin + ($pixel_per_bar * $bar_gap) + (0.3 * ($pixel_per_bar * ( 1.0 - $bar_gap))))."\" y=\"".($height)."\" font-family=\"arial\"  font-size=\"".($pixel_per_bar/ 2.0)."\"  >".(${$labels_array}[$i])."</text>\n";
	}	
	#x axis
	$SVGstring .= "<line x1=\"".$xmargin."\" y1=\"".$height."\" x2=\"".$width."\" y2=\"".$height."\"  style=\"stroke:rgb(0,0,0);stroke-width:0.5\"  />\n";
	#y axis
	$SVGstring .= "<line x1=\"".$xmargin."\" y1=\"".(0)."\" x2=\"".$xmargin."\" y2=\"".$height."\"  style=\"stroke:rgb(0,0,0);stroke-width:0.5\"  />\n";
	#ticks
	foreach my $tick (@axis_ticks) { 
		$SVGstring .= "<line x1=\"".$xmargin."\" y1=\"".($height - ($tick * $pixel_per_unit))."\" x2=\"".($xmargin - $ticklen)."\" y2=\"".($height - ($tick * $pixel_per_unit))."\"  style=\"stroke:rgb(0,0,0);stroke-width:0.5\" />\n";
		$SVGstring .= "<text x=\"".($xmargin - $ticklen)."\" alignment-baseline=\"middle\" text-anchor=\"end\" y=\"".($height - ($tick * $pixel_per_unit))."\"  font-family=\"arial\"  font-size=\"".(($xmargin -$ticklen)/ 2.0)."\" >".$tick."</text>\n";
	}
	
	return $SVGstring;
}
#	
#		
#OLD STYPE WITH MRCA CALCS
#foreach my $k ( keys %groups_to_members) { 
#	my $color = $group_colors[$i]{$k};
#	print STDERR "set of ".(scalar @{$groups_to_members{$k}})."\n";
#	my $mrca = findMRCAofmanyleaves($groups_to_members{$k},$pt);
#	print STDERR "GOT MRCA=".$mrca."\n";
#	my ($min, $max) = returnMinMaxAnglesFromMRCA($mrca,\%node_angles);
#	print STDERR "min_angle = ".$min." max=".$max."\n";
#	my $large_arc = $max - $min <= pi ? 0 : 1;
#	my $arc1_startx = $wedge_starts[$i] * cos($min);
#	my $arc1_starty =  $wedge_starts[$i] * sin($min);
#	my $arc1_endx  =  $wedge_starts[$i] * cos($max);
#	my $arc1_endy =  $wedge_starts[$i] * sin($max);
#	my $arc2_startx = ($wedge_starts[$i] + $wedge_widths[$i] ) * cos($min);
#	my $arc2_starty =  ($wedge_starts[$i] + $wedge_widths[$i])* sin($min);
#	my $arc2_endx  =  ($wedge_starts[$i] +$wedge_widths[$i] )* cos($max);
#	my $arc2_endy =  ($wedge_starts[$i] + $wedge_widths[$i]) * sin($max);
#	print SVGF "<path d=\"M".($diagram_base_radius + $space_for_bars+ $space_for_labels + $arc1_startx).",".($diagram_base_radius + $space_for_bars+ $space_for_labels + $arc1_starty);
#	print SVGF " A".$wedge_starts[$i].",".$wedge_starts[$i]." 0 ".$large_arc.",1 ";
#	print SVGF " ".($diagram_base_radius+ $space_for_bars + $space_for_labels + $arc1_endx).",".($diagram_base_radius + $space_for_bars+ $space_for_labels + $arc1_endy);
#	print SVGF " L".($diagram_base_radius + $space_for_bars+ $space_for_labels + $arc2_endx)." ".($diagram_base_radius + $space_for_bars+ $space_for_labels + $arc2_endy)." ";
#	print SVGF " A".($wedge_starts[$i] + $wedge_widths[$i]) .",".($wedge_starts[$i] + $wedge_widths[$i]) ." 0 ".$large_arc.",0 ";
#	print SVGF " ".($diagram_base_radius+ $space_for_bars + $space_for_labels + $arc2_startx).",".($diagram_base_radius + $space_for_bars+ $space_for_labels + $arc2_starty);
#	print SVGF " L".($diagram_base_radius + $space_for_bars+ $space_for_labels + $arc1_startx)." ".($diagram_base_radius + $space_for_bars+ $space_for_labels + $arc1_starty)." ";
#	print SVGF "\" style=\"stroke: #000000; stroke-width:1; fill:".$color.";\" fill-opacity=\"0.5\"/>\n";
#}
#
#


sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}


sub return_label_hash { 
	my $tree = shift;
	my %label_hash;
#	print STDERR "leaflist is ".$tree->{leaflist}."\n";
	foreach my $n (@{$tree->{leaflist}}) { 
#		print STDERR $n."\n";
		my $lab = $n->{label};
		$label_hash{ $lab } = $n;
	}
	return \%label_hash;
}

sub returnMinMaxAnglesFromMRCA {
	my $mrca = shift;
	my $nodeanglesref = shift;
	my @leaves_below;
	if ($mrca->{leaf}) { 
		push(@leaves_below,$mrca); 
	} else {
		my @nodes_to_visit = ($mrca->{child});
		while (scalar @nodes_to_visit != 0 ) { 
			my $n = shift @nodes_to_visit;
			if ($n) { 
			if ($n->{leaf} ) { 
				push(@leaves_below,$n); 
				push(@nodes_to_visit,$n->{sibling}); 
			} else {
				push(@nodes_to_visit,$n->{child});
				push(@nodes_to_visit,$n->{sibling});
			}
			}
		}
	}
	#print STDERR "found ".(scalar @leaves_below)." leaves below ".$mrca."\n";
	my @angles;
	foreach my $l (@leaves_below) { 
		push(@angles,${$nodeanglesref}{$l});
	}
	return (min(@angles),max(@angles));
}

sub findMRCAofmanyleaves {
	my $leafnodes = shift;
	my $treeobj = shift;
	my %nodehitcounter;
	foreach my $n (@{$treeobj->{nodelist}}) {
		$nodehitcounter{$n} = 0;
	}
	foreach my $n (@{$leafnodes}) { 
	#	print STDERR "walking down from ".$n."\n";
		my $p = $n;
		while ($p) { 
	#		print STDERR "now at ".$p." incrementing from ".$nodehitcounter{$p}."\n";
			$nodehitcounter{$p}++;
			$p = $p->{anc};
		}
	}
	my $arbitrary_node = ${$leafnodes}[0];
	while ($nodehitcounter{$arbitrary_node} != scalar @{$leafnodes} && $arbitrary_node != $treeobj->{rootnode} ) { 
#		print STDERR $nodehitcounter{$arbitrary_node}." != ".(scalar @{$leafnodes})."\n";
		$arbitrary_node = $arbitrary_node->{anc}; 
	}
#'	print STDERR " think I have a hit at ".$arbitrary_node." unless it is ".$treeobj->{rootnode}."\n";
	return $arbitrary_node;
}


sub visit_postorder { 
	my $node = shift;
	my $node_angles_ref = shift;
	#print STDERR "".(scalar keys %{$node_angles_ref})." angles stored:\n";

	if ($node) { 
	#	print STDERR "starting on ".$node."\n";
		visit_postorder($node->{child},$node_angles_ref);
		visit_postorder($node->{sibling},$node_angles_ref);
		my @kids;
		my $p = $node->{child};
		while ($p) { 
			if (exists $$node_angles_ref{$p} ) { 
			
			
				push(@kids,${$node_angles_ref}{$p});
			} else { 
				die( "can't find angle for ".$p." child of ".$node."\n");
			}
			$p = $p->{sibling};
		}
		if (scalar keys @kids > 0) {
			${$node_angles_ref}{$node} = (sum @kids) / (scalar keys @kids);
	#		print STDERR join("<",@kids)."\n";;
		}
	#	print STDERR "finished with: ".$node." (leaf=".$node->{label}.")\n";

	}
}

sub node_placement_preorder { 
	my $node = shift;
	my $branch_start_radius = shift;
	my $branch_end_radius = shift;
	my $arc_start  = shift;
	my $arc_end = shift;
	my $rpb = shift;
	my $nodeangles = shift;
	if ($node) { 
		if($node->{anc}) { 
			${$branch_start_radius}{$node} = ${$branch_end_radius}{$node->{anc}}; 
		} else { 
			${$branch_start_radius}{$node} = 0;
		}
	#	print STDERR $node."->".$node->{edgelength}."\n";
		if (exists $node->{edgelength} ) { 
			${$branch_end_radius}{$node} = ${$branch_start_radius}{$node} + ( $node->{edgelength} * $rpb );
		} else { 
			${$branch_end_radius}{$node} = ${$branch_start_radius}{$node};
		}
		my @descendent_angles;
		my $p = $node->{child};
		while($p) { 
			push(@descendent_angles,${$nodeangles}{$p});
			$p = $p->{sibling};
		}
		${$arc_start}{$node} = min(@descendent_angles);
		${$arc_end}{$node} = max(@descendent_angles);
		node_placement_preorder($node->{sibling},$branch_start_radius,$branch_end_radius,$arc_start,$arc_end,$rpb,$nodeangles );
		node_placement_preorder($node->{child},$branch_start_radius,$branch_end_radius,$arc_start,$arc_end,$rpb,$nodeangles );
	}
	
	
}
