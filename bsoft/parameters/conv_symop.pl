#!/usr/bin/perl -w

# Script to convert the file symop.lib to STAR format

$thetime = localtime;

print "# symop.star: Converted from symop.lib on ", $thetime, "\n\n" ;
print "# Each set of symmetry operators is written in its own data block\n";
print "# The data block tag is the space group number\n";
print "# The original spacegroup definition line is given as a comment\n";
print "\n\n";

$count = 0;
while ( <> ) {
	@list = split;
#	If the last element index in the list is greater than 5 and
#	if the line does not contain a '*', it is a spacegroup definition line
	if ( $#list > 5 & !(/\*/) ) { 
		print "\n\n############################# ", $list[0], " #############################\n";
		print "\ndata_", $list[0], "\n\n";
		print "# ", $_, "\n";
		print "_symmetry.Int_Tables_number     ", $list[0], "\n";
		print "_symmetry.space_group_name_H-M  ", $list[3], "\n";
		print "_symmetry.cell_setting          ", $list[5], "\n";
		print "\nloop_\n";
		print "_symmetry_equiv.id\n";
		print "_symmetry_equiv.pos_as_xyz\n";
		$count = 0;
	}
	if ( /[XYZ]/ ) {
		@xyz = split '\*';
		$subcount = 0;
		while ( $subcount < $#xyz ) {
			$count++;
			print $count, "  ", $xyz[$subcount], "\n";
			$subcount++;
		}
		$count++;
		print $count, "  ", $xyz[$subcount];
	}
}
