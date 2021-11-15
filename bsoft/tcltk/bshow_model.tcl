##
# @file		bshow_model.tcl
#
# @brief	Procedures to create and destroy model components and links in images.
#
# @author	Bernard Heymann
# @date		Created: 20010516
# @date		Modified: 20210223

set modelfile ""
set modeldir ""
set comp_radius 0
set link_radius 0
set comp_selected 0
set link_selected []
set show_comp 1
set show_link 1
set show_comp_labels 1
set comp_selection "Current"
set selected_color "#ffffff"
set color_selection "Selection"
set types_selected []
set line_width 1

## @brief Menu for opening and closing parameter files
#

proc setupModelMenu {} {
	global modelfile
	
	.menuBar add cascade -menu .menuBar.model -label "Model" -underline 0
	menu .menuBar.model -tearoff 0
	.menuBar.model add command -label "Read model" -underline 0 \
		-command { readModel }
	.menuBar.model add command -label "Write model" -underline 0 \
		-command { writeModel }
	.menuBar.model add command -label "Edit model" -underline 0 \
		-command { editTextFile $modelfile }
	.menuBar.model add command -label "Set views" -underline 0 \
		-command { setViews }
	.menuBar.model add command -label "Delete all components" -underline 0 \
		-command { objectsDeleteAll }
	.menuBar.model add command -label "Delete non-selected components" -underline 0 \
		-command { deleteNonSelectedObjects }
	.menuBar.model add command -label "Add or update component type" -underline 0 \
		-command { createUpdateComponentType }
	.menuBar.model add command -label "Extract marked segments" -underline 0 \
		-command { extractSegments }
	.menuBar.model add command -label "Generate particles" -underline 0 \
		-command { componentsToParticles }
	.menuBar.model add command -label "Analyze component symmetry" -underline 0 \
		-command { componentSymmetry }
#	.menuBar.model add command -label "Help" -underline 0 \
#		-command { showURL "file://$Bsoft/doc/bshow/bshow_model.html" }
}

## @brief Dialog box to manipulate componets and links on the main canvas
#

proc Model { } {
	global theimg filename project_item
	global helv12 cour12
	global comp_radius link_radius
	global selected_color color_selection line_width comp_selection

	if ![Bimage exists $theimg] {
		tk_dialog .dialog "No image in memory!" \
				{Please read in an appropriate image first.} \
				info 0 {OK}
		return
	}
	
	if ![Bmodel exists] {
		set project_item [Bmodel create]
	}
	set c [getImageCanvas $theimg]
	set modid [Bmodel get $project_item id]
	set modsel [Bmodel get $project_item selection]
	set ncomp [Bmodel component $project_item count]
	set nlink [Bmodel link $project_item count]
	set ps [.pixel_size.x get]
	if { $comp_radius < $ps } { set comp_radius [expr 2*$ps] }
	if { $link_radius < $ps } { set link_radius [expr $ps] }
	set w .wmod
	if ![winfo exists $w] {
		toplevel $w
		
		frame $w.modid
		label $w.modid.tag -text "Model" -width 15 -anchor e
		label $w.modid.str -width 50 -text $modid -relief sunken -bd 1 \
				-font $helv12 -anchor w
		label $w.modid.sel -width 5 -text $modsel -relief sunken -bd 1 \
				-font $helv12 -anchor w

		frame $w.comp
		label $w.comp.tag -text "Components" -width 15 -anchor e
		label $w.comp.num -width 10 -text $ncomp -relief sunken -bd 1 \
				-font $helv12 -anchor w
		label $w.comp.tag2 -text "Selected" -width 15 -anchor e
		label $w.comp.sel -width 10 -text $ncomp -relief sunken -bd 1 \
				-font $helv12 -anchor w

		frame $w.link
		label $w.link.tag -text "Links" -width 15 -anchor e
		label $w.link.num -width 10 -text $nlink -relief sunken -bd 1 \
				-font $helv12 -anchor w
		label $w.link.tag2 -text "Selected" -width 15 -anchor e
		label $w.link.sel -width 10 -text $nlink -relief sunken -bd 1 \
				-font $helv12 -anchor w

		frame $w.curr
		label $w.curr.tag -text "Current" -width 15 -anchor e
		label $w.curr.id -width 10 -text " " -relief sunken -bd 1 \
				-font $helv12 -anchor w
		label $w.curr.tag2 -text "FOM" -width 15 -anchor e
		label $w.curr.fom -width 10 -text " " -relief sunken -bd 1 \
				-font $helv12 -anchor w

		frame $w.dist
		label $w.dist.tag -text "Distance" -width 15 -anchor e
		label $w.dist.val -width 10 -text " " -relief sunken -bd 1 \
				-font $helv12 -anchor w

		frame $w.switches
		label $w.switches.tag -text "Show" -width 10 -anchor e
		checkbutton $w.switches.components -text "components" -variable show_comp \
				-command { objectsDrawAll }
		checkbutton $w.switches.links -text "links" -variable show_link \
				-command { objectsDrawAll }
		checkbutton $w.switches.labels -text "labels" -variable show_comp_labels \
				-command { objectsDrawAll }
		label $w.switches.widthlabel -text "Line width" -width 10 -anchor e
		tk_optionMenu $w.switches.linewidth line_width 1 2 3 4 5 6 7 8 9 10

		frame $w.select
		label $w.select.tag -text "Select" -anchor w -width 20
		tk_optionMenu $w.select.menu comp_selection \
			"Current" "Selected" "Selected Types" "None" "All"

		frame $w.compradius
		button $w.compradius.tag -text "Component radius" -width 20 \
			-anchor e -command setComponentRadius
		entry $w.compradius.entry -width 10 -validate key -vcmd { string is double %P }
		$w.compradius.entry insert 0 $comp_radius

		frame $w.linkradius
		button $w.linkradius.tag -text "Link radius" -width 20 \
				-anchor e -command setLinkRadius
		entry $w.linkradius.entry -width 10 -validate key -vcmd { string is double %P }
		$w.linkradius.entry insert 0 $link_radius

		frame $w.color
		button $w.color.button -text "Color" -anchor w -width 20 -command setColor
		canvas $w.color.well -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $selected_color
		tk_optionMenu $w.color.menu color_selection \
			"All" "Selection" "By Density" "By FOM"

		frame $w.fom
		label $w.fom.label -text "FOM cutoff" -anchor w
		scale $w.fom.scale -orient horizontal -length 200 -from 0 -to 1 \
			-command { setFOMcutoff } -resolution 0.001

		frame $w.neighbor
		button $w.neighbor.prev -text "Previous model" \
				-relief raised -command "getNeighborImageOfType2 -1"
		button $w.neighbor.next -text "Next model" \
				-relief raised -command "getNeighborImageOfType2 1"

		label $w.header -font $cour12 -anchor w \
			-text " Type       Number   Mass    FOM    Select Image File"

		frame $w.frame
		text $w.table -relief sunken -bd 2 -width 70 \
			-yscrollcommand "$w.yscroll set" \
			-setgrid 1 -highlightthickness 0 -pady 2 -padx 3 -font $cour12
		scrollbar $w.yscroll -command "$w.table yview" \
			-highlightthickness 0 -orient vertical
		grid $w.table -in $w.frame -padx 1 -pady 1 \
			-row 0 -column 0 -rowspan 1 -columnspan 1 -sticky news
		grid $w.yscroll -in $w.frame -padx 1 -pady 1 \
			-row 0 -column 1 -rowspan 1 -columnspan 1 -sticky news
		grid rowconfig    $w.frame 0 -weight 1 -minsize 0
		grid columnconfig $w.frame 0 -weight 1 -minsize 0

		
		pack $w.modid.tag $w.modid.str $w.modid.sel\
				-side left -pady 5 -padx 10
		pack $w.comp.tag $w.comp.num $w.comp.tag2 $w.comp.sel \
				-side left -pady 5 -padx 10
		pack $w.link.tag $w.link.num $w.link.tag2 $w.link.sel \
				-side left -pady 5 -padx 10
		pack $w.curr.tag $w.curr.id $w.curr.tag2 $w.curr.fom -side left -pady 5 -padx 10
		pack $w.dist.tag $w.dist.val -side left -pady 5 -padx 10
		pack $w.switches.tag $w.switches.components \
			$w.switches.links $w.switches.labels \
			$w.switches.widthlabel $w.switches.linewidth -side left -pady 5 -padx 5
		pack $w.select.tag $w.select.menu -side left -pady 5 -padx 10
		pack $w.compradius.tag $w.compradius.entry -side left -pady 5 -padx 10
		pack $w.linkradius.tag $w.linkradius.entry -side left -pady 5 -padx 10
		pack $w.color.button $w.color.well $w.color.menu -side left -pady 5 -padx 10
		pack $w.fom.label -side top -fill x
		pack $w.fom.scale -side top -fill x
		pack $w.neighbor.prev $w.neighbor.next \
			-side left -pady 5 -padx 5 -ipadx 5 -ipady 5
		pack $w.modid $w.comp $w.link $w.curr $w.dist $w.switches $w.select \
				$w.compradius $w.linkradius $w.color $w.fom $w.neighbor \
				-side top -fill x -pady 2
		pack $w.header -side top -fill x -padx 5
		pack $w.frame -expand yes -fill both -padx 5 -pady 5
		
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Model"
    wm iconname $w "Model"
	.menuBar.window entryconfigure "Model" -state normal
	bind $c <1> "componentCreateSelect %x %y"
	bind $c <2> "linkCreate %x %y"
	bind $c <Shift-1> "componentDelete %x %y"
	bind $c <Shift-2> "linkDelete %x %y"
	bind $c <B1-Motion> "componentMove %x %y"
	bind $c <B2-Motion> "componentMove %x %y"
	bind . <u> "componentMoveSelectedUpDown 1"
	bind . <d> "componentMoveSelectedUpDown -1"
	bind $w <Return> "objectsDrawAll"
	bind $w <Control-w> { destroy .wmod }
	bind $w <Command-w> { destroy .wmod }
	bind $w.color.well <1> { selectColor .wmod.color.well }
	trace add variable comp_selection write setSelection
	trace add variable line_width write objectsDrawAll
	objectsDrawAll
	updateCompTypeTable
}

proc toggleCompTypeSelection { id sel } {
	global project_item
#	puts "$id $sel"
	if { $sel } {
		Bmodel comptype $project_item deselect $id
	} else {
		Bmodel comptype $project_item select $id
	}
	updateCompTypeTable
}

## @brief Updates the component type table from the model parameters in memory.
#

proc updateCompTypeTable { } {
	global project_item types_selected
	set w .wmod
	if ![winfo exists $w] { return }
	if ![Bmodel exists] { return }
	set types [Bmodel comptype $project_item list]
#	puts "Updating type table"
	$w.table delete 0.0 end
	set types_selected []
	foreach {id fmod img_num num mass fom sel} $types {
		$w.table insert end \
			[format "%-10s %5d %9.0f %7.3f %5d %5d %s\n" $id $num $mass $fom $sel $img_num $fmod] $id
		if { $sel } {
			$w.table tag configure $id -background #aaaaff -borderwidth 1
			lappend types_selected $id
		} else {
			$w.table tag configure $id -background {}
		}
		$w.table tag bind $id <Button-1> "toggleCompTypeSelection $id $sel"
	}
#	puts $types_selected
}

## @brief Converts image coordinates to model coordinates.
#
# @param	x
# @param	y
# @param	z			Position in image coordinates.

proc convertImageToModel { x y z } {
	global theimg
	set ps [.pixel_size.x get]
	set ori [Bimage get $theimg origin [.image.scale get]]
	set x [expr ($x - [lindex $ori 0])*$ps]
	set y [expr ($y - [lindex $ori 1])*$ps]
	set z [expr ($z - [lindex $ori 2])*$ps]
#	puts "1 $x $y $z"
	return "$x $y $z"
}

## @brief Selects a component.
#
# @param	x
# @param	y			Position in canvas coordinates.

proc componentSelect { x y } {
	set w .wmod
	if ![winfo exists $w] { return }
	global theimg project_item
	global comp_selected
#	global tool
	global px py
	set c [getImageCanvas $theimg]
	set scale [.scale.scale get]
	set n [.image.scale get]
	set z [.slice.scale get]
	set px [expr floor( [$c canvasx $x] / $scale ) ]
	set py [expr floor( [$c canvasy $y] / $scale ) ]
	set pyi [yFlip $py]
	set id1 $comp_selected
	set loc [convertImageToModel $px $pyi $z]
	set id [Bmodel component $project_item select [lindex $loc 0] [lindex $loc 1] [lindex $loc 2]]
#	puts "selected comp = $id"
	set comp_selected $id
	if { $id > 0 } {
		componentDraw $id
	}
	set fom [Bmodel component $project_item fom $id]
	set dist [Bmodel component $project_item distance $id $id1]
	$w.comp.sel config -text [Bmodel component $project_item countselected]
	$w.curr.id config -text $id
	$w.curr.fom config -text [format "%8.2f" $fom]
	$w.compradius.entry delete 0 end
	$w.compradius.entry insert 0 [Bmodel component $project_item radius $id]
	if { $dist } {
		$w.dist.val config -text [format "%8.2f" $dist]
	} else {
		$w.dist.val config -text " "
	}
	return $id
}

## @brief Selects a link.
#
# @param	x
# @param	y			Position in canvas coordinates.

proc linkSelect { x y } {
	set w .wmod
	if ![winfo exists $w] { return }
	global theimg project_item
	global link_selected link_radius
#	global tool
	global px py
	set c [getImageCanvas $theimg]
	set scale [.scale.scale get]
	set n [.image.scale get]
	set z [.slice.scale get]
	set px [expr floor( [$c canvasx $x] / $scale ) ]
	set py [expr floor( [$c canvasy $y] / $scale ) ]
	set pyi [yFlip $py]
	set loc [convertImageToModel $px $pyi $z]
	set id [Bmodel link $project_item select [lindex $loc 0] [lindex $loc 1] [lindex $loc 2] $link_radius]
#	puts "selected link = $id ($loc)"
	if { [llength $id] > 1 } {
		linkDraw [lindex $id 0] [lindex $id 1]
	}
	set link_selected $id
	return $id
}

## @brief Creates or updates a component type.
#

proc createUpdateComponentType { } {
	if ![winfo exists .wmod] { return }
	global project_item
	set w .wnct
	catch {destroy $w}
	toplevel $w
	wm title $w "Component type"
	wm iconname $w "Type"
	tkwait visibility $w
	grab $w
	set type [Bmodel comptype $project_item firstselected]
	label $w.label -text "Create or update a component type" -anchor w
	frame $w.type
	label $w.type.label -text "Type name" -anchor w
	entry $w.type.entry -width 10
	$w.type.entry insert 0 [lindex $type 0]
	frame $w.mass
	label $w.mass.label -text "Mass" -anchor w
	entry $w.mass.entry -width 10 -validate key -vcmd { string is double %P }
	$w.mass.entry insert 0 [lindex $type 4]
	frame $w.file
	label $w.file.label -text "File name" -anchor w
	entry $w.file.entry -width 10
	$w.file.entry insert 0 [lindex $type 1]
	frame $w.img_num
	label $w.img_num.label -text "Image number" -anchor w
	entry $w.img_num.entry -width 10 -validate key -vcmd { string is integer %P }
	$w.img_num.entry insert 0 [lindex $type 2]
	frame $w.buttons
	button $w.buttons.set -text Set -command {
		set n [.image.scale get]
		Bmodel comptype $project_item create [.wnct.type.entry get] [.wnct.mass.entry get] \
			 [.wnct.file.entry get] [.wnct.img_num.entry get]
		destroy .wnct
	}
	button $w.buttons.cancel -text Cancel -command "destroy $w"
	pack $w.type.label $w.type.entry -side left -fill x
	pack $w.mass.label $w.mass.entry -side left -fill x
	pack $w.file.label $w.file.entry -side left -fill x
	pack $w.img_num.label $w.img_num.entry -side left -fill x
	pack $w.buttons.set $w.buttons.cancel -side left -expand 1
	pack $w.label $w.type $w.mass $w.file $w.img_num $w.buttons -side top -fill x -pady 2m
	tkwait window $w
	updateCompTypeTable
}

## @brief Creates or selects a component.
#
# @param	x
# @param	y			Component position in canvas coordinates.

proc componentCreateSelect { x y } {
	set w .wmod
	if ![winfo exists $w] { return }
	global theimg project_item
	global comp_selected
	global comp_radius selected_color
	global tool
	global px py
	set c [getImageCanvas $theimg]
	set comp_radius [$w.compradius.entry get]
	set scale [.scale.scale get]
	set n [.image.scale get]
	set z [.slice.scale get]
	set px [expr floor( [$c canvasx $x] / $scale ) ]
	set py [expr floor( [$c canvasy $y] / $scale ) ]
#	set h [Bimage get theimg height]
#	set pyi [expr $h - $py - 1]
	set pyi [yFlip $py]
	set loc [convertImageToModel $px $pyi $z]
	set id [Bmodel component $project_item select [lindex $loc 0] [lindex $loc 1] [lindex $loc 2]]
#	puts "selected ID = $id"
	if { $id < 1 && $tool == "model" } {
		if { [Bmodel comptype $project_item countselected] < 1 } { createUpdateComponentType }
		set id [Bmodel component $project_item create [lindex $loc 0] [lindex $loc 1] [lindex $loc 2]]
		Bmodel component $project_item radius $id $comp_radius
		Bmodel component $project_item color $selected_color $id
#		puts "created ID = $id  color = $selected_color"
	}
	if { $id > 0 } {
		componentDraw $id
	}
	set id1 $comp_selected
	set comp_selected $id
	set fom [Bmodel component $project_item fom $id]
	set dist [Bmodel component $project_item distance $id $id1]
	$w.comp.num config -text [Bmodel component $project_item count]
	$w.comp.sel config -text [Bmodel component $project_item countselected]
	$w.curr.id config -text $id
	$w.curr.fom config -text [format "%8.2f" $fom]
	$w.compradius.entry delete 0 end
	$w.compradius.entry insert 0 [Bmodel component $project_item radius $id]
	if { $dist } {
		$w.dist.val config -text [format "%8.2f" $dist]
	} else {
		$w.dist.val config -text " "
	}
	updateCompTypeTable
}

## @brief Creates a link if a starting component is selected
#
# @param	x
# @param	y			Ending component position in canvas coordinates.

proc linkCreate { x y } {
	set w .wmod
	if ![winfo exists $w] { return }
	global project_item comp_selected
	if { $comp_selected < 1 } { return }
	set n [.image.scale get]
	set id $comp_selected
	componentCreateSelect $x $y
	if { $id == $comp_selected } { return }
	set link_selected "$id $comp_selected"
	Bmodel link $project_item create $id $comp_selected
	linkDraw $id $comp_selected
#	set i [Bmodel link $project_item count $n]
#	$w.link.num config -text "$i"
	$w.link.num config -text [Bmodel link $project_item count]
	$w.link.num config -text [Bmodel link $project_item countselected]
}

## @brief Moves the selected component to the new coordinates.
#
# @param	x
# @param	y			Component position in canvas coordinates.

proc componentMove { x y } {
	set w .wmod
	if ![winfo exists $w] { return }
	global theimg project_item
	global comp_selected
	global px py
	if { $comp_selected < 1 } { return }
	set c [getImageCanvas $theimg]
	set scale [.scale.scale get]
	set n [.image.scale get]
	set ix [expr [$c canvasx $x] / $scale ]
	set iy [expr [$c canvasy $y] / $scale ]
	set dx [expr $ix - $px]
	set dy [expr $py - $iy]
	set px $ix
	set py $iy
#	puts "Moving $dx $dy"
	set ps [.pixel_size.x get]
	set dx [expr $dx * $ps]
	set dy [expr $dy * $ps]
	Bmodel component $project_item move $comp_selected $dx $dy
	objectsDrawAll
}

## @brief Moves the selected component one slice up or down.
#
# @param	dz			Number of slices (+ up, - down)

proc componentMoveSelectedUpDown { dz } {
	set w .wmod
	if ![winfo exists $w] { return }
	global theimg project_item comp_selected
#	puts "Moving  $comp_selected up/down $dz"
	if { $comp_selected < 1 } { return }
	set n [.image.scale get]
	set nslices [Bimage get $theimg nslices]
	set z [expr [.slice.scale get] + $dz]
	if { $z < 0 } { set z 0 }
	if { $z >= $nslices } { set z [expr $nslices - 1] }
	set ps [.pixel_size.x get]
	set dz [expr $dz * $ps]
	Bmodel component $project_item move $comp_selected 0 0 $dz
	.slice.scale set $z
}

## @brief Draws a component on the canvas.
#
# @param	id			Component identifier.

proc componentDraw { id } {
	global theimg project_item
	global show_comp_labels comp_selected line_width
	set c [getImageCanvas $theimg]
	set n [.image.scale get]
	set comp [Bmodel component $project_item img_coords $id]
	puts $comp
	set lbl ""
	set fill 0
	if { $show_comp_labels } { set lbl $id }
	drawSphere $c [lindex $comp 1] [lindex $comp 2] [lindex $comp 3] [lindex $comp 4] [lindex $comp 5] $line_width $lbl $fill component
}

proc drawSphere { c x y z rad col width lbl fill tag } {
#	puts "drawSphere $x $y $z"
	set scale [.scale.scale get]
	if { $rad < 1 } { set rad 1 }
	set zr [expr [.slice.scale get] - $z]
	if { $zr >= $rad } { return }
	set r [expr $rad*$rad - $zr*$zr]
	if { $r < 0 } { return }
	set xy_radius [expr $scale*sqrt($r)]
	set rad [expr $rad * $scale]
	set h [Bimage get theimg height]
	set x [expr $x * $scale ]
	set y [expr ($h - $y - 1) * $scale ]
#	set y [yFlip $py]
	set xmin [expr $x - $xy_radius ]
	set xmax [expr $x + $xy_radius + $scale - 1 ]
	set ymin [expr $y - $xy_radius ]
	set ymax [expr $y + $xy_radius + $scale - 1 ]
	set symbol "-outline"
	if { $fill } { set symbol "-fill" }
	$c create oval $xmin $ymin $xmax $ymax $symbol $col -width $width -tags $tag
	if { [string length $lbl] } {
		$c create text $xmax $ymax -text $lbl -fill $col -anchor w -tags $tag
	}
}

proc drawLine { c x1 y1 z1 x2 y2 z2 width col tag } {
#	puts "drawLine $x1 $y1 $z1 $x2 $y2 $z2"
	set scale [.scale.scale get]
	set width [expr int($width * $scale)]
	if { $width < 1 } { set width 1 }
	set z [.slice.scale get]
	set h [Bimage get theimg height]
	set xs [expr $x1 * $scale ]
	set ys [expr ($h - $y1 - 1) * $scale ]
	set zs $z1
	set xe [expr $x2 * $scale ]
	set ye [expr ($h - $y2 - 1) * $scale ]
	set ze $z2
	set zf -1.0
	if { abs($ze - $zs) > 0.001 } {
		set zf [expr ($z - $zs) / (1.0*($ze - $zs))]
	} else {
		if { abs($z - $zs) < 0.9 } { set zf 0.5 }
	}
	if { $xs == $xe && $ys == $ye } { set xe [expr $xe + 1] }
	if { $zf >= 0.0 && $zf <= 1.0 } {
		$c create line $xs $ys $xe $ye -fill $col -smooth yes -width $width -tags link
	}
}

## @brief Draws a link line connecting two components on the canvas.
#
# @param	id1			Component 1 identifier.
# @param	id2			Component 2 identifier.

proc linkDraw { id1 id2 } {
	global theimg project_item
	global show_comp fom_cutoff
	if { $show_comp < 1 } { return }
	set c [getImageCanvas $theimg]
	set scale [.scale.scale get]
	set n [.image.scale get]
	set linkcoor [Bmodel link $project_item img_coords $id1 $id2]
	set fom1 [Bmodel component $project_item fom $id1]
	set fom2 [Bmodel component $project_item fom $id2]
#	puts "Link:   $loc1 $fom1   $loc2 $fom2"
	if { $fom1 < $fom_cutoff || $fom2 < $fom_cutoff } { return }
	drawLine $c [lindex $linkcoor 0] [lindex $linkcoor 1] [lindex $linkcoor 2] \
		[lindex $linkcoor 3] [lindex $linkcoor 4] [lindex $linkcoor 5] \
		[lindex $linkcoor 6] [lindex $linkcoor 7] link
}

## @brief Draws all components and links on the canvas.
#

proc objectsDrawAll args {
#	puts "Drawing all objects"
	componentsDrawAll
	linksDrawAll
}

## @brief Draws all components on the canvas.
#

proc componentsDrawAll { } {
	global theimg project_item
	global show_comp fom_cutoff
	global show_comp_labels comp_selected line_width
	set c [getImageCanvas $theimg]
	$c delete component
	set w .wmod
	if ![winfo exists $w] { return }
	if { $show_comp < 1 } { return }
#	set n [.image.scale get]
#	puts "Drawing components ($project_item)"
	set comps [Bmodel component $project_item img_coords $fom_cutoff]
	foreach {id x y z rad color} $comps {
#		puts "Drawing $id"
		set lbl ""
		if { $show_comp_labels } { set lbl $id }
		set fill 0
#		puts $color
		drawSphere $c $x $y $z $rad $color $line_width $lbl $fill component
	}
	$w.modid.str config -text [Bmodel get $project_item id]
	$w.modid.sel config -text [Bmodel get $project_item selection]
	$w.comp.num config -text [Bmodel component $project_item count]
	$w.comp.sel config -text [Bmodel component $project_item countselected]
}

## @brief Draws all links on the canvas.
#

proc linksDrawAll { } {
	global theimg project_item
	global show_comp show_link fom_cutoff
#	puts "linksDrawAll"
	set c [getImageCanvas $theimg]
	$c delete link
	set w .wmod
	if ![winfo exists $w] { return }
	if { $show_comp < 1 } { return }
	if { $show_link < 1 } { return }
	set n [.image.scale get]
	set links [Bmodel link $project_item img_coords $fom_cutoff]
#	puts [llength $links]
	foreach {x1 y1 z1 x2 y2 z2 rad color} $links {
		drawLine $c $x1 $y1 $z1 $x2 $y2 $z2 $rad $color link
	}
	$w.link.num config -text [Bmodel link $project_item count]
	$w.link.sel config -text [Bmodel link $project_item countselected]
}

## @brief Deletes the selected component.
#
# @param	x
# @param	y			Cooordinates on the canvas.

proc componentDelete { x y } {
	if ![winfo exists .wmod] { return }
	global project_item
	set n [.image.scale get]
	set id [componentSelect $x $y]
#	puts "ID to be deleted = $id"
	if { $id > 0 } {
		Bmodel component $project_item delete $id
	}
	objectsDrawAll
}

## @brief Deletes the selected link.
#
# @param	x
# @param	y			Cooordinates on the canvas.

proc linkDelete { x y } {
	if ![winfo exists .wmod] { return }
	global project_item
	set n [.image.scale get]
	set id [linkSelect $x $y]
#	puts "link ID to be deleted = $id"
	if { [llength $id] > 1 } {
		Bmodel link $project_item delete [lindex $id 0] [lindex $id 1]
	}
	objectsDrawAll
}

## @brief Deletes all components and links.
#

proc objectsDeleteAll { } {
	global theimg
	set c [getImageCanvas $theimg]
	$c delete component link
	Bmodel delete
	set w .wmod
	if [winfo exists $w] {
		$w.comp.num config -text "0"
		$w.comp.sel config -text "0"
		$w.link.num config -text "0"
		$w.link.sel config -text "0"
	}
}

## @brief Deletes all components and links that are not selected.
#

proc deleteNonSelectedObjects { } {
	Bmodel delete_non_selected
	objectsDrawAll
	updateCompTypeTable
}

## @brief Sets radii of selected components.
#

proc setComponentRadius { } {
	global project_item
	global comp_radius comp_selection
	set comp_radius [.wmod.compradius.entry get]
    if { [string compare $comp_selection "Current"] == 0 } {
		set id [.wmod.curr.id cget -text]
		Bmodel set $project_item compradius $id $comp_radius
	} else {
		Bmodel set $project_item compradius $comp_radius
	}
	objectsDrawAll
}

## @brief Sets radii of selected links.
#

proc setLinkRadius { } {
	global project_item
	global link_radius
	set link_radius [.wmod.linkradius.entry get]
	Bmodel set $project_item linkradius $link_radius
	objectsDrawAll
}

## @brief Sets the selection of components.
#

proc setSelection args {
	global project_item comp_selected comp_selection types_selected
#	puts $comp_selection
    if { [string compare $comp_selection "All"] == 0 } {
		Bmodel component $project_item select all
	} elseif { [string compare $comp_selection "None"] == 0 } {
		Bmodel component $project_item select none
	} elseif { [string compare $comp_selection "Selected Types"] == 0 } {
		Bmodel component $project_item select types $types_selected
	}
	objectsDrawAll
	updateCompTypeTable
}

## @brief Sets the color for selected components.
#

proc setColor { } {
	global project_item color_selection selected_color types_selected
    if { [string compare $color_selection "All"] == 0 } {
		Bmodel component $project_item color all $selected_color
	} elseif { [string compare $color_selection "Selection"] == 0 } {
		Bmodel component $project_item color selected $selected_color
	} elseif { [string compare $color_selection "By Density"] == 0 } {
		Bmodel component $project_item color density
	} elseif { [string compare $color_selection "By FOM"] == 0 } {
		Bmodel component $project_item color fom
	}
	objectsDrawAll
}

## @brief Sets the views relative to the origin.
#

proc setViews { } {
	global project_item
	Bmodel set $project_item views
}

## @brief Extracts all panels tagged by the component list.
#

proc extractSegments { } {
	global filetypes fom_cutoff
	if ![Bmodel exists] { return }
	
	set xname [tk_getSaveFile -filetypes $filetypes \
		-initialfile "panels.tif" -defaultextension .tif]
	if { [string length $xname] < 1 } { return }
	
	set multi_level 1
	Bmodel extract_segments $xname $multi_level
}

## @brief Generates particles from the component list.
#

proc componentsToParticles { } {
	if ![Bmodel exists] { return }
	global project_item
	createMgParam
	Bmodel components_to_particles $project_item
}

## @brief Reads model from a file.
#

proc readModel { } {
	global filename
	global modelfile modeldir
	global comp_radius
#	if { ![Bimage exists $theimg] } {
#		set i [tk_dialog .dialog "No image!" \
#			{Please read an image associated with the model first.} \
#			info 0 {OK}]
#		return
#	}
	if ![Bmodel exists] { set modeldir [file dirname $filename] }
	set currdir $modeldir
	if { [string length $currdir] < 2 } { set currdir [pwd] }
	set currdir [pwd]
    set types {
		{"STAR files"		{.star .STAR}	}
		{"XML files"		{.xml .XML}	}
		{"PDB files"		{.pdb .PDB}	}
		{"CMM files"		{.cmm .CMM}	}
		{"Vega files"		{.v3d .V3D}	}
	}
	if { [string length $modelfile] < 1 } { set $modelfile "model.star" }
#	puts "$modeldir $modelfile ($currdir) ([pwd])"
#	puts [file pathtype $currdir]
#	puts [file exists $currdir]
#	puts [file system $currdir]
#	puts [file nativename $modelfile]
	set modelfile [tk_getOpenFile -title "Open model file" -filetypes $types \
			-initialfile $modelfile -initialdir $currdir -defaultextension .star]
#	puts "$modeldir $modelfile"
	if { [string length $modelfile] < 1 } { return }
	objectsDeleteAll
	loadModel $modelfile
}

# loadModel
# Read model into the model structures.
#
# @param	thefile			Model file name.

proc loadModel { thefile } {
	global filename project_item
	global modelfile modeldir
	global comp_radius
	set project_item "all"
	Bmodel read $thefile
	set ids [Bmodel component $project_item ids 0]
	foreach id $ids {
		set comp_radius [Bmodel component $project_item radius $id]
	}
	set modeldir [file dirname $thefile]
	set modelfile [file tail $thefile]
	selectAFile
	wm title . [list $filename ":" $modelfile]
	Model
}

## @brief Saves the component and link model to a file.
#

proc writeModel { } {
	global filename
	global project_item
	global modelfile modeldir
	set currdir $modeldir
	if { [string length $currdir] < 2 } { set currdir [pwd] }
#	set ncomp [Bmodel component $project_item count]
#	set nlink [Bmodel link $project_item count]
#	puts "Number of components and links: $ncomp $nlink"
	set model_name [file rootname [file tail $filename]]
	if { [string length $modelfile] < 1 } {
		set modelfile $model_name
		append modelfile "_model.star"
	}
	set modelfile [file tail $modelfile]
    set types {
		{"STAR files"		{.star .STAR}	}
		{"XML files"		{.xml .XML}	}
		{"PDB files"		{.pdb .PDB}	}
		{"CMM files"		{.cmm .CMM}	}
		{"Vega files"		{.v3d .V3D}	}
	}
	set modelfile [tk_getSaveFile -filetypes $types \
			-initialfile $modelfile -initialdir $currdir -defaultextension .star]
	if { [string length $modelfile] < 1 } { return }
	set modelfile [relativePath $modeldir $modelfile]
	if ![Bmodel exists] {
		set project_item [Bmodel create]
	}
	Bmodel set $project_item map $filename
	Bmodel write $modelfile
	wm title . [list $filename ":" $modelfile]
	writeSettingsFile
}

## @brief Edits a text parameter file.
#

proc editModel { } {
	global modelfile modeldir
	if [file exists $modelfile] {
		editTextFile $modelfile
	}
}

## @brief Dialog box for creating a shell (ring for 2D)
#

proc Shell { } {
	global helv12

	Model

	set w .wsph
	if ![winfo exists $w] {
		toplevel $w

		label $w.lab -font $helv12 -text \
			"Only two of the three parameters should be specified"

		setupEntry $w.num "Number" integer 0 "The number of components in the shell"
		setupEntry $w.rad "Radius (angstrom)" double 0 "The shell radius in angstrom"
		setupEntry $w.dist "Distance (angstrom)" double 0 "The distance between components in angstrom"

		frame $w.buttons
		button $w.buttons.do -text Do -command "doCreateShell $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons.do $w.buttons.close -side left -expand 1

		pack $w.lab $w.num $w.rad $w.dist $w.buttons -side top -fill x -pady 2m
		
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Shell"
    wm iconname $w "Shell"
#	.menuBar.window entryconfigure "Shell" -state normal
}

## @brief Analyzes model components to determine cyclic symmetry.
#

proc componentSymmetry { } {
	Model

	set w .wsym
	if ![winfo exists $w] {
		toplevel $w

		label $w.lab -text "Analyzing component symmetry"

		setupEntry $w.min "Minimum symmetry" integer 1 "The smallest cyclic symmetry to test for"
		setupEntry $w.max "Maximum symmetry" integer 20 "The largest cyclic symmetry to test for"
		setupEntry $w.ann "Annular width (pixels)" double 10 "The number of annuli in pixels to integrate over"

		frame $w.buttons
		button $w.buttons.do -text Do -command "doComponentSymmetry $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons.do $w.buttons.close -side left -expand 1

		pack $w.lab $w.min $w.max $w.ann $w.buttons -side top -fill x -pady 2m
		
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Symmetry"
    wm iconname $w "Symmetry"
}

proc doComponentSymmetry { w } {
	global env modelfile
	set min [$w.min.e get]
	set max [$w.max.e get]
	set ann [$w.ann.e get]
	if { ![file exists $modelfile] } {
		puts "$modelfile not found!"
	}
	set logfile "[file rootname $modelfile]_sym.log"
	set cmd "$env(BSOFT)/bin/bmodsym -verb 1 -all -find $min,$max -ann $ann $modelfile"
	executeCommand $cmd $logfile
}

