##
# @file		bshow_filament.tcl
#
# @brief	Procedures to create, move and destroy filament selections
#
# @author	Bernard Heymann
# @date		Created: 20020808
# @date		Modified: 20180410

set show_filaments 1
set show_node_labels 1
set show_profiles 0
set filament_width 40
set filament_node_radius 10
set boxing_interval 20
set helix_rise 0
set helix_angle 0
set node_selected []
set fil_plot_width 350
set fil_plot_height 200

## @brief Dialog box to manipulate filaments on the main canvas
#

proc Filaments { } {
	global project_item theimg imgtype
	global helv12
	global filament_width filament_node_radius boxing_interval
	global helix_rise helix_angle
	global fil_plot_width fil_plot_height

	if ![Bimage exists $theimg] {
		tk_dialog .dialog "No image in memory!" \
				{Please read in an appropriate image first.} \
				info 0 {OK}
		return
	}
	
	set c [getImageCanvas $theimg]
	if { $filament_width < 1 } { set filament_width 40 }
	if { $filament_node_radius < 1 } { set filament_node_radius [expr int($filament_width / 4)] }
	if { $filament_node_radius < 1 } { set filament_node_radius 3 }
	if { $boxing_interval < 1 } { set boxing_interval [expr int($filament_width / 2)] }
	if ![Bmg exists] {
		createMgParam
		setMgParam
	}
	puts "Launching Particles for $c: $project_item"
	
	set w .wfil
	if ![winfo exists $w] {
		toplevel $w
		menu $w.menuBar -tearoff 0
		$w.menuBar add cascade -menu $w.menuBar.fil -label "Filament" -underline 0
		menu $w.menuBar.fil -tearoff 0
		$w configure -menu $w.menuBar
		$w.menuBar.fil add command -label "Help" -underline 0 \
				-command { showURL "file://$Bsoft/doc/bshow/bshow_fil.html" }
#				-command { showHelp bshow_filament.hlp }
		$w.menuBar.fil add command -label "Write filament profile" \
				-command { filamentProfileWrite } -underline 0
		$w.menuBar.fil add command -label "Delete selected filament" \
				-command { filamentDelete } -underline 0
		$w.menuBar.fil add command -label "Delete all filaments" \
				-command { filamentDeleteAll } -underline 0
		$w.menuBar.fil add command -label "Center nodes" \
				-command { filamentCenterNodes } -underline 0
		$w.menuBar.fil add command -label "Extract filaments" \
				-command { filamentExtract } -underline 0
		$w.menuBar.fil add command -label "Generate particles" \
				-command { filamentsToParticles } -underline 0
		$w.menuBar.fil add command -label "Close" \
				-command "destroy $w" -underline 0
		
		frame $w.fil_width
		label $w.fil_width.tag -font $helv12 -text "Box width" -width 15 -anchor e
		entry $w.fil_width.entry -width 10 -validate key -vcmd { string is double %P }
		$w.fil_width.entry insert 0 $filament_width
		
		frame $w.node_radius
		label $w.node_radius.tag -font $helv12 -text "Node radius" -width 15 -anchor e
		entry $w.node_radius.entry -width 10 -validate key -vcmd { string is double %P }
		$w.node_radius.entry insert 0 $filament_node_radius

		frame $w.boxing_interval
		label $w.boxing_interval.tag -font $helv12 -text "Boxing interval" -width 15 -anchor e
		entry $w.boxing_interval.entry -width 10 -validate key -vcmd { string is double %P }
		$w.boxing_interval.entry insert 0 $boxing_interval

		frame $w.helix_rise
		label $w.helix_rise.tag -font $helv12 -text "Helix rise per unit" -width 15 -anchor e
		entry $w.helix_rise.entry -width 10 -validate key -vcmd { string is double %P }
		$w.helix_rise.entry insert 0 $helix_rise

		frame $w.helix_angle
		label $w.helix_angle.tag -font $helv12 -text "Helix angle per unit" -width 15 -anchor e
		entry $w.helix_angle.entry -width 10 -validate key \
			-vcmd { expr {[string match {[-+]} %P] || [string is double %P]} }
		$w.helix_angle.entry insert 0 $helix_angle

		frame $w.fil
		label $w.fil.tag -font $helv12 -text "Filaments" -width 15 -anchor e
		label $w.fil.num -width 10 -relief sunken -bd 1 -font $helv12 -anchor w -text "0"
		
		frame $w.node
		label $w.node.tag -font $helv12 -text "Nodes" -width 15 -anchor e
		label $w.node.num -width 10 -relief sunken -bd 1 -font $helv12 -anchor w -text "0"
		
		frame $w.fil_selected
		label $w.fil_selected.tag -font $helv12 -text "Selected filament" -width 15 -anchor e
		label $w.fil_selected.num -width 10 -relief sunken -bd 1 -font $helv12 -anchor w
		
		frame $w.node_selected
		label $w.node_selected.tag -font $helv12 -text "Selected node" -width 15 -anchor e
		label $w.node_selected.num -width 10 -relief sunken -bd 1 -font $helv12 -anchor w
		
		setupMicrographButtons $w.neighbor "micrograph"

		set g $w.plot
		canvas $g -width [expr $fil_plot_width + 10] -height [expr $fil_plot_height + 10] \
				-background white -relief sunken -borderwidth 2

		frame $w.switches
		checkbutton $w.switches.filaments -text "Show filaments" -variable show_filaments \
				-command { filamentsDrawAll }
		checkbutton $w.switches.labels -text "Show labels" -variable show_node_labels \
				-command { filamentsDrawAll }
		checkbutton $w.switches.profiles -text "Show profiles" -variable show_profiles
		
		pack $w.node_radius.tag $w.node_radius.entry -side left -pady 5 -padx 10
		pack $w.fil_width.tag $w.fil_width.entry -side left -pady 5 -padx 10
		pack $w.boxing_interval.tag $w.boxing_interval.entry -side left -pady 5 -padx 10
		pack $w.helix_rise.tag $w.helix_rise.entry -side left -pady 5 -padx 10
		pack $w.helix_angle.tag $w.helix_angle.entry -side left -pady 5 -padx 10
		pack $w.fil.tag $w.fil.num -side left -pady 5 -padx 10
		pack $w.node.tag $w.node.num -side left -pady 5 -padx 10
		pack $w.fil_selected.tag $w.fil_selected.num -side left -pady 5 -padx 10
		pack $w.node_selected.tag $w.node_selected.num -side left -pady 5 -padx 10
#		pack $w.neighbor.prev $w.neighbor.next \
#			-side left -pady 5 -padx 5 -ipadx 5 -ipady 5
		pack $w.switches.filaments $w.switches.labels $w.switches.profiles \
			 -side left -pady 5 -padx 5
		pack $w.fil_width $w.node_radius $w.boxing_interval $w.helix_rise \
			$w.helix_angle $w.fil $w.node \
			$w.fil_selected $w.node_selected $w.neighbor $w.plot $w.switches -side top -fill x -pady 2

	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Filaments"
    wm iconname $w "Filaments"
	.menuBar.window entryconfigure "Filaments" -state normal
	bind $c <1> "nodeCreateSelect %x %y"
	bind $c <2> "filamentCreate %x %y"
	bind $c <Command-1> "filamentCreate %x %y"
	bind $c <Control-1> "filamentCreate %x %y"
	bind $c <Shift-1> "nodeDelete %x %y"
	bind $c <Shift-2> "nodeDelete %x %y"
	bind $c <B1-Motion> "nodeMove %x %y"
	bind $c <B2-Motion> "nodeMove %x %y"
#	bind $c <ButtonRelease> "filamentsDrawAll"
#	bind $c <ButtonRelease> "drawSpline %W"	
	bind $w <Return> "Update 0"
	bind $w <Control-w> { destroy .wfil }
	bind $w <Command-w> { destroy .wfil }
#	updateNodeParam
	filamentsDrawAll
}

## @brief Sets the node radius from the dialog box.

proc setNodeRadius { } {
	if ![winfo exists .wfil] { return }
	global filament_node_radius project_item
	set filament_node_radius [.wfil.node_radius.entry get]
	Bmg set $project_item filament_node_radius $filament_node_radius
}

## @brief Starts a new filament
#
# @param	x
# @param	y			First node position in canvas coordinates.

proc filamentCreate { x y } {
	global node_selected
	set node_selected []
	nodeCreateSelect $x $y
}

## @brief Creates or selects a node.
#
# @param	x
# @param	y			Node position in canvas coordinates.

proc nodeCreateSelect { x y } {
	if ![winfo exists .wfil] { return }
	global theimg
	global node_selected
	global tool
	set snid $node_selected
	set node_selected [nodeSelect $x $y]
	if { [llength $node_selected] > 0 } { return }
	if { $tool != "filament" } { return }
	set wc [getControlWindow $theimg]
	setNodeRadius
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set loc [imageCoordinatesFromCanvas $x $y]
	set id [Bmg node $mg_item create [lindex $loc 0] [lindex $loc 1] [lindex $loc 2] [lindex $snid 0] [lindex $snid 1]]
	if { [llength $id] > 0 } {
		nodeDraw [lindex $id 0] [lindex $id 1]
#		puts "created ID = $id"
	}
	set node_selected $id
	filamentUpdateSelected
	filamentsDrawAll
#	filamentDraw [lindex $node_selected 0]
}

## @brief Selects a node on the canvas.
#
# @param	x
# @param	y				Node position in window coordinates.
#
# @returns	node_selected 	Index of the selected node, 0 if none selected.

proc nodeSelect { x y } {
	if ![winfo exists .wfil] { return }
	global theimg
	global node_selected
	set wc [getControlWindow $theimg]
	setNodeRadius
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set loc [imageCoordinatesFromCanvas $x $y]
#	puts $loc
	set id [Bmg node $mg_item select [lindex $loc 0] [lindex $loc 1] [lindex $loc 2]]
	if { [llength $id] > 0 } {
		nodeDraw [lindex $id 0] [lindex $id 1]
#		puts "selected ID = $id"
	}
	set node_selected $id
	filamentUpdateSelected
#	drawSelectedNode $c
	return $node_selected
}

## @brief Moves the selected node
#
# @param	x
# @param	y			Node position in canvas coordinates.

proc nodeMove { x y } {
	global theimg
	global node_selected
	global px py
	if { [llength $node_selected] < 1 } { return }
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set scale [$wc.scale.scale get]
	set n [$wc.image.scale get]
	set ix [expr [$c canvasx $x] / $scale ]
	set iy [expr [$c canvasy $y] / $scale ]
	set dx [expr $ix - $px]
	set dy [expr $py - $iy]
	set px $ix
	set py $iy
#	puts "Moving $dx $dy"
	set mg_item [micrographItem $n]
	Bmg node $mg_item move [lindex $node_selected 0] [lindex $node_selected 1] $dx $dy
	filamentsDrawAll
}

## @brief Draws a circle for a node on the canvas.
#
# @param	fid			Filament identifier.
# @param	nid			Node identifier.

proc nodeDraw { fid nid } {
	global theimg
	global node_selected show_node_labels
	global filament_color filament_first_color filament_select_color
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set z [$wc.slice.scale get]
	set mg_item [micrographItem $n]
	set loc [Bmg node $mg_item location $fid $nid]
#	puts "drawing $fid $nid: $loc"
	set scale [$wc.scale.scale get]
	set x [expr [lindex $loc 0] * $scale ]
	set y [expr [yFlip [lindex $loc 1]] * $scale ]
	set zr [expr abs([lindex $loc 2] - $z) ]
	set radius [expr [.wfil.node_radius.entry get] * $scale]
	set xmin [expr $x - $radius ]
	set xmax [expr $x + $radius + $scale - 1 ]
	set ymin [expr $y - $radius ]
	set ymax [expr $y + $radius + $scale - 1 ]
	set color $filament_color
	if { $nid == 1 } {
		set color $filament_first_color
		$c create text $xmin $y -text $fid -fill blue -anchor e -tags node
	}
	if { $fid == [lindex $node_selected 0] && $nid == [lindex $node_selected 1] } {
		set color $filament_select_color
	}
	$c create oval $xmin $ymin $xmax $ymax -outline $color -width 1 -tags node
	if { $show_node_labels } {
		$c create text $xmax $y -text $nid -fill $color -anchor w -tags node
	}
}

## @brief Draws a line between two nodes on the canvas.
#
# @param	fid			Filament identifier.
# @param	nid1		Node 1 identifier.
# @param	nid2		Node 2 identifier.

proc nodeLinkDraw { fid nid1 nid2 } {
	global theimg
	global filament_color
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set z [$wc.slice.scale get]
	set mg_item [micrographItem $n]
	set loc1 [Bmg node $mg_item location $fid $nid1]
	set loc2 [Bmg node $mg_item location $fid $nid2]
	set scale [$wc.scale.scale get]
	set x1 [expr [lindex $loc1 0] * $scale ]
	set y1 [expr [yFlip [lindex $loc1 1]] * $scale ]
	set x2 [expr [lindex $loc2 0] * $scale ]
	set y2 [expr [yFlip [lindex $loc2 1]] * $scale ]
	set fz -1.0
	set dz [expr [lindex $loc2 2] - [lindex $loc1 2]]
	if { $dz == 0 } {
		set fz [expr [lindex $loc2 2] - $z + 0.5]
	} else {
		set fz [expr ([lindex $loc2 2] - $z)/($dz*1.0)]
	}
	if { $fz >= 0.0 && $fz <= 1.0 } {
		$c create line $x1 $y1 $x2 $y2 -fill $filament_color -width 1 -tags fil
	}
}

## @brief Draws filament nodes and lines connecting nodes on the canvas
#

proc filamentsDrawAll { } {
	global theimg
	global show_filaments
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	$c delete fil node selnode
	if { $show_filaments < 1 } { return }
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set ids [Bmg node $mg_item ids]
#	puts "Node IDs: $ids"
	updateNodeParam
	set pfid 0
	set pnid 0
	foreach {fid nid} $ids {
		nodeDraw $fid $nid
		if { $pfid == $fid } {
			nodeLinkDraw $fid $pnid $nid
		}
		set pfid $fid
		set pnid $nid
	}
#	updateNodeParam
}

## @brief Updates the selected filament and node in the dialog window
#

proc filamentUpdateSelected { } {
	if ![winfo exists .wfil] { return }
	global node_selected project_item
#	puts "Project item: $project_item"
	set nfil [Bmg filament $project_item count]
	set nnode [Bmg node $project_item count]
	.wfil.fil.num config -text "$nfil"
	.wfil.node.num config -text "$nnode"
	if { [llength $node_selected] < 0 } {
		.wfil.fil_selected.num config -text ""
		.wfil.node_selected.num config -text ""
	} else {
		.wfil.fil_selected.num config -text [lindex $node_selected 0]
		.wfil.node_selected.num config -text [lindex $node_selected 1]
		filamentProfilePlot [lindex $node_selected 0] [lindex $node_selected 1]
	}
}

## @brief Deletes the selected node.
#
# @param	x
# @param	y			Node position in window coordinates.

proc nodeDelete { x y } {
	if ![winfo exists .wfil] { return }
	global theimg node_selected
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set id [nodeSelect $x $y]
#	puts "ID to be deleted = $id"
	if { [llength $id] > 0 } {
		Bmg node $mg_item delete [lindex $id 0] [lindex $id 1]
		set node_selected []
	}
	filamentsDrawAll
}

## @brief Deletes the selected filament.
#

proc filamentDelete { } {
	if ![winfo exists .wfil] { return }
	global theimg node_selected project_item
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	if { [llength $node_selected] > 0 } {
		Bmg filament $mg_item delete [lindex $node_selected 0]
		set node_selected []
	}
	filamentsDrawAll
}

## @brief Deletes all filaments.
#

proc filamentDeleteAll { } {
	global theimg project_item
	set c [getImageCanvas $theimg]
	Bmg filament $project_item delete all
	$c delete node selnode
	$c delete spline
	$c delete fil
	filamentUpdateSelected
}

## @brief Centers filament nodes.
#

proc filamentCenterNodes { } {
	global filament_width project_item
	set filament_width [.wfil.fil_width.entry get]
	puts $project_item
	Bmg filament $project_item center $filament_width
	filamentsDrawAll
}

## @brief Extracts filaments defined in the micrograph parameters.
#

proc filamentExtract { } {
	global helv12
	global filament_width project_item

	if { [Bmg node $project_item count] < 1 } {
		tk_messageBox -icon error -type ok -title "Error" -message \
			"No filaments picked!"
		return
	}
	
	set w .wfext
	catch {destroy $w}
	toplevel $w
	wm title $w "Extracting filaments"
	wm iconname $w "ExtFil"
	
	label $w.msg -font $helv12 -wraplength 4i -justify left \
		-text {Select parameters for filament extraction}
	
	frame $w.panel

	frame $w.panel.fil_width
	label $w.panel.fil_width.tag -font $helv12 -text "Filament width" -width 30 -anchor w
	entry $w.panel.fil_width.entry -width 10 -relief sunken \
		 -validate key -vcmd { string is double %P }
	$w.panel.fil_width.entry insert 0 $filament_width

	frame $w.panel.fil_axis
	label $w.panel.fil_axis.tag -font $helv12 -text "Filament axis" -width 30 -anchor w
	radiobutton $w.panel.fil_axis.x -text "x" -variable fil_axis -value 1
	radiobutton $w.panel.fil_axis.y -text "y" -variable fil_axis -value 2
	radiobutton $w.panel.fil_axis.z -text "z" -variable fil_axis -value 3
	
	frame $w.panel.path
	label $w.panel.path.tag -font $helv12 -text "Filament path" -width 10 -anchor w
	entry $w.panel.path.entry -width 50 -relief sunken
	$w.panel.path.entry insert 0 ""
	
	frame $w.check
	checkbutton $w.check.splitbutton -text "Individual output files" \
		-width 20 -anchor w -variable fil_split
	
	frame $w.buttons
	button $w.buttons.ok -text OK -command "filamentExtract_and_DestroyWindow $w"
	button $w.buttons.cancel -text Cancel -command "destroy $w"
		
	pack $w.msg -side top -padx 5 -pady 5
	
	pack $w.panel.fil_width.tag $w.panel.fil_width.entry -side left -padx 5 -pady 5
	pack $w.panel.fil_axis.tag $w.panel.fil_axis.x $w.panel.fil_axis.y \
		$w.panel.fil_axis.z -side left -padx 5 -pady 5
	pack $w.panel.path.tag $w.panel.path.entry -side left -padx 5 -pady 5
	pack $w.panel.fil_width $w.panel.fil_axis $w.panel.path -side top
	pack $w.panel -side top -fill both

	pack $w.check.splitbutton -side left -pady 5 -padx 5
	pack $w.check -side top -pady 2m
	
	pack $w.buttons.ok $w.buttons.cancel -side left -padx 5
	pack $w.buttons -side bottom -pady 2m

	bind $w <Return> "filamentExtract_and_DestroyWindow $w"
}

## @brief Extracts filaments from the current set of node coordinates.
#
# @param	w				Dialog window.

proc filamentExtract_and_DestroyWindow { w } {
	global filetypes
	global filename project_item
	global filament_width
	global fil_axis fil_split

	set xname [file root [file tail $filename]]
	append xname "_fil.mrc"
	set xname [tk_getSaveFile -filetypes $filetypes \
		-initialfile $xname -defaultextension .mrc]
	if { [string length $xname] < 1 } { return }
#	puts "Extract from $filename to $xname"
	set pathname [$w.panel.path.entry get]
	set filament_width [$w.panel.fil_width.entry get]
#	set fil_axis 0
#	puts "Filament axis = $fil_axis  split = $fil_split"
	.wfil.fil_width.entry delete 0 end
	.wfil.fil_width.entry insert 0 $filament_width
	Bmg set $project_item pixel_size [getPixelSizeEntry]
	setMgParam
	Bmg filament $project_item extract $xname $filament_width $fil_axis $fil_split $pathname
	destroy $w
}

## @brief Creates particles on filaments.
#

proc filamentsToParticles { } {
	global filament_width filament_node_radius boxing_interval
	global helix_rise helix_angle project_item
#	puts "toParticles"
	set filament_width [.wfil.fil_width.entry get]
	set filament_node_radius [.wfil.node_radius.entry get]
	set boxing_interval [.wfil.boxing_interval.entry get]
	set helix_rise [.wfil.helix_rise.entry get]
	set helix_angle [.wfil.helix_angle.entry get]
	setMgParam
	Bmg set $project_item filament_width $filament_width
	Bmg set $project_item filament_node_radius $filament_node_radius
#	puts "$filament_width $boxing_interval $helix_rise $helix_angle"
	set np [Bmg filament $project_item toparticles $filament_width $boxing_interval $helix_rise $helix_angle]
#	puts "$np particles generated"
	Particles
}

## @brief Updates node parameters from the micrograph parameters in memory.
#

proc updateNodeParam { } {
	global project_item imgtype
	global filament_width filament_node_radius
	
	set w .wfil
	if ![winfo exists $w] { return }
#	puts "Updating filament parameters"
	if ![Bmg exists] { createMgParam }
	set filename [Bmg get $project_item filename $imgtype]
#	puts "File: $filename"
	wm title $w [list "Filaments:" $filename]
	updatePixelsize
	set filament_width [$w.fil_width.entry get]
	set filament_node_radius [$w.node_radius.entry get]
	Bmg set $project_item filament_width $filament_width
	Bmg set $project_item filament_node_radius $filament_node_radius
#	puts "$filament_width $filament_node_radius"
	$w.fil_width.entry delete 0 end
	$w.fil_width.entry insert 0 $filament_width
	$w.node_radius.entry delete 0 end
	$w.node_radius.entry insert 0 $filament_node_radius
	filamentUpdateSelected
}

## @brief Plots the selected filament profile in the plot canvas
#
# @param	fid			filament identifier.
# @param	nid			node identifier.

proc filamentProfilePlot { fid nid } {
	global theimg project_item
	global filament_width
	global fil_plot_width fil_plot_height
	global show_profiles
	
	if ![winfo exists .wfil] { return }
	if { $fid < 1 } { return }
	if { $nid < 1 } { return }

	set g .wfil.plot
	$g delete filplot

	if !$show_profiles { return }

	set prof [Bmg filament $project_item profile $fid $nid $filament_width]
#	puts $prof
	set len [lindex $prof 0]
	if { $len < 1 } { return }
	
	set min [lindex $prof 1]
	set max [lindex $prof 1]
	for {set i 1} {$i <= $len} {incr i 1} {
		set v [lindex $prof $i]
		if { $min > $v } { set min $v }
		if { $max < $v } { set max $v }
	}
	if { $max <= $min } { set max [expr $min + 1] }
#	puts "profile: $min $max"
	set xoff 5
	set yoff [expr $fil_plot_height * ( (3*$max - $min) * 0.25 / ($max - $min)) + 5]
	set xscale [expr $fil_plot_width*1.0/$len]
	set yscale [expr $fil_plot_height*0.5/($max - $min)]
	
#	puts "$yoff $yscale"
	
	set fp ""
	set fps ""
	for {set i 1} {$i <= $len} {incr i 1} {
		set j [expr $i + $len]
		set x [expr ($i-1)*$xscale + $xoff]
		set y [expr $yoff - [lindex $prof $i] * $yscale]
		set ys [expr $yoff - [lindex $prof $j] * $yscale]
		append fp " " $x " " $y
		append fps " " $x " " $ys
	}
	$g create line $fp -fill blue -smooth yes -width 1 -tags filplot
	$g create line $fps -fill red -smooth yes -width 1 -tags filplot
	set orix [expr $fil_plot_width/2 + $xoff]
	set oriy $yoff
	$g create line $orix 0 $orix [expr 2*$yoff + $fil_plot_height] \
		-fill black -smooth yes -width 1 -tags filplot
	$g create line 0 $oriy [expr 2*$xoff + $fil_plot_width] $oriy \
		-fill black -smooth yes -width 1 -tags filplot
}

## @brief Writes the selected filament profile to a file
#

proc filamentProfileWrite { } {
	global theimg project_item
	global node_selected
	global filament_width
	if { [llength $node_selected] < 0 } {
		tk_dialog .dialog "No filament selected!" \
				{Please select a filament first.} \
				info 0 {OK}
		return
	}
	set filename [Bimage get $theimg filename]
	set fprfile "[file rootname [file tail $filename]].txt"
	set fprfile [tk_getSaveFile -initialfile $fprfile -defaultextension .txt]
	if { [string length $fprfile] < 1 } { return }
	set wc [getControlWindow $theimg]
	set pixel_size [$wc.pixel_size.x get]
	set fid [lindex $node_selected 0]
	set nid [lindex $node_selected 1]
	set prof [Bmg filament $project_item profile $fid $nid $filament_width]
	set len [lindex $prof 0]
	if { $len < 1 } { return }
	set ffpr [open $fprfile w]
	puts $ffpr "File name:              $filename"
	puts $ffpr "Filament:               $fid"
	puts $ffpr "Pixel size:             $pixel_size angstrom/pixel"
	puts $ffpr " "
	puts $ffpr "Pixel Density"
	for {set i 1} {$i <= $len} {incr i 1} {
		set v [lindex $prof $i]
		puts $ffpr "$i $v"
	}
	puts $ffpr " "
	close $ffpr
}

