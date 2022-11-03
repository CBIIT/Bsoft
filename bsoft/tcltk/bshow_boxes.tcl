##
# @file		bshow_boxes.tcl
#
# @brief	Procedures to create, move and destroy particle boxes
#
# @author	Bernard Heymann
# @date		Created: 20010516
# @date		Modified: 20220406

set show_boxes 1
set show_bad 1
set show_box_labels 1
set part_sel 1
set auto_renum 1
set box_color_mode "Uniform"
set fom_min 0
set fom_max 1
set sel_min 0
set sel_max 1
set box_size "100 100 100"
set bad_radius 20
set contrast 0
set box_selected 0
set scale 1
set thetemp "thetemp"
set peach 0

## @brief Dialog box to manipulate particle boxes on the main canvas
#

proc Particles { } {
	global project_item theimg imgtype
	global box_size bad_radius box_selected part_sel mg_sel
	
	if ![Bimage exists $theimg] {
		tk_dialog .dialog "No image in memory!" \
				{Please read in an appropriate image first.} \
				info 0 {OK}
		return
	}
	
	set c [getImageCanvas $theimg]
	set filename [Bimage get $theimg filename]
	if ![Bmg exists] {
		createMgParam
		setMgParam
	}
#	puts "Launching Particles for $c: $project_item"
	if { [string length $imgtype] < 1 } { set imgtype "mg" }
	if { [lindex $box_size 0] < 1 } { set box_size "100 100 100" }
	if { $bad_radius < 1 } { set bad_radius [expr int([lindex $box_size 0] / 4)] }
#	puts "Box size: $box_size"
#	puts "Bad radius: $bad_radius"

	set typelabel ""
	if { $imgtype == "mg" } {
		set typelabel "micrograph"
	} elseif { $imgtype == "frame" } {
		set typelabel "frames"
	} elseif { $imgtype == "rec" } {
		set typelabel "reconstruction"
	} elseif { $imgtype == "part" } {
		set typelabel "particle"
		if { $box_selected <= 0 } {
			set box_selected 1
		}
	}
#	puts "$imgtype: $box_selected"
		
	set w .wbox
	if ![winfo exists $w] {
		if { $imgtype != "mg" && $imgtype != "rec" && $imgtype != "part" && $imgtype != "frame" } {
			set i [tk_dialog .dialog "No micrograph, reconstruction or particle!" \
				{Please read in an appropriate image first.} \
				info 0 {Cancel} {Proceed anyway}]
			if { $i < 1 } { return }
		}
		
		toplevel $w
		
		menu $w.menuBar -tearoff 0
		$w.menuBar add cascade -menu $w.menuBar.box -label "Boxes" -underline 0
		menu $w.menuBar.box -tearoff 0
		$w configure -menu $w.menuBar
		$w.menuBar.box add command -label "Help" -underline 0 \
				-command { showURL "file://$Bsoft/doc/bshow/bshow_part.html" }
#				-command { showHelp bshow_boxes.hlp }
		$w.menuBar.box add command -label "Delete all boxes" \
				-command { boxDeleteAll } -underline 0
		$w.menuBar.box add command -label "Delete boxes in current image" \
				-command { boxDeleteCurrent } -underline 0
		$w.menuBar.box add command -label "Delete unselected boxes" \
				-command { boxDeleteUnselected } -underline 0
		$w.menuBar.box add command -label "Renumber particles" \
				-command { boxRenumber } -underline 0
		$w.menuBar.box add command -label "Center particles" \
				-command { boxCenter } -underline 0
		$w.menuBar.box add command -label "Template picker" \
				-command { templatePicker } -underline 0
		$w.menuBar.box add command -label "Variance picker" \
				-command { variancePicker } -underline 0
		$w.menuBar.box add command -label "Extract particles from current image" \
				-command { extractParticles } -underline 0
		$w.menuBar.box add command -label "Show selected particles" \
				-command { showSelectedParticles } -underline 0
		$w.menuBar.box add command -label "Close" \
				-command "destroy $w" -underline 0

		$w.menuBar add cascade -menu $w.menuBar.spaflow -label "Workflow" -underline 0
		menu $w.menuBar.spaflow -tearoff 0
		$w configure -menu $w.menuBar
		$w.menuBar.spaflow add command -label "Extract particles" \
				-command { extractParticles } -underline 0
		$w.menuBar.spaflow add command -label "Align particles" \
				-command { alignParticles } -underline 0
		$w.menuBar.spaflow add command -label "Refine particles" \
				-command { refineParticles } -underline 0
		$w.menuBar.spaflow add command -label "Align 3D particles" \
				-command { align3DParticles } -underline 0
		$w.menuBar.spaflow add command -label "Select particles" \
				-command { selectParticles } -underline 0
		$w.menuBar.spaflow add command -label "Reconstruct" \
				-command { particleReconstruct } -underline 0
		
		setupItemID $w.itemid
		
		if [regexp "Micrograph" $project_item] {
			setupMicrographThickness $w.thick
		}
		
		setupVectorEntry $w.box_size "Box size" integer $box_size "Particle box size in pixels/voxels"
		setupEntry $w.bad_radius "Bad radius" double $bad_radius "Bad area radius in pixels/voxels"

		frame $w.number
		label $w.number.tag -text "Particles picked: mg" -width 20 -anchor w
		label $w.number.num -width 6 -text "0" -relief sunken -bd 1 -anchor w
		label $w.number.tag2 -text "project" -width 10 -anchor w
		label $w.number.num2 -width 6 -text "0" -relief sunken -bd 1 -anchor w
		pack $w.number.tag $w.number.num $w.number.tag2 $w.number.num2 -side left -padx 1m

		frame $w.nbad
		label $w.nbad.tag -text "Bad regions" -width 20 -anchor w
		label $w.nbad.num -width 10 -text "0" -relief sunken -bd 1 -anchor w
		pack $w.nbad.tag $w.nbad.num -side left -padx 1m

		frame $w.sel
		label $w.sel.tag -text "Current particle" -width 20 -anchor w
		label $w.sel.id -width 6 -text "0" -relief sunken -bd 1 -anchor w
		label $w.sel.fom -width 6 -text "0" -relief sunken -bd 1 -anchor w
#		checkbutton $w.sel.set -text "Selection" -variable part_sel \
#				-command { updatePartSel }
		tk_optionMenu $w.sel.set part_sel 0 1 2 3 4 5 6 7 8 9 10 \
			11 12 13 14 15 16 17 18 19 20
		pack $w.sel.tag $w.sel.id $w.sel.fom $w.sel.set -side left -padx 1m

		frame $w.fom
		label $w.fom.label -text "FOM" -anchor w
		scale $w.fom.scale -orient horizontal -length 400 -from 0 -to 1 \
				-command { setFOMcutoff } -resolution 0.001
		pack $w.fom.label $w.fom.scale -side top -fill x -padx 1m

		frame $w.switches
		checkbutton $w.switches.boxes -text "Show boxes" \
			-variable show_boxes -command { boxesDrawAll }
		checkbutton $w.switches.labels -text "Show labels" \
			-variable show_box_labels -command { boxesDrawAll }
		checkbutton $w.switches.bad -text "Show bad areas" \
			-variable show_bad -command { boxesDrawAll }
		pack $w.switches.boxes $w.switches.labels $w.switches.bad \
			-side left -ipadx 5 -ipady 5 -padx 1m
		
		frame $w.switches2
		checkbutton $w.switches2.autorenum -text "Automatic renumbering" -variable auto_renum \
				-command { boxesDrawAll }
		tk_optionMenu $w.switches2.color box_color_mode \
			"Uniform" "FOM" "Selection"
		pack $w.switches2.autorenum $w.switches2.color \
			-side left -ipadx 5 -ipady 5 -padx 1m

		setupMicrographButtons $w.neighbor $typelabel

		pack $w.itemid -side top -fill x -pady 2
		if [regexp "Micrograph" $project_item] {
			pack $w.thick -side top -fill x -pady 2
		}
		pack $w.box_size $w.bad_radius $w.number $w.nbad $w.sel $w.fom \
			$w.switches $w.switches2 $w.neighbor -side top -fill x -pady 2

		updateBoxParam
	} else {
		wm deiconify $w
		raise $w
    }
	wm title $w [list "Particles:" $filename]
    wm iconname $w "Particles"
	.menuBar.window entryconfigure "Particles" -state normal
	bind $c <1> "boxCreateSelect %x %y part"
	bind $c <Shift-1> "boxDelete %x %y"
	bind $c <Control-1> "boxCreateSelect %x %y bad"
	bind $c <2> "boxCreateSelect %x %y bad"
	bind $c <B1-Motion> "boxMove %x %y"
	bind $c <B2-Motion> "boxMove %x %y"
	bind $w <Return> "Update 0"
	bind $w <Control-w> { destroy $w }
	bind $w <Command-w> { destroy $w }
	trace add variable part_sel write updatePartSel
	boxesDrawAll
}

## @brief Sets the 3 dimensions for a box from the entries.
#

proc setBoxSize { } {
	set w .wbox
	if ![winfo exists $w] { return }
	global theimg project_item
	global box_size bad_radius
#	puts "Setting box size for $theimg: $project_item"
#	set wc [getControlWindow $theimg]

	set box_size_x [$w.box_size.x get]
	set box_size_y [$w.box_size.y get]
	set box_size_z [$w.box_size.z get]
	set box_size "$box_size_x $box_size_y $box_size_z"
	if { [lindex $box_size 0] < 1 } { set box_size "50 50 50" }
	set bad_radius [$w.bad_radius.e get]
#	puts "setBoxSize [$wc.slice.scale cget -to]"

	if { [Bimage get $theimg nslices] < 2 } {
		set box_size [lreplace $box_size 2 2 1]
		$w.box_size.z delete 0 end
		$w.box_size.z insert 0 "1"
	}

	Bmg set $project_item box_size [lindex $box_size 0] [lindex $box_size 1] [lindex $box_size 2]
	Bmg set $project_item bad_radius $bad_radius
#	puts "setBoxSize: $box_size"
}

## @brief Selects a box or bad area near a clicked location
#
# @param	x
# @param	y			Box position in canvas coordinates.

proc boxSelect { x y } {
	if ![winfo exists .wbox] { return }
	global project_item theimg
	global box_selected
	global tool
	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set mg_item [micrographItem $img_num]
	setBoxSize
	set loc [imageCoordinatesFromCanvas $x $y]
	set id [Bmg box $mg_item select [lindex $loc 0] [lindex $loc 1] [lindex $loc 2]]
	boxDrawOldNew $id
	return $id
}

## @brief Creates or selects a box or bad area near a clicked location
#
# @param	x
# @param	y			Box position in canvas coordinates.
# @param	t			Type: "part"/"bad"

proc boxCreateSelect { x y t } {
	if ![winfo exists .wbox] { return }
	global project_item theimg
	global box_selected part_sel
	global tool
	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set mg_item [micrographItem $img_num]
	setBoxSize
	set loc [imageCoordinatesFromCanvas $x $y]
	set id [Bmg box $mg_item select [lindex $loc 0] [lindex $loc 1] [lindex $loc 2]]
#	puts "selected ID = $id"
	if { $id == 0 && $tool == "box" } {
		set id [Bmg box $mg_item create [lindex $loc 0] [lindex $loc 1] [lindex $loc 2] $t $part_sel]
#		puts "created ID = $id"
	}
	boxDrawOldNew $id
	return $id
}

proc boxDrawOldNew { id } {
	if ![winfo exists .wbox] { return }
	global box_selected
	if { $id != $box_selected } {
		set tid $box_selected
		set box_selected $id
		boxDraw $tid
	}
	if { $id != 0 } {
		boxDraw $id
		updateBoxSelected
		updateBoxCounts
#		puts "selected ID = $id"
	}
}

## @brief Moves the selected box or bad area to the new coordinates.
#
# @param	x
# @param	y			Box position in canvas coordinates.

proc boxMove { x y } {
	if ![winfo exists .wbox] { return }
	global project_item theimg
	global box_selected
	global px py
	if { $box_selected == 0 } { return }
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set mg_item [micrographItem $img_num]
	set scale [$wc.scale.scale get]
	set ix [expr [$c canvasx $x] / $scale ]
	set iy [expr [$c canvasy $y] / $scale ]
	set dx [expr $ix - $px]
	set dy [expr $py - $iy]
	set px $ix
	set py $iy
#	puts "Moving $dx $dy"
	Bmg box $mg_item move $box_selected $dx $dy
	boxesDrawAll
}


## @brief Draws a box or bad area
#
# Draws a box and circle with the current box size and number
# for a selected box, or draws a red circle with the current bad size
# for a selected bad area.
#
# @param	id			Box or bad area identifier.

proc boxDraw { id } {
	global project_item theimg
	global box_selected box_color_mode
	global show_box_labels box_color box_select_color bad_color bad_select_color
	global fom_min fom_max
	global sel_min sel_max
	if ![winfo exists .wbox] { return }
	if { $id == 0 } { return }
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set z [$wc.slice.scale get]
	set mg_item [micrographItem $n]
	set loc [Bmg box $mg_item location $id]
#	puts "Location for $mg_item $id: $loc"
	if { [lindex $loc 0] < 0 } {
		puts "ID $id not found!"
		return
	}
	set scale [$wc.scale.scale get]
	set scale_1 [expr $scale - 1]
	set yi [yFlip [lindex $loc 1]]
	set xr [expr [lindex $loc 0] * $scale ]
	set yr [expr $yi * $scale ]
	set zr [expr abs([lindex $loc 2] - $z) ]
	if { $id > 0 } {
		set box_radius_x [expr [.wbox.box_size.x get] * 0.5]
		set box_radius_y [expr [.wbox.box_size.y get] * 0.5]
		set box_radius_z [expr [.wbox.box_size.z get] * 0.5]
		if { $zr > $box_radius_z } { return }
		set zr2 1
		if { $box_radius_z > 0 } {
			set zr [expr $zr*1.0/$box_radius_z]
			set zr2 [expr sqrt(1 - $zr*$zr)]
		}
#		puts "$xr $yr"
		set box_radius_x [expr $box_radius_x * $scale ]
		set box_radius_y [expr $box_radius_y * $scale ]
		set oval_radius_x [expr $box_radius_x*$zr2]
		set oval_radius_y [expr $box_radius_y*$zr2]
		set xmin [expr $xr - $box_radius_x ]
		set xmax [expr $xr + $box_radius_x + $scale_1 ]
		set ymin [expr $yr - $box_radius_y ]
		set ymax [expr $yr + $box_radius_y + $scale_1 ]
		set cxmin [expr $xr - $oval_radius_x ]
		set cxmax [expr $xr + $oval_radius_x + $scale_1 ]
		set cymin [expr $yr - $oval_radius_y ]
		set cymax [expr $yr + $oval_radius_y + $scale_1 ]
		set idtag [format "box%d" $id]
   		if { [string compare $box_color_mode "FOM"] == 0 } {
			set fom [Bmg box $mg_item fom $id]
			set color [spectrum $fom $fom_min $fom_max]
   		} elseif { [string compare $box_color_mode "Selection"] == 0 } {
			set sel [Bmg box $mg_item select $id]
			set color [spectrum $sel $sel_min $sel_max]
		} elseif { [lsearch $box_selected $id] >= 0 } {
			set color $box_select_color;
		} else {
			set color $box_color;
		}
		$c create rect $xmin $ymin $xmax $ymax -outline $color -width 1 -tags [list box $idtag]
		$c create oval $cxmin $cymin $cxmax $cymax -outline $color -width 1 -tags [list box $idtag]
		if { $show_box_labels } {
			$c create text $xr $yr -text $id -fill $color -tags [list box $idtag]
		}
	} elseif { $id < 0 } {
		set bad_radius [.wbox.bad_radius.e get]
		if { $zr > $bad_radius } { return } 
		set bad_radius [expr $scale * sqrt($bad_radius*$bad_radius - $zr*$zr)]
		set xmin [expr $xr - $bad_radius ]
		set xmax [expr $xr + $bad_radius + $scale_1 ]
		set ymin [expr $yr - $bad_radius ]
		set ymax [expr $yr + $bad_radius + $scale_1 ]
		if { [lsearch $box_selected $id] >= 0 } {
			set color $bad_select_color;
		} else {
			set color $bad_color;
		}
		$c create oval $xmin $ymin $xmax $ymax -outline $color -width 1 -tags bad
	}
}

## @brief Redraws all boxes and bad circles.
#

proc boxesDrawAll { } {
	global project_item theimg
	global imgtype fom_cutoff
	global show_boxes show_bad
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
#	puts "Drawing boxes for $theimg: $project_item"
	$c delete box bad
	if { [string compare $imgtype "mg"] && [string compare $imgtype "rec"] } { return }
	if ![winfo exists .wbox] { return }
	setBoxSize
	updateBoxCounts
	if { $show_boxes } {
#		puts "Getting box IDs for cutoff = $fom_cutoff"
		set ids [Bmg box $mg_item ids $fom_cutoff]
		foreach id $ids {
			boxDraw $id
		}
	}
#	puts "Boxes drawn"
	if { $show_bad } {
		set ids [Bmg box $mg_item ids bad]
		foreach id $ids {
			boxDraw $id
		}
	}
#	puts "Done"
}


## @brief Deletes the selected box or bad area.
#
# @param	x
# @param	y			Box or bad area position in canvas coordinates.

proc boxDelete { x y } {
	if ![winfo exists .wbox] { return }
	global theimg project_item
	global boxes_selected
	global auto_renum
	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set mg_item [micrographItem $img_num]
	set id [boxSelect $x $y]
#	puts "ID to be deleted = $id"
	if { $id != 0 } {
		Bmg box $mg_item delete $id
		if { $auto_renum } { Bmg box $mg_item renumber }
	}
	boxesDrawAll
}

## @brief Deletes all boxes.
#

proc boxDeleteAll { } {
	global project_item theimg
	Bmg box $project_item delete all
	set c [getImageCanvas $theimg]
	$c delete box bad
	.wbox.number.num config -text "0"
	.wbox.nbad.num config -text "0"
}

## @brief Deletes all boxes in current image.
#

proc boxDeleteCurrent { } {
	global project_item theimg
	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set mg_item [micrographItem $img_num]
	Bmg box $mg_item delete current
	boxesDrawAll
}

## @brief Deletes all unselected boxes.
#

proc boxDeleteUnselected { } {
	global auto_renum fom_cutoff project_item
#	puts "boxDeleteUnselected: $project_item $fom_cutoff"
	Bmg box $project_item delete fom $fom_cutoff
	if { $auto_renum } { Bmg box $project_item renumber }
	boxesDrawAll
}

## @brief Renumbers boxes.
#

proc boxRenumber { } {
	global project_item
	Bmg box $project_item renumber
	boxesDrawAll
}

## @brief Centers boxes.
#

proc boxCenter { } {
	global project_item
	Bmg box $project_item center
	boxesDrawAll
}

## @brief Picks particles by cross-correlation with a template.
#

proc templatePicker { } {
	if ![winfo exists .wbox] { return }
	global theimg thetemp
	global window_scale hires lores
	
	set w .wtemp
	if ![winfo exists $w] {
		toplevel $w
		
		set m $w.menu
		menu $m -tearoff 0
		$m add cascade -menu $m.file -label "File" -underline 0
		menu $m.file -tearoff 0
		$w configure -menu $m
		$m.file add command -label "Open ..." -command { readTemplate } \
			-underline 0 -accelerator "Ctrl-o"
		$m.file add command -label "Save ..." -command { saveTemplate } \
			-underline 0 -accelerator "Ctrl-S"

		frame $w.frame -relief ridge -borderwidth 2 -bg black

		wm minsize $w 200 200
		set c $w.frame.c
		canvas $c -width 200 -height 200 \
			-relief sunken -borderwidth 2 \
			-xscrollcommand "$w.frame.hscroll set" \
			-yscrollcommand "$w.frame.vscroll set"
		scrollbar $w.frame.vscroll -orient vertical -command "$c yview"
		scrollbar $w.frame.hscroll -orient horizontal -command "$c xview"
			
		grid $c -in $w.frame \
			-row 0 -column 0 -rowspan 1 -columnspan 1 -sticky news
		grid $w.frame.vscroll \
			-row 0 -column 1 -rowspan 1 -columnspan 1 -sticky news
		grid $w.frame.hscroll \
			-row 1 -column 0 -rowspan 1 -columnspan 1 -sticky news
#		grid $c -in $w.frame -sticky nw -padx 2 -pady 2
		grid rowconfig    $w.frame 0 -weight 1 -minsize 0
		grid columnconfig $w.frame 0 -weight 1 -minsize 0

		image create photo $thetemp
		$c create image 0 0 -image $thetemp -anchor nw -tags template
		
		frame $w.scale
		label $w.scale.label -text "Scale" -width 6 -anchor w
		scale $w.scale.scale -orient horizontal -from 0.1 -to 10 \
			-command { setTemplateScale } -resolution 0.1
		pack $w.scale.label -side left
		pack $w.scale.scale -side bottom -fill x
		$w.scale.scale set $window_scale

		setupResolutionEntry $w

		labelframe $w.fomcut -text "FOM range"
		setupEntry $w.fomcut.min "Minimum" double 0 \
			"The minimum FOM cutoff to select picked particles\nIf it is zero, an automatic cutoff will be calculated"
		setupEntry $w.fomcut.max "Maximum" double 1e30 \
			"The maximum FOM cutoff to select picked particles"
		pack $w.fomcut.min $w.fomcut.max -side top -expand 1

		setupEntry $w.excld "Exclusion distance" double 0 "Minimum distance between picked particles"
		setupEntry $w.bin "Binning" integer 1 "Binning the image to accelerate picking"
		
		frame $w.buttons
		button $w.buttons.update -text Update -command "updateBoxTemplate $w" -relief raised 
		button $w.buttons.pick -text Pick -command "doTemplatePicker $w" -relief raised 
		button $w.buttons.cancel -text Close -command  " destroy $w " -relief raised 
		pack $w.buttons.update $w.buttons.pick $w.buttons.cancel -side left -expand 1
		
		pack $w.frame $w.scale $w.res $w.fomcut $w.excld $w.bin $w.buttons -side top -fill both -expand yes
		
		wm title $w "Template Picker"
		wm iconname $w "Template"
	} else {
		wm deiconify $w
		raise $w
    }
	
	setTemplateScale 1
	updateBoxTemplate $w
	
	bind $w <Control-w> { destroy $w }
	bind $w <Command-w> { destroy $w }
}

proc setTemplateScale { scale } {
	global thetemp
#	global window_width window_height window_scale
	set w .wtemp
	if ![winfo exists $w] { return }
	if ![Bimage exists thetemp] { return }
	set c $w.frame.c
	set width [Bimage get $thetemp width]
	set height [Bimage get $thetemp height]
	set size "0 0 [expr $width * $scale ] [expr $height * $scale ]"
	$c configure -scrollregion $size
#	$c configure -width $window_width -height $window_height -scrollregion $size
#	puts "set scale to $scale"
	
	set oz [expr [Bimage get $thetemp nslices] / 2]
	
#	Bimage template $thetemp 0 $oz $scale 1
	Bimage template $thetemp $oz $scale
}

proc readTemplate { } {
	global imgdir thetemp
	global filetypes
	set w .wtemp
	if ![winfo exists $w] { return }
	set currdir $imgdir
	if { [string length $currdir] < 2 } { set currdir [pwd] }
	set new_filename [tk_getOpenFile -title "Open image file" \
			-filetypes $filetypes -initialdir $currdir]
	if { [string length $new_filename] < 1 } { return }
	if { [catch { set nimages [Bimage open $thetemp $new_filename 0] } err] } {
		bgerror $err
		return
	}
	set oz [expr [Bimage get $thetemp nslices] / 2]
	set scale [$w.scale.scale get]
#	Bimage template $thetemp 0 $oz 1 1
	Bimage template $thetemp $oz $scale
}

proc saveTemplate { } {
	global imgdir thetemp
	global filetypes filename
	set currdir $imgdir
	if { [string length $currdir] < 2 } { set currdir [pwd] }
	set new_filename [tk_getSaveFile -filetypes $filetypes \
		-initialfile [file tail $filename] -initialdir $currdir \
		-defaultextension .mrc]
	if { [string length $new_filename] < 1 } { return }
    Bimage set $thetemp pixel_size [getPixelSizeEntry]
	Bimage save $thetemp $new_filename
	writeSettingsFile
}

## @brief Updates and shows the particle template.
#
# @param	w 			the window.

proc updateBoxTemplate { w } {
	if ![winfo exists .wbox] { return }
	global theimg thetemp
	
	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set mg_item [micrographItem $img_num]

	Bmg box $mg_item template

	if ![winfo exists $w] {
		templatePicker
	}
	
#	set z [$wc.slice.scale get]
	if { [Bimage exists $thetemp] } {
		set oz [expr [Bimage get $thetemp nslices] / 2]
		set scale [$w.scale.scale get]
#		Bimage template $thetemp 0 $oz $scale 1
		Bimage template $thetemp $oz $scale
	}
}

## @brief Picks particles using cross-correlation.
#
# @param	w 			the window.

proc doTemplatePicker { w } {
	global project_item theimg
	global hires lores fom_cutoff

	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set mg_item [micrographItem $img_num]

	set hires [$w.res.hi.e get]
	set lores [$w.res.lo.e get]
	set fom_cutoff [$w.fomcut.min.e get]
	set fommax [$w.fomcut.max.e get]
	set excld [$w.excld.e get]
	set bin [$w.bin.e get]

	Bmg box $mg_item pickcc $hires $lores $fom_cutoff $fommax $excld $bin

	updateBoxParam
	boxesDrawAll
}

## @brief Picks particles from a variance map.
#

proc variancePicker { } {
	if ![winfo exists .wbox] { return }
	global theimg thetemp
	
	set w .wvar
	if ![winfo exists $w] {
		toplevel $w
		
		setupEntry $w.avg "Average kernel" integer 21 "The edge size of the averaging kernel"
		setupEntry $w.var "Variance kernel" integer 301 "The edge size of the variance kernel"
		setupEntry $w.excld "Exclusion distance" double 100 "Minimum distance between picked particles"
		setupEntry $w.bin "Binning" integer 1 "Binning the image to accelerate picking"
#		setupEntry $w.nsig "Relative to sigma" double 5 "Number of standard deviations from the average in the variance map"
		setupEntry $w.fmin "Minimum cutoff" double 0.1 "Minimum to accept in the variance map"
		setupEntry $w.fmax "Maximum cutoff" double 0.2 "Maximum to accept in the variance map"
		
		frame $w.buttons
		button $w.buttons.pick -text Pick -command "doVariancePicker $w" -relief raised 
		button $w.buttons.cancel -text Close -command  " destroy $w " -relief raised 
		pack $w.buttons.pick $w.buttons.cancel -side left -expand 1
		
		pack $w.avg $w.var $w.excld $w.bin $w.fmin $w.fmax $w.buttons -side top -fill both -expand yes
		
		wm title $w "Variance Picker"
		wm iconname $w "Variance"
	} else {
		wm deiconify $w
		raise $w
    }
	
	bind $w <Control-w> { destroy $w }
	bind $w <Command-w> { destroy $w }
}

## @brief Picks particles using variance map.
#
# @param	w 			the window.

proc doVariancePicker { w } {
	global project_item theimg

	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set mg_item [micrographItem $img_num]

	set avg [$w.avg.e get]
	set var [$w.var.e get]
	set excld [$w.excld.e get]
	set bin [$w.bin.e get]
#	set nsig [$w.nsig.e get]
	set fmin [$w.fmin.e get]
	set fmax [$w.fmax.e get]

	Bmg box $mg_item pickvar $avg $var $excld $bin $fmin $fmax

	updateBoxParam
	boxesDrawAll
}

proc updatePartSel args {
	global project_item
	global box_selected part_sel
#	puts "Particle $box_selected: $part_sel"
	Bmg box $project_item select $box_selected $part_sel
#	boxesDrawAll
}

## @brief Updates box and bad area parameters.
#

proc updateBoxParam { } {
	global theimg imgtype
	global project_item
	global box_size bad_radius
	global fom_min fom_max
	global sel_min sel_max
	set w .wbox
	if ![winfo exists $w] { return }
	if { $imgtype != "mg" && $imgtype != "rec" && $imgtype != "part" && $imgtype != "frame" } { return }
#	puts "Updating box parameters for $theimg: $project_item"
	if ![Bmg exists] {
		createMgParam
		setBoxSize
	}
#	puts "Project item $project_item $idx"
	updatePixelsize
	set temp_size [Bmg get $project_item box_size]
	if { [lindex $temp_size 0] > 0 && [lindex $temp_size 1] > 0 } {
		set box_size $temp_size
	}
	if { [Bmg get $project_item bad_radius] > 0 } {
		set bad_radius [Bmg get $project_item bad_radius]
	}
	set filename [Bmg get $project_item filename $imgtype]
	wm title $w [list "Particles:" $filename]
	updateItemID $w.itemid
	updateMicrographThickness $w.thick
	$w.box_size.x delete 0 end
	$w.box_size.x insert 0 [lindex $box_size 0]
	$w.box_size.y delete 0 end
	$w.box_size.y insert 0 [lindex $box_size 1]
	$w.box_size.z delete 0 end
	$w.box_size.z insert 0 [lindex $box_size 2]
	$w.bad_radius.e delete 0 end
	$w.bad_radius.e insert 0 $bad_radius
	updateBoxCounts
	updateBoxSelected
	set fom_min [Bmg box all fom_min]
	set fom_max [Bmg box all fom_max]
	set sel_min [Bmg box all select_min]
	set sel_max [Bmg box all select_max]
}

## @brief Updates box and bad area counts.
#

proc updateBoxCounts { } {
	global project_item theimg
	global fom_max
	if ![winfo exists .wbox] { return }
#	puts "Updating box counts for $project_item"
	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set mg_item [micrographItem $img_num]
	set nbox [Bmg box $mg_item count]
	set nbad [Bmg box $mg_item count bad]
	set nboxall [Bmg box all count]
	set w .wbox
	if [winfo exists $w] {
		$w.number.num config -text "$nbox"
		$w.number.num2 config -text "$nboxall"
		$w.nbad.num config -text "$nbad"
		set fom_max [Bmg box $project_item fom_max]
		if { $fom_max < 0.001 } { set fom_max 0.001 }
		set sres 0.00001
		if { $fom_max > 0.005 } {
			set sres 0.0001
		} elseif { $fom_max > 0.05 } {
			set sres 0.001
		} elseif { $fom_max > 0.5 } {
			set sres 0.01
		}
		$w.fom.scale config -to $fom_max -tickinterval [expr $fom_max / 5] -resolution $sres
	}
}

## @brief Updates selected box and bad area counts.
#

proc updateBoxSelected { } {
	global theimg imgtype project_item
	global box_selected part_sel
	set w .wbox
	if ![winfo exists $w] { return }
	set wc [getControlWindow $theimg]
	if { $imgtype == "part" } {
		set box_selected [expr [$wc.image.scale get] + 1]
#		puts $box_selected
	}
	if { $box_selected > 0 } {
		$w.sel.id config -text "$box_selected"
		set fom [Bmg box $project_item fom $box_selected]
		$w.sel.fom config -text "$fom"
		set sel [Bmg box $project_item select $box_selected]
		if { $sel >= 0 } {
			set part_sel $sel
		}
	}
#	puts "$box_selected $part_sel"
}

## @brief Extracts particles from the current set of box coordinates.
#

proc extractParticles { } {
	global project_item theimg
	global scale

	if { [Bmg box $project_item count] < 1 } {
		tk_messageBox -icon error -type ok -title "Error" -message \
			"No particles picked!"
		return
	}
	
	set w .wpext
	catch {destroy $w}
	toplevel $w
	wm title $w "Extracting particles"
	wm iconname $w "ExtPart"
	
	label $w.msg -wraplength 4i -justify left \
		-text {Select parameters for particle extraction}
	
	frame $w.fill
	label $w.fill.tag -text "Fill value" -width 12 -anchor w
	radiobutton $w.fill.avg -text "Average" -variable filltype -value 1
	radiobutton $w.fill.back -text "Background" -variable filltype -value 2
	radiobutton $w.fill.val -text "Value" -variable filltype -value 0
	entry $w.fill.entry -width 10 -relief sunken \
		 -validate key -vcmd { string is double %P }
	$w.fill.entry insert 0 [Bimage get $theimg background]
#	$w.fill.entry insert 0 0
	
	setupEntry $w.scale "Scale" double $scale "Particle extraction scale relative to the micrograph sampling"

	setupEntry $w.mask_width "Filament mask width" double 0.0 "Mask width to erase areas outside filament"
	
	frame $w.path
	label $w.path.tag -text "Particle path" -width 10 -anchor w
	entry $w.path.entry -width 50 -relief sunken
	$w.path.entry insert 0 ""
	
	frame $w.check
#	checkbutton $w.check.oddbutton -text "Odd size" -width 20 -anchor w -variable odd
	checkbutton $w.check.backbutton -text "Correct background" -width 20 -anchor w -variable back
	checkbutton $w.check.normbutton -text "Normalize" -width 14 -anchor w -variable norm
	checkbutton $w.check.splitbutton -text "Individual output files" \
		-width 20 -anchor w -variable part_split
	
	frame $w.buttons
	button $w.buttons.ok -text OK -command "doExtractParticles $w"
	button $w.buttons.cancel -text Cancel -command "destroy $w"
		
	pack $w.msg -side top -padx 5 -pady 5
	
	pack $w.fill.tag $w.fill.avg $w.fill.back $w.fill.val \
		$w.fill.entry -side left -fill both -expand true
	pack $w.path.tag $w.path.entry -side left -padx 5 -pady 5
	pack $w.fill $w.scale $w.mask_width $w.path -side top

	pack $w.check.backbutton $w.check.normbutton $w.check.splitbutton -side left -pady 5 -padx 2
	pack $w.check -side top -pady 2m
	
	pack $w.buttons.ok $w.buttons.cancel -side left -padx 5
	pack $w.buttons -side bottom -pady 2m

	bind $w <Return> "doExtractParticles $w"
}

## @brief Extracts particles from the current set of box coordinates.
#
# @param	w				Dialog window.

proc doExtractParticles { w } {
	global filetypes project_item
	global filename imgdir paramdir
	global scale filltype
	global back norm part_split
	set xdir $paramdir
	if { [string length $xdir] < 2 } { set xdir $imgdir }	
	if { [string length $xdir] < 2 } { set xdir [pwd] }	
	set xname [file root [file tail $filename]]
	append xname "_part.mrc"
	set xname [tk_getSaveFile -filetypes $filetypes \
		-initialfile $xname -initialdir $xdir -defaultextension .mrc]
	if { [string length $xname] < 1 } { return }
#	puts "Extract from $filename to $xname"
	set xname [relativePath [pwd] $xname]
	set pathname [$w.path.entry get]
	set scale [$w.scale.e get]
	set mask_width [$w.mask_width.e get]
	set fill [$w.fill.entry get]
	Bmg set $project_item pixel_size [getPixelSizeEntry]
	setMgParam
	Bmg box $project_item extract $xname $scale $back $norm $filltype $fill $mask_width $part_split $pathname
	destroy $w
}

## @brief Aligns particles by projection matching to a reference.
#

proc alignParticles { } {
	global paramfile
	global hires lores 
	global angle annmin annmax shift_limit
	global peach
	checkParamFile
	set w .waln
	if ![winfo exists $w] {
		toplevel $w
		
		set angle 3
		set shift_limit 10
		set annmin 1
		set annmax 100
		
		set newfile [file rootname $paramfile]
		append newfile "_aln.star"
		set logfile [file rootname $paramfile]
		append logfile "_aln.log"
		
		frame $w.head
		label $w.head.lab -text \
			"Aligning will be done on the currently saved parameter file: $paramfile"
		pack $w.head.lab

		setupNewFileEntry $w.file "Output parameter file" $newfile \
			"New parameter file written after alignment"

		setupScrollEntry $w.ref "Reference map" "" "Template or reference map to align to"

		setupEntry $w.sym "Symmetry" alnum "C1" "Point group symmetry"

		setupResolutionEntry $w

		setupEntry $w.ang "Angular step" double $angle \
			"The grid search step intervals in degrees"
		setupEntry $w.shift "Shift limit" double $shift_limit \
			"The maximum number of pixels of shift allowed"

		labelframe $w.ann -text "Annuli"
		setupEntry $w.ann.min "Minimum" double $annmin \
			"The innermost annulus to include"
		setupEntry $w.ann.max "Maximum" double $annmax \
			"The outermost annulus to include"
		pack $w.ann.min $w.ann.max -side top -expand 1

		checkbutton $w.peach -text "Peach" -variable peach

		pack $w.head $w.file $w.ref $w.sym $w.res $w.ang \
			$w.shift $w.ann -side top -pady 2 -padx 10

		frame $w.buttons
		button $w.buttons.do -text Do -command "doAlignParticles $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doAlignParticles $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Align"
    wm iconname $w "Align"
}

## @brief Executes marker tracking.
#
# @param	w 			dialog window.

proc doAlignParticles { w } {
	global env paramfile
	global hires lores 
	global peach
	set newfile [$w.file.e get]
	set ref [$w.ref.e get]
	set sym [$w.sym.e get]
	set hires [$w.res.hi.e get]
	set lores [$w.res.lo.e get]
	set angle [$w.ang.e get]
	set shift [$w.shift.e get]
	set annmin [$w.ann.min.e get]
	set annmax [$w.ann.max.e get]
	setMgParam
	if { ![file exists $paramfile] } {
		puts "$paramfile not found!"
	}
	set logfile "[file rootname $newfile].log"
#	set psfile "[file rootname $newfile].ps"
#	puts "logfile = $logfile"
    if { $peach } {
        set cmd "psubmit /stachraid/bin/bomrun4 -mode ccc -CTF -sym $sym -resol $hires,$lores -ang $angle -shift $shift -ann $annmin,$annmax -ref $ref $paramfile -ss"
   } else {
        set cmd "$env(BSOFT)/bin/borient -verb 1 -mode ccc -CTF -sym $sym -resol $hires,$lores -ang $angle -shift $shift -ann $annmin,$annmax -ref $ref -out $newfile $paramfile"
        executeCommand $cmd $logfile
    }
}

## @brief Refines particle orientations by central section matching to a reference.
#

proc refineParticles { } {
	global paramfile
	global hires lores 
	global angini angacc stpini stpacc def mag
	checkParamFile
	set w .wref
	if ![winfo exists $w] {
		toplevel $w
		
		set angini 3
		set angacc 0.5
		set stpini 4
		set stpacc 0.5
		set def 0.0
		set mag 0.0
		
		set newfile [file rootname $paramfile]
		append newfile "_aln.star"
		set logfile [file rootname $paramfile]
		append logfile "_aln.log"
		
		frame $w.head
		label $w.head.lab -text \
			"Aligning will be done on the currently saved parameter file: $paramfile"
		pack $w.head.lab

		setupNewFileEntry $w.file "Output parameter file" $newfile \
			"New parameter file written after alignment"

		setupScrollEntry $w.ref "Reference map" "" "Template or reference map to align to"

		labelframe $w.grd -text "Angular grid search"
		setupEntry $w.grd.ini "Initial angle" double $angini "Initial angluar step size in degrees"
		setupEntry $w.grd.acc "Angular accuracy" double $angacc "Target angular accuracy in degrees"
		pack $w.grd.ini $w.grd.acc -side top -expand 1

		labelframe $w.stp -text "Translational search"
		setupEntry $w.stp.ini "Initial step" double $stpini "Initial translational step size in pixels"
		setupEntry $w.stp.acc "Translational accuracy" double $stpacc "Target translational accuracy in pixels"
		pack $w.stp.ini $w.stp.acc -side top -expand 1

		setupResolutionEntry $w
		
		setupEntry $w.def "Defocus refinement" double $def \
			"Defocus step size to refine defocus in micrometer"
		setupEntry $w.mag "Magnification refinement" double $mag \
			"The maximum number of pixels of shift allowed"

		checkbutton $w.peach -text "Peach" -variable peach

		pack $w.head $w.file $w.ref $w.grd $w.stp $w.res $w.def $w.mag \
			-side top -pady 2 -padx 10

		frame $w.buttons
		button $w.buttons.do -text Do -command "doRefineParticles $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doRefineParticles $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Refine"
    wm iconname $w "Refine"
}

## @brief Executes particle orientation refinement.
#
# @param	w 			dialog window.

proc doRefineParticles { w } {
	global env paramfile
	global hires lores 
	global angini angacc stpini stpacc def mag
	set newfile [$w.file.e get]
	set ref [$w.ref.e get]
	set angini [$w.grd.ini.e get]
	set angacc [$w.grd.acc.e get]
	set stpini [$w.stp.ini.e get]
	set stpacc [$w.stp.acc.e get]
	set hires [$w.res.hi.e get]
	set lores [$w.res.lo.e get]
	set def [$w.def.e get]
	set mag [$w.mag.e get]
	setMgParam
	if { ![file exists $paramfile] } {
		puts "$paramfile not found!"
	}
	set opt "-grid $angini,$angacc "
	if { $stpacc } { append opt "-step $stpini,$stpacc " }
	append opt "-resol $hires,$lores "
	if { $def } { append opt "-defocus $def " }
	if { $mag } { append opt "-magnif $mag " }
	
	set logfile "[file rootname $newfile].log"
#	set psfile "[file rootname $newfile].ps"
#	puts "logfile = $logfile"
	set cmd "$env(BSOFT)/bin/brefine -verb 1 $opt -ref $ref -out $newfile $paramfile"
	executeCommand $cmd $logfile
}

## @brief Aligns 3D particles to a reference.
#

proc align3DParticles { } {
	global paramfile
	global hires lores
	global angle angacc shift_limit
	global peach
	checkParamFile
	set w .waln3
	if ![winfo exists $w] {
		toplevel $w
		
		set angle 5
		set angacc 1
		set shift_limit 10
		
		set newfile [file rootname $paramfile]
		append newfile "_aln.star"
		set logfile [file rootname $paramfile]
		append logfile "_aln.log"
		
		frame $w.head
		label $w.head.lab -text \
			"Aligning will be done on the currently saved parameter file: $paramfile"
		pack $w.head.lab

		setupNewFileEntry $w.file "Output parameter file" $newfile \
			"New parameter file written after alignment."

		setupScrollEntry $w.ref "Reference map" "" "Template or reference map to align to."

		setupScrollEntry $w.mask "Mask" "" "Frequency space mask defining missing regions."

		setupEntry $w.sym "Symmetry" alnum "C1" "Point group symmetry"

		setupResolutionEntry $w

		labelframe $w.grd -text "Angular grid search"
		setupEntry $w.grd.ini "Initial angle" double $angle "Initial angluar step size in degrees"
		setupEntry $w.grd.acc "Angular accuracy" double $angacc "Target angular accuracy in degrees"
		pack $w.grd.ini $w.grd.acc -side top -expand 1

		setupEntry $w.shift "Shift limit" double $shift_limit \
			"The maximum number of pixels of shift allowed"

		checkbutton $w.peach -text "Peach" -variable peach

		pack $w.head $w.file $w.ref $w.mask $w.sym $w.res $w.grd \
			$w.shift -side top -pady 2 -padx 10

		frame $w.buttons
		button $w.buttons.do -text Do -command "doAlign3DParticles $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doAlign3DParticles $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Align"
    wm iconname $w "Align"
}

## @brief Executes marker tracking.
#
# @param	w 			dialog window.

proc doAlign3DParticles { w } {
	global env paramfile
	global hires lores
	global angle angacc shift_limit
	global peach
	set newfile [$w.file.e get]
	set mask [$w.mask.e get]
	set ref [$w.ref.e get]
	set sym [$w.sym.e get]
	set hires [$w.res.hi.e get]
	set lores [$w.res.lo.e get]
	set angle [$w.grd.ini.e get]
	set angacc [$w.grd.acc.e get]
	set shift [$w.shift.e get]
	setMgParam
	if { ![file exists $paramfile] } {
		puts "$paramfile not found!"
	}
	set logfile "[file rootname $newfile].log"
#	set psfile "[file rootname $newfile].ps"
#	puts "logfile = $logfile"
    if { $peach } {
   } else {
        set cmd "$env(BSOFT)/bin/3daln.pl -verb 1 -global 2 -mode refine -sym $sym -resol $hires,$lores -ang $angle -acc $angacc -shift $shift -Temp $ref -Mask $mask -out $newfile $paramfile"
		executeCommand $cmd $logfile
    }
}

## @brief Selects particles after an alignment run.
#

proc selectParticles { } {
	global paramfile theimg project_item
	global selall fom_check cv_check top_check ran_check defadj
	showSelectedParticles
	set w .wspasel
	if ![winfo exists $w] {
		toplevel $w

		set newfile [file rootname $paramfile]
		append newfile "_sel.star"
		
		frame $w.head
		label $w.head.lab -text \
			"Particle seletion will be done using the currently saved parameter file: $paramfile"
		
		setupNewFileEntry $w.file "Output parameter file" $newfile \
			"New parameter file written after selection"

		checkbutton $w.all -text "Select all" -variable selall
		setupCheckEntry $w.fom "FOM" list fom fom_check "Selection based on the figure-of-merit"
		setupCheckEntry $w.cv "CV" list cv cv_check "Selection based on cross-validation"
		setupCheckEntry $w.top "Top" list top top_check "Selection based on the top scorers"
		checkbutton $w.defadj -text "Adjust for defocus" -variable defadj
		setupCheckEntry $w.ran "Random" list ran ran_check "Random selection"

		pack $w.head.lab
		pack $w.head $w.file $w.all $w.fom $w.cv $w.top $w.defadj $w.ran -side top -pady 2 -padx 15

		frame $w.buttons
		button $w.buttons.do -text Do -command "doSelectParticles $w"
		button $w.buttons.reset -text Reset -command "showSelectedParticles"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.reset $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doSelectParticles $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Select"
    wm iconname $w "Select"
}

## @brief Shows the current particle selection.
#

proc showSelectedParticles { } {
	global env paramfile
	set logfile "[file rootname $paramfile].log"
	set cmd "$env(BSOFT)/bin/bpartsel -verb 7 $paramfile"
	executeCommand $cmd $logfile
}

## @brief Executes particle selection.
#
# @param	w 			dialog window.

proc doSelectParticles { w } {
	global env paramfile
	global selall fom_check cv_check top_check ran_check defadj
	set newfile [$w.file.e get]
	set fom [$w.fom.e get]
	set cv [$w.cv.e get]
	set top [$w.top.e get]
	set ran [$w.ran.e get]
	setMgParam
	if { ![file exists $paramfile] } {
		puts "$paramfile not found!"
	}
	set ad ""
	if { $defadj } { set ad ",1" }
	set opt ""
	if { $selall } { append opt "-all " }
	if { $fom_check } { append opt "-fom $fom$ad " }
	if { $cv_check } { append opt " -cv $cv$ad " }
	if { $top_check } { append opt "-top $top$ad " }
	if { $ran_check } { append opt "-random $ran " }
	set logfile "[file rootname $newfile].log"
	set cmd "$env(BSOFT)/bin/bpartsel -verb 7 $opt -out $newfile $paramfile"
	executeCommand $cmd $logfile
}

## @brief Reconstructs a single particle.
#

proc particleReconstruct { } {
	global paramfile theimg project_item
	global hires lores
	set w .wsparec
	if ![winfo exists $w] {
		toplevel $w

		set sz [Bmg get $project_item box_size]
		set x [lindex $sz 0]
		set y [lindex $sz 1]
		set z [lindex $sz 0]

		set recfile [file rootname $paramfile]
		append recfile "_rec.mrc"
		set newfile [file rootname $paramfile]
		append newfile "_rec.star"
		
		frame $w.head
		label $w.head.lab -text \
			"Reconstruction will be done using the currently saved parameter file: $paramfile"
		
		frame $w.rec
		label $w.rec.l -text "Output reconstruction" -width 20 -anchor w
		entry $w.rec.e -width 40
		$w.rec.e insert 0 $recfile
		pack $w.rec.l $w.rec.e -side left -pady 5

		setupNewFileEntry $w.file "Output parameter file" $newfile \
			"New parameter file written after reconstruction"

		setupEntry $w.cls "Classes" alnum "all" "Classes of particles to select for reconstruction"

		setupEntry $w.interp "Type" integer 0 "Reconstruction interpolation type: 0=nearest neighbor, 1=weighted, 2=trilinear"

		labelframe $w.size -text "Reconstruction size"
		setupEntry $w.size.x "x" double $x "Reconstruction x size"
		setupEntry $w.size.y "y" double $y "Reconstruction y size"
		setupEntry $w.size.z "z" double $z "Reconstruction z size"
		pack $w.size.x $w.size.y $w.size.z -side top -expand 1

		setupEntry $w.scale "Scale" double 1 "Reconstruction scale"
		setupEntry $w.res "High resolution" double $hires "Spatial frequency limit for reconstruction"
		setupEntry $w.ctf "CTF correction" alnum "flip" "Contrast transfer function correction algorithm"
		setupEntry $w.sym "Symmetry" alnum "C1" "Point group symmetry"

		pack $w.head.lab
		pack $w.head $w.rec $w.file $w.cls $w.interp $w.size $w.scale $w.res $w.ctf $w.sym -side top -fill x -pady 2 -padx 15

		frame $w.buttons
		button $w.buttons.do -text Do -command "doReconstruct $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doReconstruct $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Reconstruct"
    wm iconname $w "Reconstruct"
}

## @brief Executes tomogram reconstruction.
#
# @param	w 			dialog window.

proc doReconstruct { w } {
	global env paramfile
	global hires
	set recfile [$w.rec.e get]
	set newfile [$w.file.e get]
	set cls [$w.cls.e get]
	set interp [$w.interp.e get]
	set x [$w.size.x.e get]
	set y [$w.size.y.e get]
	set z [$w.size.z.e get]
	set scale [$w.scale.e get]
	set hires [$w.res.e get]
	set ctf [$w.ctf.e get]
	set sym [$w.sym.e get]
	if { $ctf == "none" } {
		set ctf ""
	} else {
		set ctf "-CTF $ctf"
	}
	setMgParam
	if { ![file exists $paramfile] } {
		puts "$paramfile not found!"
	}
	set logfile "[file rootname $newfile].log"
	set cmd "$env(BSOFT)/bin/breconstruct -verb 1,tim -classes $cls -interp $interp -size $x,$y,$z -scale $scale -resol $hires $ctf -sym $sym -rescale 0,1 -full -half -rec $recfile -out $newfile $paramfile"
	executeCommand $cmd $logfile
}

