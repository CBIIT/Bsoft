##
# @file		bshow_xtal.tcl
#
# @brief	Procedures to manage crystallographic images
#
# @author	Bernard Heymann
# @date		Created: 20061104
# @date		Modified: 20181030	

global unitcell
set show_spot 1
set show_spot_indices 1
set spot_radius 5
set spot_resolution 20
set spot_selected []

## @brief Dialog box to manage crystallographic information
#

proc Crystallography { } {
	global theimg imgtype project_item
	global helv12
	global spot_radius

	if ![Bimage exists $theimg] {
		tk_dialog .dialog "No image in memory!" \
				{Please read in an appropriate image first.} \
				info 0 {OK}
		return
	}
	
	if { [string first "ft" $imgtype] < 0 } {
		if { [string first "ps" $imgtype] < 0 } {
			set i [tk_dialog .dialog "No power spectrum!" \
				{No Fourier transform or power spectrum was found. Please transform the image first.} \
				info 0 {Cancel} {Proceed anyway}]
			if { $i < 1 } { return }
			if { [Bimage get $theimg channels] > 1 } {
				set imgtype "ft"
			} else {
				set imgtype "ps"
			}
		}
	}
	
	set c [getImageCanvas $theimg]

	if ![Bmg exists] {
		createMgParam
	}
	
	set cen [Bimage get $theimg center]
	Bmg set $project_item origin [lindex $cen 0] [lindex $cen 1] [lindex $cen 2]
	
	set w .wxtal
	if ![winfo exists $w] {
		toplevel $w
		menu $w.menuBar -tearoff 0
		$w.menuBar add cascade -menu $w.menuBar.tomo -label "Crystallography" -underline 0
		menu $w.menuBar.tomo -tearoff 0
		$w configure -menu $w.menuBar
		$w.menuBar.tomo add command -label "Help" -underline 0 \
				-command { showURL "file://$Bsoft/doc/bshow/bshow_xtal.html" }
#				-command { showHelp bshow_xtal.hlp }
		$w.menuBar.tomo add command -label "Set spot indices" \
				-command { setSpotIndices } -underline 0
		$w.menuBar.tomo add command -label "Calculate unit cell vectors" \
				-command { calcUnitCellVectors } -underline 0
		$w.menuBar.tomo add command -label "Generate spots" \
				-command { spotGenerate } -underline 0
		$w.menuBar.tomo add command -label "Mask spots" \
				-command { spotMask } -underline 0
		$w.menuBar.tomo add command -label "Delete all spots" \
				-command { spotsDeleteAll } -underline 0
		$w.menuBar.tomo add command -label "Close" \
				-command { destroy $w } -underline 0
		
		frame $w.spot_radius
		label $w.spot_radius.tag -font $helv12 -text "Structure factor radius" -width 20 -anchor e
		entry $w.spot_radius.entry -width 10 -validate key -vcmd { string is double %P }
		$w.spot_radius.entry insert 0 $spot_radius
		
		frame $w.number
		label $w.number.tag -font $helv12 -text "Structure factor count" -width 20 -anchor e
		label $w.number.num -width 10 -text "0" -relief sunken -bd 1 \
				-font $helv12 -anchor w

		frame $w.selected
		label $w.selected.tag -font $helv12 -text "Selected" -width 15 -anchor e
		label $w.selected.h -width 10 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w
		label $w.selected.k -width 10 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w
		label $w.selected.l -width 10 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w
		label $w.selected.fom -width 10 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w

		frame $w.switches
		checkbutton $w.switches.spots -text "Show spots" -variable show_spot \
				-command { spotsDrawAll }
		checkbutton $w.switches.labels -text "Show indices" -variable show_spot_indices \
				-command { spotsDrawAll }
		
		frame $w.unitcell
		label $w.unitcell.tag -font $helv12 -text "Unit cell" -width 10 -anchor e
		label $w.unitcell.hx -width 6 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w
		label $w.unitcell.hy -width 6 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w
		label $w.unitcell.kx -width 6 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w
		label $w.unitcell.ky -width 6 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w
		
		frame $w.frame
		text $w.table -relief sunken -bd 2 -width 40\
			-yscrollcommand "$w.yscroll set" \
			-setgrid 1 -highlightthickness 0 -pady 2 -padx 3
		scrollbar $w.yscroll -command "$w.table yview" \
			-highlightthickness 0 -orient vertical
		grid $w.table -in $w.frame -padx 1 -pady 1 \
			-row 0 -column 0 -rowspan 1 -columnspan 1 -sticky news
		grid $w.yscroll -in $w.frame -padx 1 -pady 1 \
			-row 0 -column 1 -rowspan 1 -columnspan 1 -sticky news
		grid rowconfig    $w.frame 0 -weight 1 -minsize 0
		grid columnconfig $w.frame 0 -weight 1 -minsize 0

		pack $w.spot_radius.tag $w.spot_radius.entry -side left -pady 5 -padx 10
		pack $w.number.tag $w.number.num -side left -pady 5 -padx 10
		pack $w.selected.tag $w.selected.h $w.selected.k $w.selected.l \
				$w.selected.fom -side left -pady 5 -padx 10
		pack $w.switches.spots $w.switches.labels -side left -ipadx 5 -ipady 5 -pady 5 -padx 5
		pack $w.unitcell.tag $w.unitcell.hx $w.unitcell.hy \
				$w.unitcell.kx $w.unitcell.ky -side left -pady 5 -padx 5
		pack $w.spot_radius $w.number $w.selected $w.switches \
				$w.unitcell -side top -pady 5 -padx 10
#		pack $w.frame -expand yes -fill both -padx 5 -pady 5
		
#		updateSFTable
		updateSpotParam		
		spotsDrawAll
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Crystallography"
    wm iconname $w "Crystallography"
	.menuBar.window entryconfigure "Crystallography" -state normal
	bind $c <1> "spotCreateSelect %x %y"
	bind $c <2> "spotCreateSelect %x %y"
	bind $c <Shift-1> "spotDelete %x %y"
	bind $c <Shift-2> "spotDelete %x %y"
	bind $c <B1-Motion> "spotMove %x %y"
	bind $c <B2-Motion> "spotMove %x %y"
	bind $w <Return> "Update 0"
	bind $w <Control-w> { destroy $w }
	bind $w <Command-w> { destroy $w }
#	updateSFTable
}

## @brief Dialog box to set Miller indices for the selected structure factor

proc setSpotIndices { } {
	global helv12
	global spot_selected
	set h 0
	set k 0
	set l 0
	if { [llength $spot_selected] > 0 } {
		set h [lindex $spot_selected 0]
		set k [lindex $spot_selected 1]
		set l [lindex $spot_selected 2]
	}
	set w .wind
	if ![winfo exists $w] {
		toplevel $w
		
		frame $w.index
		label $w.index.tag -font $helv12 -text "Indices" -width 14
		entry $w.index.hentry -width 6 -validate key \
			-vcmd { expr {[string match {[-+]} %P] || [string is int %P]} }
		entry $w.index.kentry -width 6 -validate key \
			-vcmd { expr {[string match {[-+]} %P] || [string is int %P]} }
		entry $w.index.lentry -width 6 -validate key \
			-vcmd { expr {[string match {[-+]} %P] || [string is int %P]} }
		$w.index.hentry insert 0 $h
		$w.index.kentry insert 0 $k
		$w.index.lentry insert 0 $l

		pack $w.index.tag $w.index.hentry $w.index.kentry $w.index.lentry -side left -padx 5 -pady 5
		pack $w.index -side top -fill x -pady 2 -padx 15
		
		frame $w.buttons
		button $w.buttons.ok -text OK -command "setSpotIndices_and_DestroyWindow $w"
		button $w.buttons.cancel -text Cancel -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.ok $w.buttons.cancel -side left -expand 1
	
		bind $w <Return> "setSpotIndices_and_DestroyWindow $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Indices"
    wm iconname $w "Indices"
	tkwait window $w
}

## @brief Sets the indices for a structure factor.
#
# @param	w 			the window to destroy.

proc setSpotIndices_and_DestroyWindow { w } {
	global theimg
	global spot_selected
	set wc [getControlWindow $theimg]
	set h [$w.index.hentry get]
	set k [$w.index.kentry get]
	set l [$w.index.lentry get]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	if { [llength $spot_selected] > 0 } {
		set spot_selected [Bmg spot $mg_item reindex [lindex $spot_selected 0] \
			[lindex $spot_selected 1] [lindex $spot_selected 2] $h $k $l]
	} else {
		set spot_selected "$h $k $l"
	}
	.wxtal.selected.h config -text [lindex $spot_selected 0]
	.wxtal.selected.k config -text [lindex $spot_selected 1]
	.wxtal.selected.l config -text [lindex $spot_selected 2]
	destroy $w
	spotsDrawAll
}

## @brief Calculates the unit cell vectors from the current structure factors.
#

proc calcUnitCellVectors { } {
	global project_item
	setMgParam
	Bmg spot $project_item unitcell
	set unitcell [Bmg get $project_item unitcell]
	.wxtal.unitcell.hx config -text [lindex $unitcell 0]
	.wxtal.unitcell.hy config -text [lindex $unitcell 1]
	.wxtal.unitcell.kx config -text [lindex $unitcell 2]
	.wxtal.unitcell.ky config -text [lindex $unitcell 3]
}

## @brief Select a structure factor.
#
# @param	x
# @param	y			Structure factor position in canvas coordinates.

proc spotSelect { x y } {
	if ![winfo exists .wxtal] { return }
	global theimg
	global spot_selected
	global tool
	global px py
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set scale [$wc.scale.scale get]
	set n [$wc.image.scale get]
	set z [$wc.slice.scale get]
	set px [expr floor( [$c canvasx $x] / $scale ) ]
	set py [expr floor( [$c canvasy $y] / $scale ) ]
	set pyi [yFlip $py]
	set mg_item [micrographItem $n]
	set id [Bmg spot $mg_item select $px $pyi $z]
#	puts "selected ID = $id"
	if { [llength $id] > 0 } {
		spotDraw [lindex $id 0] [lindex $id 1] [lindex $id 2]
		set fom [Bmg spot $mg_item fom [lindex $id 0] [lindex $id 1] [lindex $id 2]]
		.wxtal.selected.h config -text [lindex $id 0]
		.wxtal.selected.k config -text [lindex $id 1]
		.wxtal.selected.l config -text [lindex $id 2]
		.wxtal.selected.fom config -text $fom
	}
	set spot_selected $id
	return $id
}

## @brief Creates or selects a structure factor.
#
# @param	x
# @param	y			Structure factor position in canvas coordinates.

proc spotCreateSelect { x y } {
	if ![winfo exists .wxtal] { return }
	global theimg
	global spot_selected
	global tool
	global px py
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set id [spotSelect $x $y]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	if { [llength $id] < 1 && $tool == "xtal" } {
		setSpotIndices
#		puts "After indices set: $spot_selected"
		if { [llength $spot_selected] < 1 } { return }
		set scale [$wc.scale.scale get]
		set z [$wc.slice.scale get]
		set px [expr floor( [$c canvasx $x] / $scale ) ]
#		set pyi [expr $h - $py - 1]
		set pyi [yFlip $py]
		set id [Bmg spot $mg_item create [lindex $spot_selected 0] [lindex $spot_selected 1] \
			[lindex $spot_selected 2] $px $pyi $z]
#		puts "created ID = $id"
	}
	if { [llength $id] > 0 } {
		spotDraw [lindex $id 0] [lindex $id 1] [lindex $id 2]
	}
	set spot_selected $id
	set i [Bmg spot $mg_item count]
	.wxtal.number.num config -text "$i"
}

## @brief Moves the selected structure factor.
#
# @param	x
# @param	y			SF position in canvas coordinates.

proc spotMove { x y } {
	if ![winfo exists .wxtal] { return }
	global theimg
	global spot_selected
	global px py
	if { [llength $spot_selected] < 1 } { return }
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set scale [$wc.scale.scale get]
	set n [$wc.image.scale get]
	set ix [expr [$c canvasx $x] / $scale ]
	set iy [expr [$c canvasy $y] / $scale ]
	set dx [expr $ix - $px]
	set dy [expr $py - $iy]
	set dz 0
	set px $ix
	set py $iy
#	puts "Moving $dx $dy"
	set mg_item [micrographItem $n]
	Bmg spot $mg_item move [lindex $spot_selected 0] [lindex $spot_selected 1] \
		[lindex $spot_selected 2] $dx $dy $dz
	spotsDrawAll
}

proc spotDisplayLocation { loc } {
	global theimg
	if ![winfo exists .wxtal] { return }
#	set xm [expr [Bimage get $theimg width] / 2]
#	set ym [expr [Bimage get $theimg height] / 2]
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set scale [$wc.scale.scale get]
#	set x [expr ([lindex $loc 0] + $xm) * $scale]
#	set y [expr [lindex $loc 1] + $ym]
	set x [expr [lindex $loc 0] * $scale]
	set y [lindex $loc 1]
	set y [expr [yFlip $y] * $scale ]
#	set z [$wc.slice.scale get]
	set z [lindex $loc 2]
	set loc "$x $y $z"
	return $loc
}

## @brief Draws a structure factor on the canvas.
#
# @param	h				SF index.
# @param	k				SF index.
# @param	l				SF index.

proc spotDraw { h k l } {
	global theimg
	global spot_selected
	global show_spot_indices
	global spot_color spot_select_color
	if ![winfo exists .wxtal] { return }
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set z [$wc.slice.scale get]
	set scale [$wc.scale.scale get]
	set mg_item [micrographItem $n]
	set loc [Bmg spot $mg_item location $h $k $l]
	set loc [spotDisplayLocation $loc]
#	puts $loc
#	set x [expr [lindex $loc 0] * $scale ]
#	set y [expr [yFlip [lindex $loc 1]] * $scale ]
	set x [lindex $loc 0]
	set y [lindex $loc 1]
	set zr [expr abs([lindex $loc 2] - $z) ]
	set radius [expr [.wxtal.spot_radius.entry get] * $scale]
	set xmin [expr $x - $radius ]
	set xmax [expr $x + $radius + $scale - 1 ]
	set ymin [expr $y - $radius ]
	set ymax [expr $y + $radius + $scale - 1 ]
	set color $spot_color
	if { $h == [lindex $spot_selected 0] && $k == [lindex $spot_selected 1]  && $l == [lindex $spot_selected 2] } {
		set color $spot_select_color
	}
	$c create oval $xmin $ymin $xmax $ymax -outline $color -width 1 -tags spot
	if { $show_spot_indices } {
		$c create text $xmax $y -text "$h,$k" -fill $color -anchor w -tags spot
	}
}

## @brief Draws all structure factors on the canvas.
#

proc spotsDrawAll { } {
	global theimg
	global show_spot fom_cutoff
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	$c delete spot
	if { $show_spot < 1 } { return }
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set ids [Bmg spot $mg_item ids $fom_cutoff]
	set nsel 0
	foreach { h k l} $ids {
		spotDraw $h $k $l
		incr nsel 1
	}
	.wxtal.number.num config -text "$nsel"
}

## @brief Deletes the selected structure factor on the canvas.
#
# @param	x
# @param	y			SF particle position in window coordinates.

proc spotDelete { x y } {
	global theimg
	global spot_selected
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set id [spotSelect $x $y]
	puts "Deleting $id"
	if { [llength $id] > 0 } {
		set mg_item [micrographItem $n]
		Bmg spot $mg_item delete [lindex $id 0] [lindex $id 1] [lindex $id 2]
		set spot_selected []
	}
	spotsDrawAll
}

## @brief Deletes all structure factors.
#

proc spotsDeleteAll { } {
	global theimg
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	Bmg spot $mg_item delete all
	$c delete spot
	.wxtal.number.num config -text "0"
}

## @brief Generates structure factors from unit cell vectors.
#

proc spotGenerate { } {
	global spot_resolution
	set w .wxtalres
	catch {destroy $w}
	toplevel $w
	wm title $w "Set the resolution limit"
	wm iconname $w "Resolution"
	tkwait visibility $w
	grab $w
	frame $w.buttons
	button $w.buttons.ok -text OK -command "spotGenerate_and_DestroyWindow $w"
	button $w.buttons.cancel -text Cancel -command "destroy $w"
	frame $w.res
	label $w.res.tag -text "Resolution limit" -anchor e
	entry $w.res.entry -width 10 -validate key -vcmd { string is double %P }
	$w.res.entry insert 0 "$spot_resolution"
	pack $w.res.tag $w.res.entry -side left -padx 5 -pady 5
	pack $w.res -side top -fill x -pady 2
	pack $w.buttons -side bottom -fill x -pady 2m
	pack $w.buttons.ok $w.buttons.cancel -side left -expand 1
	bind $w <Return> "spotGenerate_and_DestroyWindow $w"
	tkwait window $w
}

## @brief Changes the data type of an image.
#
# @param	w 			the window to destroy.

proc spotGenerate_and_DestroyWindow { w } {
	global theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set resolution [$w.res.entry get]
#	set pixel_size [$wc.pixel_size.entry get]
#	set nx [Bimage get $theimg width]
#	set hkl_lim [expr int($nx*$pixel_size/$resolution)]
#	puts "HKL limit = $hkl_lim"
	set mg_item [micrographItem $n]
	Bmg spot $mg_item generate $resolution
	destroy $w
	updateSpotParam
	spotsDrawAll
}

## @brief Masks structure factors.
#

proc spotMask { } {
	global theimg spot_radius
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	Bmg spot $mg_item mask $spot_radius
	Update 1
}

## @brief Updates structure factor parameters from the micrograph parameters in memory.
#

proc updateSpotParam { } {
	global project_item imgtype
	global spot_radius spot_selected
	set w .wxtal
	if ![winfo exists $w] { return }
	if ![Bmg exists] { createMgParam }
	set filename [Bmg get $project_item filename $imgtype]
	wm title $w [list "Xtal:" $filename]
#	puts "Updating spot parameters"
	set pixel_size [Bmg get $project_item pixel_size]
	Bmg set $project_item sf_radius $spot_radius
#	set spot_radius [Bmg get $project_item spot_radius]
	updatePixelsize
	set spot_selected []
	$w.spot_radius.entry delete 0 end
	$w.spot_radius.entry insert 0 $spot_radius
	set unitcell [Bmg get $project_item unitcell]
	$w.unitcell.hx config -text [lindex $unitcell 0]
	$w.unitcell.hy config -text [lindex $unitcell 1]
	$w.unitcell.kx config -text [lindex $unitcell 2]
	$w.unitcell.ky config -text [lindex $unitcell 3]
}

