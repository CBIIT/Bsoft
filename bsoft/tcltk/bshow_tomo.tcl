##
# @file		bshow_tomo.tcl
#
# @brief	Procedures to manage tomographic images
#
# @author	Bernard Heymann
# @date		Created: 20021112
# @date		Modified: 20200325

set show_markers 1
set show_marker_errors 1
set show_marker_labels 1
set marker_radius 10
set markers_selected []
set refop "m"
set hires 10
set lores 1000
set samrat 1
set norm_type "Simple"
set shift_limit 500
set cc_type 0

## @brief Dialog box to manage tomographic information
#

proc Tomography { } {
#	puts "In Tomography"
	global tool project_item theimg imgtype
	global helv12 cour12
	global project_item
	global thickness tilt_axis marker_radius

	if ![Bimage exists $theimg] {
		tk_dialog .dialog "No image in memory!" \
				{Please read in an appropriate image first.} \
				info 0 {OK}
		return
	}
	
	set c [getImageCanvas $theimg]
	if ![Bmg exists] {
		set imgtype "mg"
		createMgParam
		setMgParam
		Bmg set $project_item axis $tilt_axis
		Bmg set $project_item marker_radius $marker_radius
#		puts "Setting tilt axis in the project"
	} else {
		set axis [expr [Bmg get $project_item axis] ]
		if { $axis > -1000 } {
			set tilt_axis $axis
#			puts "Getting tilt axis from the project"
		}
		set r [Bmg get $project_item marker_radius]
		if $r {
			set marker_radius $r
		} else {
#			set marker_radius 10
			Bmg set $project_item marker_radius $marker_radius
		}
	}

	set field [Bmg get $project_item field]
	set nmg [Bmg get $project_item number_of_mg]
	set nrec [Bmg get $project_item number_of_rec]
	
#	puts "Number of micrographs: $nmg"
#	puts "Number of reconstructions: $nrec"
#	puts "Field: $field"
#	puts "Marker radius: $marker_radius"
#	puts "Tilt axis: $tilt_axis"
	
	set w .wtomo
	if ![winfo exists $w] {
		toplevel $w
		menu $w.menuBar -tearoff 0
		$w.menuBar add cascade -menu $w.menuBar.tomo -label "Tomography" -underline 0
		menu $w.menuBar.tomo -tearoff 0
		$w configure -menu $w.menuBar
		$w.menuBar.tomo add command -label "Help" -underline 0 \
				-command { showURL "file://$Bsoft/doc/bshow/bshow_tomo.html" }
#				-command { showHelp bshow_tomo.hlp }
		$w.menuBar.tomo add command -label "Set tilt angles" \
				-command { Tilts } -underline 0
		$w.menuBar.tomo add command -label "Read tilt angles file" \
				-command { readTlt } -underline 0
#		$w.menuBar.tomo add command -label "Find markers in all images" \
#				-command { markerFind -1 } -underline 0
		$w.menuBar.tomo add command -label "Renumber markers in current image" \
				-command { markerRenumber } -underline 0
		$w.menuBar.tomo add command -label "Center markers in current image" \
				-command { markerCenter } -underline 0
		$w.menuBar.tomo add command -label "Extract markers" \
				-command { markerExtract } -underline 0
		$w.menuBar.tomo add command -label "Show marker table" \
				-command { markerTable } -underline 0
		$w.menuBar.tomo add command -label "Delete selected markers" \
				-command { markerDeleteSelected } -underline 0
		$w.menuBar.tomo add command -label "Delete selected markers in all images" \
				-command { markerDeleteSelectedAll } -underline 0
		$w.menuBar.tomo add command -label "Delete all markers" \
				-command { markerDeleteAll } -underline 0
		$w.menuBar.tomo add command -label "Close" \
				-command "destroy $w" -underline 0
		
		$w.menuBar add cascade -menu $w.menuBar.tomoflow -label "Workflow" -underline 0
		menu $w.menuBar.tomoflow -tearoff 0
		$w configure -menu $w.menuBar
		$w.menuBar.tomoflow add command -label "Align and sum frames" \
				-command { alignSumFrames } -underline 0
		$w.menuBar.tomoflow add command -label "Normalize images" \
				-command { normalizeImages } -underline 0
		$w.menuBar.tomoflow add command -label "Align micrographs" \
				-command { alignMicrographs } -underline 0
		$w.menuBar.tomoflow add command -label "Find markers in current image" \
				-command { markerFind [.image.scale get] } -underline 0
		$w.menuBar.tomoflow add command -label "Generate markers from seed" \
				-command { markerGenerate } -underline 0
		$w.menuBar.tomoflow add command -label "Find tilt axis" \
				-command { findTiltAxis } -underline 0
		$w.menuBar.tomoflow add command -label "Transfer markers to second series" \
				-command { markerTransfer } -underline 0
		$w.menuBar.tomoflow add command -label "Track markers" \
				-command { markerTrack } -underline 0
		$w.menuBar.tomoflow add command -label "Refine alignment" \
				-command { markerRefine } -underline 0
		$w.menuBar.tomoflow add command -label "Calculate and fit power spectra" \
				-command { tomoPowerSpectra } -underline 0
		$w.menuBar.tomoflow add command -label "Estimate resolution" \
				-command { tomoResolution } -underline 0
		$w.menuBar.tomoflow add command -label "Reconstruct" \
				-command { tomoReconstruct } -underline 0
		$w.menuBar.tomoflow add command -label "Denoise" \
				-command { tomoDenoise } -underline 0
		
		frame $w.switches
		checkbutton $w.switches.markers -text "Show markers" -variable show_markers \
				-command { markerDrawAll }
		checkbutton $w.switches.errors -text "Show errors" -variable show_marker_errors \
				-command { markerDrawAll }
		checkbutton $w.switches.labels -text "Show labels" -variable show_marker_labels \
				-command { markerDrawAll }
		pack $w.switches.markers $w.switches.errors $w.switches.labels -side left -ipadx 2 -ipady 5 -pady 5 -padx 1m
		
		setupEntry $w.thick "Thickness       " double $thickness "The tomogram thickness estimate set by the user or determined using the effective mean free path."
		label $w.thick.unit -text "Å" -width 3 -font $helv12 -anchor w
#		button $w.thick.update -text "Fit intensities" \
#				-relief raised -command "fitIntensities"
		label $w.thick.update -text "Update" -relief raised
		set m [menu $w.thick.pup]
		$m add command -label "Fit intensities" -command "fitIntensities"
		$m add command -label "Marker range" -command "markerRange"
		pack $w.thick.unit $w.thick.update -side left -padx 1m
		bind $w.thick.update <1> [list tk_popup "$w.thick.pup" %X %Y]
		
		setupEntry $w.tilt_axis "Tilt axis       " double $tilt_axis "The general tilt axis angle starting from the x-axis and proceeding counter-clockwise"
		checkbutton $w.tilt_axis.show -text "Show" -variable show_tilt_axis \
				-command { drawTiltAxis }
		bind $w.tilt_axis.e <Return> { setTiltAxis }
		pack $w.tilt_axis.show -side left -padx 1m

		setupEntry $w.marker_radius "Marker radius   " double $marker_radius "The typical radius of fidicial markers in pixels"
		
		frame $w.marker
		label $w.marker.tag -text "Markers" -width 20 -anchor w
		label $w.marker.num -width 10 -text "0" -relief sunken -bd 1  -anchor w
		button $w.marker.tabbut -text "Marker table" \
				-relief raised -command "markerTable"
		pack $w.marker.tag $w.marker.num $w.marker.tabbut -side left -padx 1m

		frame $w.selected
		label $w.selected.numtag -text "Selected marker" -width 20 -anchor w
		label $w.selected.num -width 4 -text "" -relief sunken -bd 1 -anchor w
		label $w.selected.restag -text "Residual" -width 10 -anchor e
		label $w.selected.res -width 6 -text "" -relief sunken -bd 1 -anchor w
		label $w.selected.fomtag -text "FOM" -width 5 -anchor e
		label $w.selected.fom -width 6 -text "" -relief sunken -bd 1 -anchor w
		button $w.selected.refbut -text "Refine" \
				-relief raised -command "refineMarker"
		pack $w.selected.numtag $w.selected.num $w.selected.restag \
				$w.selected.res $w.selected.fomtag $w.selected.fom \
				$w.selected.refbut -side left -padx 1m

		frame $w.resid
		label $w.resid.tag -text "Residual" -width 20 -anchor w
		label $w.resid.val -width 10 -text "" -relief sunken -bd 1 -anchor w
		button $w.resid.update -text "Update" \
				-relief raised -command "updateTomoTable"
		pack $w.resid.tag $w.resid.val $w.resid.update -side left -padx 1m

		frame $w.fom
		label $w.fom.label -text "FOM cutoff" -width 20 -anchor w
		scale $w.fom.scale -orient horizontal -length 300 -from 0 -to 1 \
			-command { setFOMcutoff } -resolution 0.001
		pack $w.fom.label $w.fom.scale -side left -padx 1m

		frame $w.field
		label $w.field.tag -text "Field" -width 20 -anchor w
		label $w.field.str -width 50 -text $field -relief sunken -bd 1 -anchor w
		pack $w.field.tag $w.field.str -side left -padx 1m

		setupMicrographButtons $w.neighbor "micrograph"

		label $w.header -font $helv12 -anchor w \
			-text "   Mg          Tilt           Axis         Level       OriginX    OriginY    ScaleX  ScaleY    Residual"

		frame $w.frame
		text $w.table -relief sunken -bd 2 -width 70\
			-yscrollcommand ".wtomo.yscroll set" \
			-setgrid 1 -highlightthickness 0 -pady 2 -padx 3 -font $cour12
		scrollbar $w.yscroll -command ".wtomo.table yview" \
			-highlightthickness 0 -orient vertical
		grid $w.table -in $w.frame -padx 1 -pady 1 \
			-row 0 -column 0 -rowspan 1 -columnspan 1 -sticky news
		grid $w.yscroll -in $w.frame -padx 1 -pady 1 \
			-row 0 -column 1 -rowspan 1 -columnspan 1 -sticky news
		grid rowconfig    $w.frame 0 -weight 1 -minsize 0
		grid columnconfig $w.frame 0 -weight 1 -minsize 0

		pack $w.switches $w.thick $w.tilt_axis $w.marker_radius $w.marker $w.selected \
				$w.resid $w.fom $w.field $w.neighbor -side top -pady 2 -padx 5 -anchor w
		pack $w.header -side top -fill x -padx 5
		pack $w.frame -expand yes -fill both -padx 5 -pady 5
		
#		puts "Tomo window setup done"

		updateTomoTable
		updateMarkerParam
#		puts [getPixelSizeEntry]
		markerDrawAll
	
	} else {
		wm deiconify $w
		raise $w
    }
	
    wm title $w "Tomography"
    wm iconname $w "Tomography"
	if { $tool == "marker_group" } {
		set n [.image.scale get]
		set mg_item [micrographItem $n]
		set markers_selected [Bmg marker $mg_item selectrect -10000 -10000 10000 10000]
	}

	.menuBar.window entryconfigure "Tomography" -state normal
	bind $c <1> "markerCreateSelect %x %y"
	bind $c <2> "markerCreateSelect %x %y"
	bind $c <Shift-1> "markerDelete %x %y"
	bind $c <Shift-2> "markerDelete %x %y"
	bind $c <B1-Motion> "markerMove %x %y"
	bind $c <B2-Motion> "markerMove %x %y"
	bind $c <ButtonRelease-1> "markerSetSelection %x %y"
	bind $c <ButtonRelease-2> "markerSetSelection %x %y"
	bind $w <Return> "Update 0"
	bind $w <Control-w> { destroy .wtomo }
	bind $w <Command-w> { destroy .wtomo }
	bind . <Key-a> "markerAcceptModelLocation" 
#	bind . <Key-i> "invertSelection $c" 
	bind $c <Control-1> "markerSelectAppend %x %y"
	bind $c <Command-1> "markerSelectAppend %x %y"
}

## @brief Sets the tilt axis angle
#

proc setTiltAxis { } {
	global tilt_axis project_item
	set tilt_axis [.wtomo.tilt_axis.e get]
	if [Bmg exists] {
		Bmg set $project_item axis $tilt_axis
		updateTomoTable
	}
}

## @brief Dialog box to set the tilts for a tomographic series

proc Tilts { } {
	global helv12
	global tilt_axis project_item
	set tilt_axis [.wtomo.tilt_axis.e get]
	set tilt_step 1
	set tilt_start -60
	if [Bmg exists] {
		set mg_item [micrographItem 0]
		set tilt_axis [expr [Bmg get $mg_item axis]]
		set tilt_start [expr [Bmg get $mg_item tilt]]
		set mg_item [micrographItem 1]
		set tilt_step [expr [Bmg get $mg_item tilt] - $tilt_start]
	}
	set w .wtilt
	if ![winfo exists $w] {
		toplevel $w
		
		setupEntry $w.tilt_step "Tilt step     " double $tilt_step \
			"The angular increment between micrographs"

		setupEntry $w.tilt_start "Starting angle" double $tilt_start \
			"The tilt angle of the first micrograph"

		pack $w.tilt_start $w.tilt_step -side top -fill x -pady 2 -padx 2
		
		frame $w.buttons
		button $w.buttons.ok -text OK -command "setTilts .wtilt"
		button $w.buttons.cancel -text Cancel -command "destroy .wtilt"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.ok $w.buttons.cancel -side left -expand 1
	
		bind $w <Return> "setTilts .wtilt"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Tilts"
    wm iconname $w "Tilts"
}

## @brief Generates the tilt list for a tomographic series.
#
# @param	w 			the window to destroy.

proc setTilts { w } {
	global tilt_axis project_item
	set step [expr [$w.tilt_step.e get]]
	set start [expr [$w.tilt_start.e get]]
	destroy $w
	set nimg [Bmg get $project_item number_of_mg]
#	puts "Number of micrographs = $nimg"
	setTiltAxis
	for {set i 0} {$i < $nimg} {incr i 1} {
		set tilt [expr $start + $i * $step]
		set mg_item [micrographItem $i]
		Bmg set $mg_item tilt $tilt
	}
	Bmg update_matrices
	updateTomoTable
}

## @brief Reads tilt angles from a rawtlt file
#

proc readTlt { } {
	global project_item
    set types {
		{"tlt files"		{.tlt .TLT}	}
		{"rawtlt files"		{.rawtlt .RAWTLT}	}
	}
	set rawtltfile [tk_getOpenFile -filetypes $types]
	if { [string length $rawtltfile] < 1 } { return }
	set rawtltID [open $rawtltfile r]
	set nimg [Bmg get $project_item number_of_mg]
	set i 0
	foreach line [split [read $rawtltID] \n] {
		if { [scan $line "%g" tilt ] > 0 } {
			if { $i < $nimg } {
				Bmg set $project_item tilt $tilt $i
			}
			incr i 1
		}
	}
	close $rawtltID
	set tilt_axis [.wtomo.tilt_axis.e get]
	Bmg update_matrices
	updateTomoTable
}

## @brief Draws the tilt axis on the canvas
#

proc drawTiltAxis { } {
	global theimg PI
	global show_tilt_axis project_item
	if ![Bmg exists] { return }
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	$c delete tilt_axis
	if { [$wc.slice.scale get] > 1 } { set show_tilt_axis 0 }
	if { $show_tilt_axis < 1 } { return }
	set img_num [$wc.image.scale get]
	set mg_item [micrographItem $img_num]
	set tax [expr [Bmg get $mg_item axis] * $PI/180.0]
	set ang [expr [Bmg get $mg_item tilt] * $PI/180.0]
#	puts "$mg_item $tax $ang"
	set img_num [$wc.image.scale get]
	set scale [$wc.scale.scale get]
	set height [Bimage get $theimg height]
	set radius [expr $height * $scale / 2]
	set origin [Bmg get $mg_item origin]
	set ox [expr [lindex $origin 0] * $scale]
	set oy [expr ($height - [lindex $origin 1] - 1) * $scale]
	set ct [expr cos($tax)]
	set st [expr sin($tax)]
	set dx [expr $radius*$ct]
	set dy [expr $radius*$st]
	set xs [expr $ox - $dx]
	set ys [expr $oy + $dy]
	set xe [expr $ox + $dx]
	set ye [expr $oy - $dy]
	$c create line $xs $ys $xe $ye -fill LightBlue -dash - -width 1.5 -tags tilt_axis
#	puts "$tax $xs $ys $xe $ye"
	# Note that this line points to the tomogram z-axis
	set vlen [expr 0.1 * $height * $scale * sin($ang)]
	set xe [expr $ox + $vlen*$st]
	set ye [expr $oy + $vlen*$ct]
	$c create line $ox $oy $xe $ye -fill LightBlue -width 1.5 -tags tilt_axis
#	puts "$vlen $ox $oy $xe $ye"
}

## @brief Selects a marker
#
# @param	x
# @param	y			Marker position in canvas coordinates.

proc markerSelect { x y } {
	if ![winfo exists .wtomo] { return }
	global theimg project_item
	global marker_radius markers_selected
	global px py
	global mx my
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set marker_radius [.wtomo.marker_radius.e get]
	Bmg set $project_item marker_radius $marker_radius
#	puts $marker_radius
	set n [$wc.image.scale get]
	set z [$wc.slice.scale get]
	set scale [$wc.scale.scale get]
	set px [expr floor( [$c canvasx $x] / $scale ) ]
	set py [expr floor( [$c canvasy $y] / $scale ) ]
	set mx [expr int($px)]
	set my [expr int($py)]
	set myi [yFlip $my]
	set mg_item [micrographItem $n]
	set id [Bmg marker $mg_item select $mx $myi $z]
	set markers_selected $id
	set res [Bmg marker $mg_item residual $id]
	set fom [Bmg marker $mg_item fom $id]
	.wtomo.selected.num config -text "$id"
	.wtomo.selected.res config -text "$res"
	.wtomo.selected.fom config -text "$fom"
	return $id
}

## @brief Appends marker at a location to list of selected markers.
#
# @param	x
# @param	y			Cursor position in window coordinates.

proc markerSelectAppend { x y } {
	global markers_selected
	set id [markerSelect $x $y]
	if { $id < 1 } { return }
	if { [lsearch $markers_selected $id] < 0 } {
		lappend markers_selected $id
	}
}

## @brief The marker selection is reset and set to within a rectangle.
#
# @param	x
# @param	y			Marker position in canvas coordinates.

proc markerSetSelection { x y } {
	global theimg project_item
	global markers_selected
	global tool
	global px py
#	set markers_selected []
	if { $tool != "marker_group" } { return }
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set scale [$wc.scale.scale get]
	set n [$wc.image.scale get]
	set x1 $px
	set y1 [yFlip $py]
	set x2 [expr [$c canvasx $x] / $scale ]
	set y2 [yFlip [expr [$c canvasx $y] / $scale ]]
#	puts "$x1 $y1 $x2 $y2"
	set mg_item [micrographItem $n]
	set markers_selected [Bmg marker $mg_item selectrect $x1 $y1 $x2 $y2]
#	puts "selected IDS = $markers_selected"
	markerDrawAll
}

## @brief Selects or creates a marker
#
# @param	x
# @param	y			Marker position in canvas coordinates.

proc markerCreateSelect { x y } {
	if ![winfo exists .wtomo] { return }
	global theimg project_item
	global marker_radius markers_selected
	global tool
	global px py
	global mx my
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set marker_radius [.wtomo.marker_radius.e get]
	Bmg set $project_item marker_radius $marker_radius
#	puts $marker_radius
	set n [$wc.image.scale get]
	set z [$wc.slice.scale get]
	set scale [$wc.scale.scale get]
	set px [expr floor( [$c canvasx $x] / $scale ) ]
	set py [expr floor( [$c canvasy $y] / $scale ) ]
	set mx [expr int($px)]
	set my [expr int($py)]
	set myi [yFlip $my]
	set mg_item [micrographItem $n]
	set id [Bmg marker $mg_item select $mx $myi $z]
#	puts "selected ID = $id"
	if { $id < 1 } {
		$c delete selection
		if { $tool == "marker" } {
			set id [Bmg marker $mg_item create $mx $myi $z]
#			puts "created ID = $id"
		}
		if { $tool == "marker_group" } {
			set markers_selected []
		}
	}
	if { $tool == "marker_group" } {
		if { $id > 0 && [llength $markers_selected] == 0 } {
			set markers_selected [Bmg marker $mg_item ids 0]
		}
		return
	}
	if { $id > 0 } {
		set markers_selected $id
#		Bmg marker $mg_item snap $id
		set res [Bmg marker $mg_item residual $id]
		set fom [Bmg marker $mg_item fom $id]
		.wtomo.selected.num config -text "$id"
		.wtomo.selected.res config -text "$res"
		.wtomo.selected.fom config -text "$fom"
	}
#	puts "Selected IDs = $markers_selected"
	markerDrawAll
#	Update 0
}

## @brief Generates markers from a seed.
#

proc markerGenerate { } {
	global theimg
	global show_tilt_axis project_item
	set show_tilt_axis 1
	Bmg marker $project_item generate
	updateOrigins
	Update 0
}

## @brief Extracts markers from a micrograph.
#

proc markerExtract { } {
	global filetypes
	global filename project_item theimg
	if ![winfo exists .wtomo] { return }
	set wc [getControlWindow $theimg]
	set xname [file root [file tail $filename]]
	set xdir [file dirname $filename]
	append xname "_marker.mrc"
	set xname [tk_getSaveFile -filetypes $filetypes \
		-initialfile $xname  -initialdir $xdir -defaultextension .mrc]
	if { [string length $xname] < 1 } { return }
	setMgParam
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	Bmg marker $mg_item extract $xname
}

## @brief Moves the selected markers
#
# @param	x
# @param	y			Mouse position in canvas coordinates.

proc markerMove { x y } {
	if ![winfo exists .wtomo] { return }
	global theimg
	global markers_selected project_item
	global tool
	global px py
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set scale [$wc.scale.scale get]
	set cx [$c canvasx $x]
	set cy [$c canvasy $y]
	if { $tool == "marker_group" && [llength $markers_selected] == 0 } {
		$c delete selection
		set cpx [expr int($px * $scale)]
		set cpy [expr int($py * $scale) ]
		$c create rectangle $cpx $cpy $cx $cy -tags selection -outline red
		return
	}
	set n [$wc.image.scale get]
	set ix [expr $cx / $scale ]
	set iy [expr $cy / $scale ]
	set dx [expr $ix - $px]
	set dy [expr $py - $iy]
	set px $ix
	set py $iy
	set mg_item [micrographItem $n]
	foreach id $markers_selected {
		Bmg marker $mg_item move $id $dx $dy
	}
	markerDrawAll
}

## @brief Draws a marker with an error line and label.
#
# @param	id			Marker identifier.

proc markerDraw { id } {
	global theimg project_item
	global show_marker_errors show_marker_labels
	global marker_color marker2_color marker_select_color error_color error_select_color
	global marker_radius markers_selected
	if ![winfo exists .wtomo] { return }
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set z [$wc.slice.scale get]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set locerr [Bmg marker $mg_item location $id]
#	puts $locerr
	set selflag [Bmg marker $mg_item getflag $id]
#	puts $selflag
	set zr [expr abs([lindex $locerr 2] - $z)]
	set marker_radius [.wtomo.marker_radius.e get]
	if { $zr > $marker_radius } { return }
	set zr2 1
	if { $zr > 0.1 && $marker_radius > 0 } {
		set zr2 [expr sqrt(1 - $zr*$zr*1.0/($marker_radius*$marker_radius))]
	}	
	set scale [$wc.scale.scale get]
	set x [expr [lindex $locerr 0] * $scale ]
	set y [expr [yFlip [lindex $locerr 1]] * $scale ]
	set radius [expr $marker_radius * $scale * $zr2]
	set ex [expr $x + [lindex $locerr 3] * $scale ]
	set ey [expr $y - [lindex $locerr 4] * $scale ]
#	puts "$id $locerr"
	set xmin [expr $x - $radius ]
	set xmax [expr $x + $radius ]
	set ymin [expr $y - $radius ]
	set ymax [expr $y + $radius ]
	if { [lsearch $markers_selected $id] >= 0 } {
		set markcolorshow $marker_select_color;
		set errcolorshow $error_select_color;
	} elseif { $selflag == 2 } {
		set markcolorshow $marker2_color;
		set errcolorshow $error_color;
	} else {
		set markcolorshow $marker_color;
		set errcolorshow $error_color;
	}
	set idtag [format "m%d" $id]
	$c create oval $xmin $ymin $xmax $ymax -outline $markcolorshow -width 1 -tags [list marker $idtag]
	if { $show_marker_errors } {
		$c create line $x $y $ex $ey -fill $errcolorshow -width 2 -tags [list marker $idtag]
	}
	if { $show_marker_labels } {
		$c create text $xmax $y -text $id -fill $markcolorshow -anchor w -tags [list marker $idtag]
	}
}

## @brief Redraws all markers on the canvas
#

proc markerDrawAll {} {
	global theimg project_item
	global fom_cutoff
	global show_markers
	if ![winfo exists .wtomo] { return }
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	$c delete marker
	if { $show_markers < 1 } { return }
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
#	puts "markerDrawAll $mg_item"
	set ids [Bmg marker $mg_item ids $fom_cutoff]
#	puts "fom_cutoff = $fom_cutoff"
#	puts "IDs = $ids"
	foreach id $ids {
#		puts "drawing $id"
		markerDraw $id
	}
	set nsel [Bmg marker $project_item count]
	.wtomo.marker.num config -text "$nsel"
}

## @brief Deletes the selected marker.
#
# @param	x
# @param	y			Marker position in canvas coordinates.

proc markerDelete { x y } {
	global project_item theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set id [markerSelect $x $y]
	set mg_item [micrographItem $n]
	if { $id > 0 } {
		Bmg marker $mg_item delete $id
	}
	markerDrawAll
}

## @brief Deletes the selected marker.
#
# @param	id			Marker identifier.

proc markerDeleteID { id } {
	global project_item theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	if { $id > 0 } {
		Bmg marker $mg_item delete $id
	}
	markerDrawAll
}

## @brief Deletes all the selected markers.
#

proc markerDeleteSelected { } {
	global markers_selected project_item theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
#	puts "IDs = $markers_selected"
	set mg_item [micrographItem $n]
	foreach id $markers_selected {
		Bmg marker $mg_item delete $id
	}
	markerDrawAll
}

## @brief Deletes all the selected markers in all images.
#

proc markerDeleteSelectedAll { } {
	global markers_selected project_item
#	puts "IDs = $markers_selected"
	foreach id $markers_selected {
		Bmg marker all delete $id
	}
	markerDrawAll
}

## @brief Deletes all markers.
#

proc markerDeleteAll { } {
	global theimg project_item
	set c [getImageCanvas $theimg]
	Bmg marker all delete all
	$c delete marker
	.wtomo.marker.num config -text "0"
}

## @brief Renumbers markers in the current micrograph.
#

proc markerRenumber { } {
	global project_item theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	Bmg marker $mg_item renumber
	markerDrawAll
}

## @brief Centers markers in the current micrograph.
#

proc markerCenter { } {
	global project_item theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	Bmg marker $mg_item center
	markerDrawAll
}


## @brief Finds markers in a file.
#
# @param	n				Image number.

proc markerFind { n } {
	global theimg
	global marker_radius
	global fom_cutoff project_item
	if ![winfo exists .wtomo] { return }
	set marker_radius [.wtomo.marker_radius.e get]
	set mg_item [micrographItem $n]
	set nmark [Bmg marker $mg_item find $marker_radius]
	if { $nmark < 1 } {
		tk_messageBox -icon error -type ok -title "Error" -message \
			"No markers found!"
		return
	}
	
	set max_fom [Bmg marker $mg_item fom_maximum]
	set fom_cutoff [expr $max_fom / 2]
	set fom_inc [expr $max_fom/100]
	set h [Bimage get $theimg height]
	set maxedge [expr int($h / 4)]
	markerDrawAll
	
	if { $nmark } {
		set w .wtomfom
		catch {destroy $w}
		toplevel $w
		wm title $w "Set the FOM cutoff"
		wm iconname $w "FOMcut"
		tkwait visibility $w
		grab $w
		frame $w.fom
		label $w.fom.label -text "Choose a FOM cutoff and close this window before doing anything else" -anchor w
		scale $w.fom.scale -orient horizontal -length 500 -from 0 -to $max_fom \
			-command { setFOMcutoff } -resolution $fom_inc
		label $w.fom.edgelabel -text "Choose an edge width in pixels to exclude markers" -anchor w
		scale $w.fom.edge -orient horizontal -length 500 -from 0 -to $maxedge \
			-command { setEdge } -resolution 1
		button $w.fom.done -text "Done" -relief raised -command { destroy .wtomfom }
		pack $w.fom.label -side top -fill x
		pack $w.fom.scale -side top -fill x
		pack $w.fom.edgelabel -side top -fill x
		pack $w.fom.edge -side top -fill x
		pack $w.fom.done -side right
		$w.fom.scale set $fom_cutoff
		pack $w.fom -side top -fill x -pady 2
		tkwait window $w
	}
	
	Bmg marker $mg_item delete fom $fom_cutoff
	Bmg marker $mg_item renumber
	markerDrawAll
	set fom_cutoff 0.0
}

proc setEdge { edge } {
	global project_item theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	Bmg marker $mg_item edge $edge
	Update 0
}

## @brief Centers the selected marker.
#
# @param	x
# @param	y			Marker position in canvas coordinates.

proc markerSnap { x y } {
	global project_item theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set id [markerSelect $x $y]
	set mg_item [micrographItem $n]
	if { $id > 0 } {
		Bmg marker $mg_item snap $id
	}
	markerDrawAll
}

## @brief Updates the residual calculation which also sets all marker residuals.
#

proc updateResidual { } {
	global project_item
	if ![winfo exists .wtomo] { return }
#	puts "Updating residuals"
	set resid [Bmg marker $project_item residual]
	.wtomo.resid.val config -text $resid
}

## @brief Updates the tomographic table from the micrograph parameters in memory.
#

proc updateTomoTable { } {
	global tilt_axis project_item
	if ![winfo exists .wtomo] { return }
	if ![Bmg exists] { return }
	set act [Bmg get all active]
#	puts "Updating tomo table with active = $act"
	set nmg [Bmg get $project_item number_of_mg]
	if { $nmg < 1 } { return }
#	puts "Number of micrographs: $nmg"
	Bmg set all active 0
	updateResidual
	set wt 0.0
	set tax 0.0
	set table ""
	for { set i 0 } { $i < $nmg } { incr i 1 } {
		set mg_item [micrographItem $i]
#		puts "Project: $project_item  mg: $mg_item"
		set axis [expr [Bmg get $mg_item axis]]
		set tilt [expr [Bmg get $mg_item tilt]]
		set level [expr [Bmg get $mg_item level]]
		set origin [Bmg get $mg_item origin]
		set scale [Bmg get $mg_item scale]
		set ox [lindex $origin 0]
		set oy [lindex $origin 1]
		set sx [lindex $scale 0]
		set sy [lindex $scale 1]
		set fom [Bmg get $mg_item fom]
#		puts "$i $tilt $axis $level $ox $oy $sx $sy $fom"
		append table [format " %3d  %6.2f %6.2f %6.2f  %6.1f %6.1f  %5.3f %5.3f  %6.3f\n" $i $tilt $axis $level $ox $oy $sx $sy $fom]
		if { [expr abs($tilt)] > 2 } {
			set w1 [expr 1 - exp(-$tilt*$tilt/0.5)]
			set tax [expr $tax + $axis]
			set wt [expr $wt + $w1]
		}
	}
	Bmg set all active $act
	if { $wt } { set tilt_axis [expr $tax / $wt] }
	.wtomo.tilt_axis.e delete 0 end
	.wtomo.tilt_axis.e insert 0 $tilt_axis
	.wtomo.field.str config -text [Bmg get $project_item field]
	.wtomo.table delete 1.0 end
	.wtomo.table insert 1.0 $table
	drawTiltAxis
}

## @brief Updates marker parameters from the micrograph parameters in memory.
#

proc updateMarkerParam { } {
	global marker_radius project_item
	if ![Bmg exists] { return }
#	puts "Updating marker parameters"
	updatePixelsize
	set r [Bmg get $project_item marker_radius]
	if { $r } {
		set marker_radius $r
	}
	Bmg set $project_item marker_radius $r
	set w .wtomo
	if [winfo exists $w] {
#		updateItemID $w.itemid
		$w.marker_radius.e delete 0 end
		$w.marker_radius.e insert 0 $marker_radius
	}
}

## @brief Sets the selection flag for the selected markers.
#

proc markerSetSelectFlag { } {
	global markers_selected project_item theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	foreach id $markers_selected {
		Bmg marker $mg_item setflag $id
	}
	updateMarkerTable
}

## @brief Clears the selection flag for the selected markers.
#

proc markerClearSelectFlag { } {
	global markers_selected project_item theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	foreach id $markers_selected {
		Bmg marker $mg_item clearflag $id
	}
	updateMarkerTable
}

## @brief Sets up the window for the marker table and updates it.
#

proc markerTable { } {
	if ![winfo exists .wtomo] { return }
	if ![Bmg exists] { return }
	global cour12 helv12
	global refop hires lores
	set w .wmark
	if ![winfo exists $w] {
		toplevel $w
		
		frame $w.frame
		text $w.table -relief sunken -bd 2 -width 40\
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

		label $w.header -font $helv12 -anchor w \
			-text "     Mg       Marker     Residual        FOM     Select"

		pack $w.header -side top -fill x -padx 5
		pack $w.frame -expand yes -fill both -padx 5 -pady 5
		
		frame $w.butres
		label $w.butres.label -font $helv12 -anchor w -text Residual -width 8
		button $w.butres.min -text "|<" -command "markerNextResidual $w min"
		button $w.butres.lower -text "<" -command "markerNextResidual $w lower"
		button $w.butres.higher -text ">" -command "markerNextResidual $w higher"
		button $w.butres.max -text ">|" -command "markerNextResidual $w max"
		entry $w.butres.entry -width 10 -validate key -vcmd [list checkType double %P]

		frame $w.butfom
		label $w.butfom.label -font $helv12 -anchor w -text FOM -width 8
		button $w.butfom.min -text "|<" -command "markerNextFOM $w min"
		button $w.butfom.lower -text "<" -command "markerNextFOM $w lower"
		button $w.butfom.higher -text ">" -command "markerNextFOM $w higher"
		button $w.butfom.max -text ">|" -command "markerNextFOM $w max"
		entry $w.butfom.entry -width 10 -validate key -vcmd [list checkType double %P]

		frame $w.flag
		button $w.flag.clear -text "Clear selection" -command "markerClearSelectFlag"
		button $w.flag.set -text "Set selection" -command "markerSetSelectFlag"
		
		frame $w.buttons
		button $w.buttons.update -text Update -command "updateMarkerTable"
		button $w.buttons.close -text Close -command "destroy $w"
		
		pack $w.butres.label $w.butres.min $w.butres.lower \
			$w.butres.entry $w.butres.higher $w.butres.max -side left -expand 1
		pack $w.butfom.label $w.butfom.min $w.butfom.lower \
			$w.butfom.entry $w.butfom.higher $w.butfom.max -side left -expand 1
#		pack $w.flag.clear $w.flag.set -side left -expand 1
		pack $w.buttons.update $w.buttons.close -side left -expand 1

		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.flag -side bottom -fill x -pady 2m
		pack $w.butfom -side bottom -fill x -pady 2m
		pack $w.butres -side bottom -fill x -pady 2m
	
		$w.table tag configure selmark -background #a0b7ce
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Marker list"
    wm iconname $w "Markers"
	updateMarkerTable
}

## @brief Updates the marker table.
#

proc updateMarkerTable {} {
	global markers_selected project_item
	if ![winfo exists .wmark] { return }
#	puts "Updating tomography table parameters"
	updateResidual
	set m [Bmg marker $project_item list]
	set table ""
	foreach {img_num id res fom sel} $m {
		append table [format " %4d  %4d  %7.3f  %7.3f  %4d\n" $img_num $id $res $fom $sel]
	}
	.wmark.table delete 1.0 end
	.wmark.table insert 1.0 $table
}

## @brief Centers the viewable canvas on the location in the image.
#
# @param	img_num 	image number.
# @param	loc			location list.

proc moveCanvasTo { img_num loc } {
	global theimg
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	$wc.image.scale set $img_num
	set scale [$wc.scale.scale get]
	set wd [Bimage get $theimg width]
	set ht [Bimage get $theimg height]
	set ws [expr [winfo width $c] / (2.0 * $scale)]
	set hs [expr [winfo height $c] / (2.0 * $scale)]
	set fx [expr ([lindex $loc 0] - $ws) / $wd ]
	set fy [expr ($ht - [lindex $loc 1] - 1 - $hs) / $ht ]
	$c xview moveto $fx
	$c yview moveto $fy
}

## @brief Goes to a marker.
#
# @param	w 			dialog window.
# @param	id			marker identifier
# @param	isel		offset in table ?
# @param	img			micrograph index

proc goToMarker { w id isel img } {
	global markers_selected project_item
	if { $id <= 0 } { return }
	set markers_selected $id
	set mg_item [micrographItem $img]
	set locerr [Bmg marker $mg_item location $id]
	moveCanvasTo $img $locerr
	$w.butres.entry delete 0 end
	$w.butres.entry insert 0 [Bmg marker $mg_item residual $id]
	$w.butfom.entry delete 0 end
	$w.butfom.entry insert 0 [Bmg marker $mg_item fom $id]
	$w.table tag add selmark $isel.0 $isel.end
	set isel [expr $isel - 3]
	if { $isel < 0 } { set isel 0 }
	$w.table yview $isel
	markerDrawAll
}

## @brief Goes to the next marker based on its residual.
#
# @param	w 			dialog window.
# @param	dir			direction: min, lower, higher, max

proc markerNextResidual { w dir } {
	global markers_selected
	global theimg project_item
	set wc [getControlWindow $theimg]
	$w.table tag remove selmark 1.0 end
	set sgn 1
	set n -1
	set previd -1
	set prevres -1
	set ressel -1
	if { $dir == "min" || $dir == "higher" } {
		set sgn -1
		set ressel 1e10
	}
	if { $dir == "max" } { set prevres 1e10 }
	if { $dir == "min" || $dir == "max" } {
		set markers_selected []
	}
	if { [llength $markers_selected] > 0 } {
		set n [$wc.image.scale get]
		set mg_item [micrographItem $n]
		set previd [lindex $markers_selected 0]
		set prevres [Bmg marker $mg_item residual $previd]
#		puts "1. $n $previd $prevres"
	}
	set m [Bmg marker $project_item list]
	set i 0
	set isel 0
	set img 0
	set idsel 0
	foreach {img_num id res fom sel} $m {
		incr i
		if { $sel } {
			set dr [expr $sgn * ($ressel - $res)]
			set pdr [expr $sgn * ($res - $prevres)]
			if { !($img_num == $n && $id == $previd) && $dr < 0 && $pdr < 0 } {
				set img $img_num
				set idsel $id
				set ressel $res
				set isel $i
			}
		}
	}
	if { $idsel > 0 } {
		goToMarker $w $idsel $isel $img
	}
#	puts "2. $img $idsel $ressel"
}

## @brief Goes to next marker based on its FOM.
#
# @param	w 			dialog window.
# @param	dir			direction: min, lower, higher, max

proc markerNextFOM { w dir } {
	global markers_selected
	global theimg project_item
	set wc [getControlWindow $theimg]
	$w.table tag remove selmark 1.0 end
	set sgn 1
	set n -1
	set previd -1
	set prevfom -1
	set fomsel -1
	if { $dir == "min" || $dir == "higher" } {
		set sgn -1
		set fomsel 1e10
	}
	if { $dir == "max" } { set prevfom 1e10 }
	if { $dir == "min" || $dir == "max" } {
		set markers_selected []
	}
	if { [llength $markers_selected] > 0 } {
		set n [$wc.image.scale get]
		set mg_item [micrographItem $n]
		set previd [lindex $markers_selected 0]
		set prevfom [Bmg marker $mg_item fom $previd]
#		puts "1. $n $previd $prevfom"
	}
	set m [Bmg marker $project_item list]
	set i 0
	set isel 0
	set img 0
	set idsel 0
	foreach {img_num id res fom sel} $m {
		incr i
		if { $sel } {
			set df [expr $sgn * ($fomsel - $fom)]
			set pdf [expr $sgn * ($fom - $prevfom)]
			if { !($img_num == $n && $id == $previd) && $df < 0 && $pdf < 0 } {
				set img $img_num
				set idsel $id
				set fomsel $fom
				set isel $i
			}
		}
	}
	if { $idsel > 0 } {
		goToMarker $w $idsel $isel $img
	}
#	puts "2. $img $idsel $fomsel"
}

## @brief Accepts the model locations of selected markers.
#

proc markerAcceptModelLocation { } {
	global markers_selected
	global theimg project_item
	if { [llength $markers_selected] < 1 } { return }
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	foreach id $markers_selected {
		Bmg marker $mg_item accept $id
	}
	markerDrawAll
	updateMarkerTable
}

##### Workflow #####

## @brief Aligns a series of frames from dose fractionation.
#

proc alignSumFrames {} {
	global paramfile counts
	checkParamFile
	set w .wframes
	if ![winfo exists $w] {
		toplevel $w
		
		set newfile [newParamFile "_framesum.star"]
		if ![string length $newfile] { destroy $w; return }

		frame $w.head
		label $w.head.lab -text \
			"Frame alignment and summation will be done using the currently saved parameter file:\n$paramfile"
		
		setupNewFileEntry $w.file "Output parameter file" $newfile \
			"New parameter file written after summing the frames"

		radiobutton $w.count -text "Convert to counts" -variable counts -value 1
		
		setupEntry $w.refnum "Reference frame" integer 3 "Number of reference frame to use for the initial progressive alignment"

		labelframe $w.res -text "Resolution limits"
		setupEntry $w.res.hi "High" double 20 "High resolution limit in Å"
		setupEntry $w.res.lo "Low" double 1000 "Low resolution limit in Å"
		pack $w.res.hi $w.res.lo -side top -fill x

		setupEntry $w.bin "Binning        " integer 2 "Binning to speed up alignment"

		setupEntry $w.shift "Shift limit    " double 20 "Limit to cross-correlation shift allowed"

		setupEntry $w.dose "Dose/frame     " double 1 "Dose rate in electrons per frame per Ångstrom squared"

		pack $w.head.lab
		pack $w.head $w.file -side top -pady 2 -padx 2 -fill x
		pack $w.head $w.count -side top -pady 2 -padx 2 -anchor w
		pack $w.head $w.refnum $w.res $w.bin $w.shift $w.dose \
			-side top -pady 2 -padx 2 -fill x

		frame $w.buttons
		button $w.buttons.do -text Do -command "doAlignSumFrames $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doAlignSumFrames $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "AlignFrames"
    wm iconname $w "AlignFrames"
}

## @brief Aligns and averages frames.
#
# @param	w 			dialog window.

proc doAlignSumFrames { w } {
	global env paramfile counts
	set newfile [$w.file.e get]
	set refnum [$w.refnum.e get]
	set hires [$w.res.hi.e get]
	set lores [$w.res.lo.e get]
	set bin [$w.bin.e get]
	set shift [$w.shift.e get]
	set dose [$w.dose.e get]
	setMgParam
	if { ![file exists $paramfile] } {
#		puts "$paramfile not found!"
		tk_dialog .dialog "$paramfile not found!" \
			{Please save an appropriate parameter file first.} \
			info 0 {OK}
		return
	}
	set c ""
	if { $counts } {
		set c "-counts"
	}
	set logfile "[file rootname $newfile].log"
	set cmd "$env(BSOFT)/bin/bseries -verb 1 -frames $c -rate $dose -align $refnum -resol $hires,$lores -shift $shift -bin $bin -write sum -out $newfile $paramfile"
	executeCommand $cmd $logfile
}

## @brief Fit the intensities of a tilt series.
#

proc fitIntensities {} {
	global paramfile project_root
	global thickness emfp adjust_tilt
 	set w .wthick
	if ![winfo exists $w] {
		toplevel $w
		
		panelEMFP $w

		checkbutton $w.flag -text "Apply tilt adjustment" \
			-variable adjust_tilt -relief raised

		set psfile [newParamFile "_int.ps"]
		setupPostscriptEntry $w.file "Postscript file" $psfile \
			"Postscript file with the fitted intensity-tilt plot"

		pack $w.flag $w.file -side top -pady 2 -padx 10

		frame $w.buttons
		button $w.buttons.do -text Fit -command "doFitIntensities $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w.emfp.e <Return> "doFitIntensities $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Fit thickness"
    wm iconname $w "Fit"
	return $thickness
}


## @brief Fits the intensities of a tilt series.
#

proc doFitIntensities { w } {
	global paramfile thickness emfp adjust_tilt
	set emfp [$w.emfp.e get]
	set psfile [$w.file.e get]
	set thickness [Bmg thickness $emfp $adjust_tilt $psfile]
	updateThickness .wtomo
	updateTomoTable
}

## @brief Calcul;ates the thickness from the marker range.
#

proc markerRange { } {
	global project_item thickness
	set thickness [Bmg marker $project_item z_range]
	updateThickness .wtomo
	updateTomoTable
}

## @brief Normalizes a series of images based on simple statistics.
#

proc normalizeImages {} {
	global paramfile
	global norm_type
	checkParamFile
	set w .wnorm
	if ![winfo exists $w] {
		toplevel $w
		
		set newfile [newParamFile "_norm.star"]
		if ![string length $newfile] { destroy $w; return }
		
		frame $w.head
		label $w.head.lab -text \
			"The normalization will be done using the currently saved parameter file:\n$paramfile"
		
		setupNewFileEntry $w.file "Output parameter file" $newfile \
			"New parameter file written after normalizing the micrographs"

		frame $w.type
		label $w.type.l -text "Normalization type"
		tk_optionMenu $w.type.v norm_type "Simple" "Gaussian" "Poisson"
		pack $w.type.l $w.type.v -side left -ipadx 2 -fill x

		pack $w.head.lab
		pack $w.head $w.file $w.type -side top -pady 2 -padx 2 -anchor w

		frame $w.buttons
		button $w.buttons.do -text Do -command "doNormalizeImages $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doNormalizeImages $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Normalize"
    wm iconname $w "Normalize"
}

## @brief Normalizes images.
#
# @param	w 			dialog window.

proc doNormalizeImages { w } {
	global env paramfile
	global norm_type
	set newfile [$w.file.e get]
	setMgParam
	if { ![file exists $paramfile] } {
#		puts "$paramfile not found!"
		tk_dialog .dialog "$paramfile not found!" \
			{Please save an appropriate parameter file first.} \
			info 0 {OK}
		return
	}
	set logfile "[file rootname $newfile].log"
	set cmd "$env(BSOFT)/bin/bnorm -v 1 -data float -type $norm_type -out $newfile $paramfile"
	executeCommand $cmd $logfile
}

## @brief Aligns a tilt series by cross correlation.
#

proc alignMicrographs {} {
	global paramfile
	checkParamFile
	set w .wframes
	if ![winfo exists $w] {
		toplevel $w
		
		set newfile [newParamFile "_aln.star"]
		if ![string length $newfile] { destroy $w; return }

		frame $w.head
		label $w.head.lab -text \
			"Micrograph alignment will be done using the currently saved parameter file:\n$paramfile"
		
		setupNewFileEntry $w.file "Output parameter file" $newfile \
			"New parameter file written after aligning micrographs"

		setupEntry $w.iter "Iterations         " integer 3 "Maximum number of iterations for alignment refinement"

		setupEntry $w.stop "Stopping condition " double 1.0 "Stop iterations when the average change in shift drops below this number of pixels"

		setupEntry $w.adj "Adjcent micrographs" integer 3 "Number of adjacent micrographs to include in each subtomogram reconstruction"

		setupResolutionEntry $w

		pack $w.head.lab
		pack $w.head $w.file $w.iter $w.stop $w.adj $w.res \
			-side top -pady 2 -padx 2 -fill x

		frame $w.buttons
		button $w.buttons.do -text Do -command "doAlignMicrographs $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doAlignMicrographs $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "AlignMicrographs"
    wm iconname $w "AlignMicrographs"
}

## @brief Aligns Micrographs.
#
# @param	w 			dialog window.

proc doAlignMicrographs { w } {
	global env paramfile emfp
	set newfile [$w.file.e get]
	set iter [$w.iter.e get]
	set stop [$w.stop.e get]
	set adj [$w.adj.e get]
	set hires [$w.res.hi.e get]
	set lores [$w.res.lo.e get]
	setMgParam
	if { ![file exists $paramfile] } {
#		puts "$paramfile not found!"
		tk_dialog .dialog "$paramfile not found!" \
			{Please save an appropriate parameter file first.} \
			info 0 {OK}
		return
	}
	set logfile "[file rootname $newfile].log"
	set cmd "$env(BSOFT)/bin/btomaln -verb 1 -align $iter,$stop,$adj -resol $hires,$lores -edge 20,3 -emfp $emfp -out $newfile $paramfile"
	executeCommand $cmd $logfile
}

# findTiltAxis
# Locates the tilt axis.
#

proc findTiltAxis { } {
	global paramfile
	global hires lores
	checkParamFile
	set w .wtax
	if ![winfo exists $w] {
		toplevel $w
		
		set newfile [newParamFile "_findaxis.star"]
		if ![string length $newfile] { destroy $w; return }

		frame $w.head
		label $w.head.lab -text \
			"The search for the tilt axis angle will use the currently saved parameter file:\n$paramfile"
		
		setupNewFileEntry $w.file "Output parameter file" $newfile \
			"New parameter file written after searching for the tilt axis"

		setupEntry $w.tlt "Tilt angle" double 15 "Tilt angle to find tilt axis"
		
		labelframe $w.tax -text "Tilt axis search range"
		setupEntry $w.tax.start "Start" double -180 "Starting angle to test for the tilt axis"
		setupEntry $w.tax.end "End" double 180 "Ending angle to test for the tilt axis"
		setupEntry $w.tax.inc "Increment" double 1 "Angle increment to test for the tilt axis"
		pack $w.tax.start $w.tax.end $w.tax.inc -side top -fill x
		
		setupResolutionEntry $w

		pack $w.head.lab
		pack $w.head $w.file $w.tlt $w.tax $w.res -side top -fill x -pady 2 -padx 2

		frame $w.buttons
		button $w.buttons.do -text Do -command "doFindTiltAxis $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doFindTiltAxis $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "FindAxis"
    wm iconname $w "FindAxis"
}

## @brief Executes the tilt axis search.
#
# @param	w 			dialog window.

proc doFindTiltAxis { w } {
	global env paramfile
	global tilt_axis hires lores
	set newfile [$w.file.e get]
	set tlt [$w.tlt.e get]
	set start [$w.tax.start.e get]
	set end [$w.tax.end.e get]
	set step [$w.tax.inc.e get]
	set hires [$w.res.hi.e get]
	set lores [$w.res.lo.e get]
	setMgParam
	if { ![file exists $paramfile] } {
#		puts "$paramfile not found!"
		tk_dialog .dialog "$paramfile not found!" \
			{Please save an appropriate parameter file first.} \
			info 0 {OK}
		return
	}
	set logfile "[file rootname $newfile].log"
#	Bmg findaxis $axis $start $end $step $hires $lores
	set cmd "$env(BSOFT)/bin/btrack -verb 1 -findaxis $tlt,$step,$start,$end -resol $hires,$lores -out $newfile $paramfile"
	executeCommand $cmd $logfile
}

# markerTransfer
# Transfers markers from one tilt series to a second.
#

proc markerTransfer { } {
	global paramfile
	global hires lores refine
	checkParamFile
	set w .wtrans
	if ![winfo exists $w] {
		toplevel $w

		labelframe $w.ang -text "Rotation angle search range"
		setupEntry $w.ang.start "Start" double -180 "Starting angle to find transfer rotation"
		setupEntry $w.ang.end "End" double 180 "Ending angle to find transfer rotation"
		setupEntry $w.ang.inc "Increment" double 1 "Angle increment to find transfer rotation"
		pack $w.ang.start $w.ang.end $w.ang.inc -side top -fill x
		
		setupResolutionEntry $w

		setupEntry $w.shift "Shift limit" double 100 "Limit on how much the second image may shift to fit the first"

		checkbutton $w.refine -text "Refine markers" \
			-variable refine -anchor w

		pack $w.ang $w.res $w.shift -side top -fill x -pady 2 -padx 2
		pack $w.ang $w.refine -side top -pady 2 -padx 2

		frame $w.buttons
		button $w.buttons.do -text Do -command "doMarkerTransfer $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doMarkerTransfer $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Transfer"
    wm iconname $w "Transfer"
}

## @brief Executes the marker transfer.
#
# @param	w 			dialog window.

proc doMarkerTransfer { w } {
	global env paramfile
	global tilt_axis hires lores refine
	set start [$w.ang.start.e get]
	set end [$w.ang.end.e get]
	set step [$w.ang.inc.e get]
	set hires [$w.res.hi.e get]
	set lores [$w.res.lo.e get]
	set shift [$w.shift.e get]
	setMgParam
	if { ![file exists $paramfile] } {
#		puts "$paramfile not found!"
		tk_dialog .dialog "$paramfile not found!" \
			{Please save an appropriate parameter file first.} \
			info 0 {OK}
		return
	}
	Bmg transfer $start $end $step $hires $lores $shift $refine
	updateOrigins
	updateTomoTable
	markerDrawAll
}

## @brief Tracks the markers through the tilt series.
#

proc markerTrack { } {
	global paramfile
	global hires lores recenter cc_type refine
	global shift_limit thickness
	checkParamFile
	set w .wtrk
	if ![winfo exists $w] {
		toplevel $w
		
		set newfile [newParamFile "_trk.star"]
		if ![string length $newfile] { destroy $w; return }
		
		set logfile [file rootname $paramfile]
		append logfile "_trk.log"
		
		frame $w.head
		label $w.head.lab -text \
			"Tracking will be done on the currently saved parameter file:\n $paramfile"
		pack $w.head.lab

		setupNewFileEntry $w.file "Output parameter file" $newfile \
			"New parameter file written after tracking"

		setupEntry $w.track "Iterations         " integer 5 "Number of iterative tracking cycles"
		setupEntry $w.term "Target residual    " double 1 \
			"The desired residual to terminate before the set number of cycles"

		setupResolutionEntry $w

		setupEntry $w.shift "Shift limit        " double $shift_limit \
			"The maximum number of pixels of shift allowed"
		setupEntry $w.thick "Thickness          " double $thickness \
			"The estimated maximum thickness of the tomogram"
			
		frame $w.cb
		checkbutton $w.cb.cc_type -text "Cross-correlate markers" -anchor e -variable cc_type
		checkbutton $w.cb.recenter -text "Recenter" -anchor e -variable recenter
		checkbutton $w.cb.refine -text "Refine markers" -variable refine
		pack $w.cb.cc_type $w.cb.recenter $w.cb.refine -side left -pady 2 -padx 10

		pack $w.head $w.file $w.track $w.term $w.res $w.shift \
			$w.thick $w.cb -side top -pady 2 -padx 10 -fill x
#		pack $w.cc_type $w.recenter -side top -pady 2 -padx 10

		frame $w.buttons
		button $w.buttons.do -text Do -command "doMarkerTrack $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doMarkerTrack $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Track"
    wm iconname $w "Track"
}

## @brief Executes marker tracking.
#
# @param	w 			dialog window.

proc doMarkerTrack { w } {
	global env paramfile
	global hires lores recenter cc_type refine
	global shift_limit thickness
	set newfile [$w.file.e get]
	set iter [$w.track.e get]
	set term [$w.term.e get]
	set hires [$w.res.hi.e get]
	set lores [$w.res.lo.e get]
	set shift_limit [$w.shift.e get]
	set thickness [$w.thick.e get]
	setMgParam
	if { ![file exists $paramfile] } {
#		puts "$paramfile not found!"
		tk_dialog .dialog "$paramfile not found!" \
			{Please save an appropriate parameter file first.} \
			info 0 {OK}
		return
	}
	set cct ""
	set rec ""
	set ref ""
	if { $cc_type } { set cct "-cross" }
	if { $recenter } { set rec "-recenter" }
	if { $refine } { set ref "-refine m" }
	if { $term < 0.1 } { set term 1 }
	set logfile "[file rootname $newfile].log"
	set psfile "[file rootname $newfile].ps"
#	puts "logfile = $logfile"
	set cmd "$env(BSOFT)/bin/btrack -verb 1 -reset -update -track $iter,$term -resol $hires,$lores -shift $shift_limit -thick $thickness $cct $rec $ref -Post $psfile -out $newfile $paramfile"
	executeCommand $cmd $logfile
}

## @brief Refines marker positions and micrograph geometry.
#

proc markerRefine { } {
	global helv12
	global hires lores ref_iter
	set w .wref
	if ![winfo exists $w] {
		toplevel $w
		
		frame $w.mark
		button $w.mark.b -text "Refine marker locations" -command "doMarkerRefine $w m"
		pack $w.mark.b -side left -pady 5

		setupResolutionEntry $w

		frame $w.aln
		button $w.aln.b -text "Refine alignment" -command "doMarkerRefine $w a"
		label $w.aln.tag -text Iterations
		tk_optionMenu $w.aln.iter ref_iter 1 2 3 4 5 6 7 8 9 10
		pack $w.aln.b $w.aln.tag $w.aln.iter -side left -padx 1m

		frame $w.z
		checkbutton $w.z.z -text "Marker z-coordinates" -variable ref_z
		
		frame $w.mg
		label $w.mg.g -text "Micrographs"
		checkbutton $w.mg.v -text "Views" -variable ref_view
		checkbutton $w.mg.o -text "Origins" -variable ref_ori
		checkbutton $w.mg.s -text "Scales" -variable ref_scale

		frame $w.iter

		pack $w.z.z -side left -pady 5
		pack $w.mg.g $w.mg.v $w.mg.o $w.mg.s -side left -pady 5
		pack $w.mark $w.res $w.aln $w.z $w.mg $w.iter \
			-side top -fill x -pady 2 -padx 2
		
		frame $w.buttons
#		button $w.buttons.do -text Refine -command "doMarkerRefine $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doMarkerRefine $w a"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Refine"
    wm iconname $w "Refine"
}

## @brief Does a refinement of markers or views.
#
# @param	w 			dialog window.

proc doMarkerRefine { w refop } {
	global ref_z ref_view ref_ori ref_scale
	global hires lores ref_iter
	setMgParam
	if { $refop == "m" } {
		set hires [$w.res.hi.e get]
		set lores [$w.res.lo.e get]
		Bmg refine $refop $hires $lores
	} else {
		set refop ""
		if $ref_z { append refop "z" }
		if $ref_view { append refop "v" }
		if $ref_ori { append refop "o" }
		if $ref_scale { append refop "s" }
		puts "doMarkerRefine: $refop $ref_iter"
#		Bmg refine $refop $ref_view $ref_ori $ref_scale
		Bmg refine $refop $ref_iter
	}
	updateOrigins
	updateTomoTable
	markerDrawAll
}

## @brief Refines one marker position.
#

proc refineMarker { } {
	global markers_selected
	global hires lores
	set hires 20
	set lores 1000
	set w .wref
	if [winfo exists $w] {
		set hires [$w.res.hi get]
		set lores [$w.res.lo get]
	}
	Bmg refine_one $markers_selected $hires $lores
	updateTomoTable
	markerDrawAll
}

## @brief Dialog box for calculating tilt-adjusted power spectra
#

proc tomoPowerSpectra { } {
	global paramfile theimg
	global volt Cs amp_fac defocus
	global tilex tiley tilez
	global hires lores
	checkParamFile
	set w .wps
	if ![winfo exists $w] {
		toplevel $w

		set newfile [newParamFile "_ps.star"]
		if ![string length $newfile] { destroy $w; return }

		frame $w.head
		label $w.head.lab -text \
			"The power spectra will be calculated using the currently saved parameter file:\n$paramfile"
		pack $w.head.lab

		setupNewFileEntry $w.file "Output parameter file" $newfile \
			"New parameter file written after calculating the power spectra"

		setupEntry $w.volt "Acceleration voltage" double $volt "Acceleration voltage in kV"
		setupEntry $w.cs "Spherical aberration" double $Cs "Spherical aberration in mm"
		setupEntry $w.amp "Amplitude fraction  " double $amp_fac "Amplitude contrast as a fraction"
		setupEntry $w.def "Defocus             " double $defocus "Average defocus in um"

		labelframe $w.tiles -text "Tile size:"
		setupEntry $w.tiles.x "x" integer $tilex "X-size for tiles to calculate power spectra"
		setupEntry $w.tiles.y "y" integer $tiley "Y-size for tiles to calculate power spectra"
		pack $w.tiles.x $w.tiles.y -side top -fill x

		labelframe $w.res -text "Resolution limits"
		setupEntry $w.res.hi "High" double 8 "High resolution limit in Å"
		setupEntry $w.res.lo "Low" double 30 "Low resolution limit in Å"
		pack $w.res.hi $w.res.lo -side top -fill x
		
		labelframe $w.deflim -text "Defocus limits"
		setupEntry $w.deflim.lo "Low" double 1.0 "Low defocus limit in µm"
		setupEntry $w.deflim.hi "High" double 5.0 "High defocus limit in µm"
		pack $w.deflim.lo $w.deflim.hi -side top -fill x
		
		pack $w.head $w.file $w.volt $w.cs $w.amp $w.def $w.tiles $w.res \
			$w.deflim -side top -fill x -pady 2 -padx 2

		frame $w.buttons
		button $w.buttons.do -text Do -command "doTomoPowerSpectra $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doTomoResolution $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Power spectrum"
    wm iconname $w "Power spectrum"
}

## @brief Calculating tilt-adjusted power spectra
#

proc doTomoPowerSpectra { w } {
	global env paramfile
	global volt Cs amp_fac defocus
	global tilex tiley tilez
	global hires lores recenter
	set newfile [$w.file.e get]
	set volt [$w.volt.e get]
	set Cs [$w.cs.e get]
	set amp_fac [$w.amp.e get]
	set defocus [$w.def.e get]
	set tilex [$w.tiles.x.e get]
	set tiley [$w.tiles.y.e get]
	set hires [$w.res.hi.e get]
	set lores [$w.res.lo.e get]
	set lodef [$w.deflim.lo.e get]
	set hidef [$w.deflim.hi.e get]
	setMgParam
	if { ![file exists $paramfile] } {
#		puts "$paramfile not found!"
		tk_dialog .dialog "$paramfile not found!" \
			{Please save an appropriate parameter file first.} \
			info 0 {OK}
		return
	}
	set logfile "[file rootname $newfile].log"
	set cmd "$env(BSOFT)/bin/bctf -verb 1,tim -act prepfit -tile $tilex,$tiley,1 -Volt $volt -Cs $Cs -Amp $amp_fac -Defocus $defocus -Range $lodef,$hidef -resol $hires,$lores -fitastig 2 -out $newfile $paramfile"
	executeCommand $cmd $logfile
}

## @brief Estimates the resolution of images in a tilt series.
#

proc tomoResolution { } {
	global paramfile theimg
	global hires samrat
	checkParamFile
	set w .wres
	if ![winfo exists $w] {
		toplevel $w

		frame $w.head
		label $w.head.lab -text \
			"Resolution will be estimated using the currently saved parameter file: $paramfile"
		pack $w.head.lab

		set psfile [newParamFile "_res.ps"]
		if ![string length $psfile] { destroy $w; return }
		setupPostscriptEntry $w.ps "Postscript file" $psfile \
			"Postscript file with the tilt-resolution plot"

		setupEntry $w.fast "Fast angle     " double 10 "Limit images added to each reconstruction"
		setupEntry $w.res "High resolution" double $hires "Spatial frequency limit for reconstructions"
		setupEntry $w.rat "Sampling ratio " double $samrat "Ratio to average adjacent frequency shells (≥1)"
		setupEntry $w.cut "FRC cutoff     " double 0.3 "Fourier ring correlation threshold to report resolution"
		setupEntry $w.ctf "CTF correction " alnum "none" "Contrast transfer function correction algorithm"

		pack $w.head $w.ps $w.fast $w.res $w.rat $w.cut $w.ctf -side top -fill x -pady 2 -padx 2

		frame $w.buttons
		button $w.buttons.do -text Do -command "doTomoResolution $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doTomoResolution $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Resolution"
    wm iconname $w "Resolution"
}

## @brief Executes estimation of the resolution of the tilt series.
#
# @param	w 			dialog window.

proc doTomoResolution { w } {
	global env paramfile
	global hires samrat
	set psfile [$w.ps.e get]
	set fast [$w.fast.e get]
	set hires [$w.res.e get]
	set samrat [$w.rat.e get]
	set cut [$w.cut.e get]
	set ctf [$w.ctf.e get]
	if { $ctf == "none" } {
		set ctf ""
	} else {
		set ctf "-CTF $ctf"
	}
	setMgParam
	if { ![file exists $paramfile] } {
#		puts "$paramfile not found!"
		tk_dialog .dialog "$paramfile not found!" \
			{Please save an appropriate parameter file first.} \
			info 0 {OK}
		return
	}
	set logfile "[file rootname $psfile].log"
	set cmd "$env(BSOFT)/bin/btomres -verb 1 -fast $fast -resol $hires -ratio $samrat $ctf -cut $cut -Post $psfile $paramfile"
	executeCommand $cmd $logfile
}

## @brief Reconstructs a tomogram.
#

proc tomoReconstruct { } {
	global paramfile theimg project_item
	global hires lores
	checkParamFile
	set w .wrec
	if ![winfo exists $w] {
		toplevel $w

		set sz [Bmg get $project_item size]
		set x [lindex $sz 0]
		set y [lindex $sz 1]
		set z [expr $x/10]
	
		set recfile [newParamFile "_rec.mrc"]
		if ![string length $recfile] { destroy $w; return }
		set newfile [newParamFile "_rec.star"]
		if ![string length $newfile] { destroy $w; return }
		
		frame $w.head
		label $w.head.lab -text \
			"Reconstruction will be done using the currently saved parameter file:\n$paramfile"
		
		setupFileEntry $w.rec "Output reconstruction " $recfile \
			"New tomographic reconstruction file name"

		setupNewFileEntry $w.file "Output parameter file " $newfile \
			"New parameter file written after reconstruction"

		setupEntry $w.interp "Type                 " integer 0 "Reconstruction interpolation type: 0=nearest neighbor, 1=weighted, 2=trilinear"

		labelframe $w.size -text "Reconstruction size"
		setupEntry $w.size.x "x" double $x "Reconstruction x size"
		setupEntry $w.size.y "y" double $y "Reconstruction y size"
		setupEntry $w.size.z "z" double $z "Reconstruction z size"
		pack $w.size.x $w.size.y $w.size.z -side top -fill x

		setupEntry $w.scale "Scale                 " double 1 "Reconstruction scale"
		setupEntry $w.res "High resolution       " double $hires "Spatial frequency limit for reconstruction"
		setupEntry $w.remark "Marker deletion radius" double 0 "If greater than 0, discs of this radius will be deleted from micrographs"
		setupEntry $w.edge "Edge smoothing width  " double 3 "If greater than 0, the micrograph edges and extraneous areas will be deleted"
		setupEntry $w.ctf "CTF correction        " alnum "none" "Contrast transfer function correction algorithm"

		frame $w.mem
		button $w.mem.b -text "Memory required" -command "tomoRecMemory $w"
		label $w.mem.v -text "" -relief sunken -bd 1 -width 16 -anchor w
		pack $w.mem.b $w.mem.v -side left

		pack $w.head.lab
#		pack $w.file.lab $w.file.entry $w.file.load -side left -pady 5
		pack $w.head $w.rec $w.file $w.interp $w.size $w.scale $w.res $w.remark $w.edge $w.ctf $w.mem -side top -fill x -pady 2 -padx 2

		frame $w.buttons
		button $w.buttons.do -text Do -command "doTomoReconstruct $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doTomoReconstruct $w"
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

proc tomoRecMemory { w } {
	set x [$w.size.x.e get]
	set y [$w.size.y.e get]
	set z [$w.size.z.e get]
	set mem [expr 20 * $x * $y * $z]
	set p [expr round ($mem * 100.0 / [systemMemory] ) ]
	set mem [ expr round($mem * 1.0e-9) ]
	$w.mem.v config -text "$mem Gb ($p %)"
}

## @brief Executes tomogram reconstruction.
#
# @param	w 			dialog window.

proc doTomoReconstruct { w } {
	global env paramfile
	global hires
	set recfile [$w.rec.e get]
	set newfile [$w.file.e get]
	set interp [$w.interp.e get]
	set x [$w.size.x.e get]
	set y [$w.size.y.e get]
	set z [$w.size.z.e get]
	set scale [$w.scale.e get]
	set hires [$w.res.e get]
	set remark [$w.remark.e get]
	set edge [$w.edge.e get]
	set ctf [$w.ctf.e get]
	set remove ""
	if { $remark > 0 } { set remove "-remove $remark" }
	set smooth ""
	if { $edge > 0 } { set smooth "-edge $edge" }
	if { $ctf == "none" } {
		set ctf ""
	} else {
		set ctf "-CTF $ctf"
	}
	setMgParam
	if { ![file exists $paramfile] } {
#		puts "$paramfile not found!"
		tk_dialog .dialog "$paramfile not found!" \
			{Please save an appropriate parameter file first.} \
			info 0 {OK}
		return
	}
	set logfile "[file rootname $newfile].log"
	set cmd "$env(BSOFT)/bin/btomrec -verb 1,tim -interp $interp -size $x,$y,$z -scale $scale -resol $hires $ctf -rescale 0,1 -trans full $remove $smooth -rec $recfile -out $newfile $paramfile"
	executeCommand $cmd $logfile
}

## @brief Denoises a tomogram.
#

proc tomoDenoise { } {
	global paramfile theimg project_item
	global hires lores
	checkParamFile
	set w .wrec
	if ![winfo exists $w] {
		toplevel $w

		set denfile [file rootname [Bmg get $project_item filename rec]]
#		set denfile [file rootname $paramfile]
		append denfile "_den.mrc"
		set newfile [newParamFile "_den.star"]
		if ![string length $newfile] { destroy $w; return }
		
		frame $w.head
		label $w.head.lab -text \
			"Denoising will be done using the currently saved parameter file: \n$paramfile"
		
		frame $w.den
		label $w.den.l -text      "Output denoised map  " -anchor w
		entry $w.den.e -width 40
		$w.den.e insert 0 $denfile
		pack $w.den.l $w.den.e -side left -pady 5

		setupNewFileEntry $w.file "Output parameter file" $newfile \
			"New parameter file written after denoising"

		setupEntry $w.iter "Iterations      " integer 100 "Number of iterative denoising cycles"
		setupEntry $w.out "Output frequency" integer 10 "Number of iterations between writing intermediate denoised images"
		setupEntry $w.slab "Slab size       " integer 100 "Denoising slab size for parallel processing"
		
		pack $w.head.lab
#		pack $w.log.lab $w.log.entry -side left -pady 5
		pack $w.head $w.den $w.file $w.iter $w.out $w.slab -side top -fill x -pady 2 -padx 2

		frame $w.buttons
		button $w.buttons.do -text Do -command "doTomoDenoise $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	
		bind $w <Return> "doTomoDenoise $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Denoise"
    wm iconname $w "Denoise"
}

## @brief Executes tomogram denoising.
#
# @param	w 			dialog window.

proc doTomoDenoise { w } {
	global env paramfile project_item
	global hires
	set recfile [Bmg get $project_item filename rec]
	set denfile [$w.den.e get]
	set newfile [$w.file.e get]
	set iter [$w.iter.e get]
	set out [$w.out.e get]
	set slab [$w.slab.e get]
	setMgParam
	if { ![file exists $paramfile] } {
#		puts "$paramfile not found!"
		tk_dialog .dialog "$paramfile not found!" \
			{Please save an appropriate parameter file first.} \
			info 0 {OK}
		return
	}
	set logfile "[file rootname $newfile].log"
#	set cmd "$env(BSOFT)/bin/bnad -verb 1 -iter $iter -out $out -slab $slab -map $denfile -out $newfile $paramfile"
	set cmd "$env(BSOFT)/bin/bnad -verb 1 -iter $iter -out $out -slab $slab $recfile $denfile"
	executeCommand $cmd $logfile
}


