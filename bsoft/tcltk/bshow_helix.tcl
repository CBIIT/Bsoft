##
# @file		bshow_helix.tcl
#
# @brief	Procedures to manage helical diffraction patterns
#
# @author	Bernard Heymann
# @date		Created: 20070516
# @date		Modified: 20160203

set show_ll 1
set show_ll_indices 1
set ll_length 100
set ll_width 10
set ll_resolution 20
set ll_selected -1
set ll_plot_width 300
set ll_plot_height 100

## @brief Dialog box to manage helix information
#

proc Helix { } {
	global project_item theimg
	global helv12
	global ll_width ll_length
	global ll_plot_width ll_plot_height

	if ![Bimage exists $theimg] {
		tk_dialog .dialog "No image in memory!" \
				{Please read in an appropriate image first.} \
				info 0 {OK}
		return
	}
	
	set c [getImageCanvas $theimg]

	if ![Bmg exists] {
		createMgParam
	}
	set rad [Bmg get $project_item helix_radius]
	set axis [Bmg get $project_item helix_axis]
#	puts "axis = $axis"
	if { $rad < 1 } { set rad 100 }
	if { $axis < 0 } { set axis 90 }
	if ![winfo exists .whel] {
		toplevel .whel
		menu .whel.menuBar -tearoff 0
		.whel.menuBar add cascade -menu .whel.menuBar.tomo -label "Helix" -underline 0
		menu .whel.menuBar.tomo -tearoff 0
		.whel configure -menu .whel.menuBar
		.whel.menuBar.tomo add command -label "Help" \
				-command { showURL "file://$Bsoft/doc/bshow/bshow_helix.html" }
#				-command { showHelp bshow_xtal.hlp } -underline 0
		.whel.menuBar.tomo add command -label "Set layer line index and order" \
				-command { setLayerLineIndex } -underline 0
		.whel.menuBar.tomo add command -label "Generate layer lines" \
				-command { layerLinesGenerate } -underline 0
		.whel.menuBar.tomo add command -label "Mask layer lines" \
				-command { layerLinesMask } -underline 0
		.whel.menuBar.tomo add command -label "Delete all layer lines" \
				-command { layerLineDeleteAll } -underline 0
		.whel.menuBar.tomo add command -label "Close" \
				-command "destroy .whel" -underline 0
		
		frame .whel.switches
		checkbutton .whel.switches.ll -text "Show layer lines" -variable show_ll \
				-command { layerLinesDrawAll }
		checkbutton .whel.switches.labels -text "Show indices" -variable show_ll_indices \
				-command { layerLinesDrawAll }
		
		frame .whel.axis
		label .whel.axis.tag -font $helv12 -text "Helix axis angle" -width 20 -anchor e
		entry .whel.axis.entry -width 10 -validate key -vcmd { string is double %P }
		.whel.axis.entry insert 0 $axis

		frame .whel.radius
		label .whel.radius.tag -font $helv12 -text "Helix radius" -anchor e
		entry .whel.radius.entry -width 4 -validate key -vcmd { string is double %P }
		.whel.radius.entry insert 0 $rad
		
		frame .whel.ll_length
		label .whel.ll_length.tag -font $helv12 -text "Layer line length" -width 20 -anchor e
		entry .whel.ll_length.entry -width 10 -validate key -vcmd { string is double %P }
		.whel.ll_length.entry insert 0 $ll_length
		
		frame .whel.ll_width
		label .whel.ll_width.tag -font $helv12 -text "Layer line width" -width 20 -anchor e
		entry .whel.ll_width.entry -width 10 -validate key -vcmd { string is double %P }
		.whel.ll_width.entry insert 0 $ll_width
		
		frame .whel.number
		label .whel.number.tag -font $helv12 -text "Layer line count" -width 20 -anchor e
		label .whel.number.num -width 10 -text "0" -relief sunken -bd 1 \
				-font $helv12 -anchor w

		frame .whel.selected
		label .whel.selected.tag -font $helv12 -text "Selected" -width 10
		label .whel.selected.index -width 4 -text "" -relief sunken -bd 1 -font $helv12
		label .whel.selected.order_tag -width 5 -font $helv12 -text "order"
		entry .whel.selected.order_entry -width 4
		label .whel.selected.fom_tag -width 5 -font $helv12 -text "FOM"
		label .whel.selected.fom -width 10 -text "" -relief sunken -bd 1 -font $helv12

		frame .whel.frame
		text .whel.table -relief sunken -bd 2 -width 40\
			-yscrollcommand ".whel.yscroll set" \
			-setgrid 1 -highlightthickness 0 -pady 2 -padx 3
		scrollbar .whel.yscroll -command ".whel.table yview" \
			-highlightthickness 0 -orient vertical
		grid .whel.table -in .whel.frame -padx 1 -pady 1 \
			-row 0 -column 0 -rowspan 1 -columnspan 1 -sticky news
		grid .whel.yscroll -in .whel.frame -padx 1 -pady 1 \
			-row 0 -column 1 -rowspan 1 -columnspan 1 -sticky news
		grid rowconfig    .whel.frame 0 -weight 1 -minsize 0
		grid columnconfig .whel.frame 0 -weight 1 -minsize 0

		set g .whel.plot
		canvas $g -width [expr $ll_plot_width + 10] -height [expr $ll_plot_height + 10] \
				-background white -relief sunken -borderwidth 2

		pack .whel.axis.tag .whel.axis.entry -side left -pady 5 -padx 10
		pack .whel.radius.tag .whel.radius.entry -side left -pady 5 -padx 10
		pack .whel.ll_length.tag .whel.ll_length.entry -side left -pady 5 -padx 10
		pack .whel.ll_width.tag .whel.ll_width.entry -side left -pady 5 -padx 10
		pack .whel.number.tag .whel.number.num -side left -pady 5 -padx 10
		pack .whel.selected.tag .whel.selected.index .whel.selected.order_tag \
				.whel.selected.order_entry .whel.selected.fom_tag .whel.selected.fom \
				-side left -pady 5 -padx 5
		pack .whel.switches.ll .whel.switches.labels -side left -ipadx 5 -ipady 5 -pady 5 -padx 5
		pack .whel.switches .whel.axis .whel.radius .whel.ll_length .whel.ll_width .whel.number \
				.whel.selected .whel.plot -side top -pady 5 -padx 5
#		pack .whel.frame -expand yes -fill both -padx 5 -pady 5
		
#		updateLayerLineTable
		updateLayerLineParam
	} else {
		wm deiconify .whel
		raise .whel
    }
    wm title .whel "Helix"
    wm iconname .whel "Helix"
	.menuBar.window entryconfigure "Helix" -state normal
	bind $c <1> "layerLineCreateSelect %x %y"
	bind $c <2> "layerLineCreateSelect %x %y"
	bind $c <Shift-1> "layerLineDelete %x %y"
	bind $c <Shift-2> "layerLineDelete %x %y"
	bind $c <B1-Motion> "layerLineMove %x %y"
	bind $c <B2-Motion> "layerLineMove %x %y"
	bind .whel <Return> "Update 0"
	bind .whel <Control-w> { destroy .whel }
	bind .whel <Command-w> { destroy .whel }
#	updateLayerLineTable
	layerLinesDrawAll
	drawHelixAxis
}

## @brief Dialog box to set the index for the selected layer line

proc setLayerLineIndex { } {
	global helv12 theimg
	global ll_selected
	set wc [getControlWindow $theimg]
	set w .wlli
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]

	if ![winfo exists $w] {
		toplevel $w
		
		frame $w.index
		label $w.index.tag -font $helv12 -text "Index" -width 14
		entry $w.index.entry -width 6 -validate key \
			-vcmd { expr {[string match {[-+]} %P] || [string is int %P]} }

		frame $w.order
		label $w.order.tag -font $helv12 -text "Order" -width 14
		entry $w.order.entry -width 6 -validate key -vcmd { string is int %P }

		if { $ll_selected >= 0 } {
			$w.index.entry insert 0 $ll_selected
			$w.order.entry insert 0 [Bmg layerline $mg_item order $ll_selected]
		}

		pack $w.index.tag $w.index.entry -side left -padx 5 -pady 5
		pack $w.order.tag $w.order.entry -side left -padx 5 -pady 5
		pack $w.index $w.order -side top -fill x -pady 2 -padx 15
		
		frame $w.buttons
		button $w.buttons.ok -text OK -command "setLayerLineIndex_and_DestroyWindow $w"
		button $w.buttons.cancel -text Cancel -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.ok $w.buttons.cancel -side left -expand 1
	
		bind $w <Return> "setLayerLineIndex_and_DestroyWindow $w"
	} else {
		wm deiconify $w
		raise $w
    }
	
    wm title $w "Index"
    wm iconname $w "Index"
	tkwait window $w
}

## @brief Sets the index for a layer line.
#
# @param	w 			the window to destroy.

proc setLayerLineIndex_and_DestroyWindow { w } {
	global theimg
	global ll_selected
	set wc [getControlWindow $theimg]
	set num [$w.index.entry get]
	set order [$w.order.entry get]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	if { $ll_selected < 0 } {
		set ids [Bmg layerline $mg_item ids]
		foreach id $ids {
			if { $id == $num } {
				tk_messageBox -icon error -type ok -title "Error" -message \
					"Layer line $num already exists!"
				return
			}
		}
		set ll_selected $num
	} else {
		set ll_selected [Bmg layerline $mg_item renumber $ll_selected $num]
	}
	.whel.selected.index config -text "$ll_selected"
	.whel.selected.order_entry delete 0 end
	.whel.selected.order_entry insert 0 $order
	Bmg layerline $mg_item order $ll_selected $order
	destroy $w
	layerLinesDrawAll
}

## @brief Creates or selects a layer line.
#
# @param	x
# @param	y			Layer line position in canvas coordinates.

proc layerLineCreateSelect { x y } {
	if ![winfo exists .whel] { return }
	global theimg
	global ll_selected
	global tool
	global px py
	set ll_selected [layerLineSelect $x $y]
#	puts "Layer line selected = $ll_selected"
	if { $ll_selected > -1 } { return }
	if { $tool != "helix" } { return }
	setLayerLineIndex
	if { $ll_selected < 0 } { return }
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set scale [$wc.scale.scale get]
	set px [expr floor( [$c canvasx $x] / $scale ) ]
	set py [expr floor( [$c canvasy $y] / $scale ) ]
	set pyi [yFlip $py]
#	puts "Layer line created = $ll_selected"
	set mg_item [micrographItem $n]
#	puts "Bmg layerline $mg_item create $ll_selected $px $pyi"
	set ll_selected [Bmg layerline $mg_item create $ll_selected $px $pyi]
	puts "Layer line created = $ll_selected"
	layerLineDraw $ll_selected
	.whel.number.num config -text [Bmg layerline $mg_item count]
}

## @brief Selects a layer line.
#
# @param	x
# @param	y			Layer line position in canvas coordinates.

proc layerLineSelect { x y } {
	if ![winfo exists .whel] { return }
	global theimg
	global ll_selected
	global fom_cutoff ll_width
	global px py
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set ll_selected -1
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	if { [Bmg layerline $mg_item count] < 1 } { return $ll_selected }
	set scale [$wc.scale.scale get]
	set px [expr floor( [$c canvasx $x] / $scale ) ]
	set py [expr floor( [$c canvasy $y] / $scale ) ]
	set pyi [yFlip $py]
	set ll_selected [Bmg layerline $mg_item select $px $pyi $ll_width]
	if { $ll_selected > -1 } {
		set fom [Bmg layerline $mg_item fom $ll_selected]
		set order [Bmg layerline $mg_item order $ll_selected]
		.whel.selected.index config -text "$ll_selected"
		.whel.selected.fom config -text "$fom"
		.whel.selected.order_entry delete 0 end
		.whel.selected.order_entry insert 0 $order
		layerLineDraw $ll_selected
		layerLinePlot $ll_selected $order [.whel.radius.entry get]
	}
	return $ll_selected
}

## @brief Moves the selected layer line location.
#
# @param	x
# @param	y			Layer line position in canvas coordinates.

proc layerLineMove { x y } {
	if ![winfo exists .whel] { return }
	global theimg
	global ll_selected
	global px py
	if { $ll_selected < 0 } { return }
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set scale [$wc.scale.scale get]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set ix [expr [$c canvasx $x] / $scale ]
	set iy [expr [$c canvasy $y] / $scale ]
	set dx [expr $ix - $px]
	set dy [expr $py - $iy]
	set px $ix
	set py $iy
#	puts "Moving $dx $dy"
	Bmg layerline $mg_item move $ll_selected $dx $dy
	layerLinesDrawAll
}

## @brief Draws a line along the helix axis on the canvas
#

proc drawHelixAxis { } {
	global theimg
	if ![winfo exists .whel] { return }
	global project_item
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set scale [$wc.scale.scale get]
	set axis [.whel.axis.entry get]
	set width [Bimage get $theimg width]
	set height [Bimage get $theimg height]
	set ox [getOriginX $c]
	set oy [getOriginY $c]
	set len $width
	if { $len < $height } { set len $height }
	set len [expr $scale * $len/2]
	set cosax [expr cos($axis)]
	set sinax [expr sin($axis)]
	set xmin [expr $ox - $len*$cosax]
	set ymin [expr $oy + $len*$sinax]
	set xmax [expr $ox + $len*$cosax]
	set ymax [expr $oy - $len*$sinax]
	$c delete helix_axis
	$c create line $xmin $ymin $xmax $ymax -fill LightBlue -width 1.5 -tags helix_axis
	Bmg set $project_item helix_axis $axis
}

## @brief Draws a layer line in the original image on the canvas
#
# @param	id			Layer line index.

proc layerLineDraw { id } {
	if ![winfo exists .whel] { return }
	global theimg
	global show_ll_indices
	global ll_length
	global ll_selected ll_color ll_select_color
	if ![winfo exists .whel] { return }
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	if { $id < 0 } { return }
	if { [llength $id] < 1 } {
		puts "No ID specified!"
		return
	}
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set d [Bmg layerline $mg_item distance $id]
	if { [llength $d] < 1 } {
		puts "No distance for ID = $id!"
		return
	}
	set scale [$wc.scale.scale get]
	set d [expr $d * $scale ]
#	puts "Drawing $id at $d"
	set ox [getOriginX $c]
	set oy [getOriginY $c]
	set axis [.whel.axis.entry get]
	set cosax [expr cos($axis)]
	set sinax [expr sin($axis)]
	set ll_length [.whel.ll_length.entry get]
	set x [expr $ox + $d * $cosax]
	set y [expr $oy - $d * $sinax]
	set len [expr $scale * $ll_length]
	set xmin [expr $x - $len*$sinax]
	set ymin [expr $y - $len*$cosax]
	set xmax [expr $x + $len*$sinax]
	set ymax [expr $y + $len*$cosax]
#	puts "$xmin $ymin $xmax $ymax"
	set color $ll_color
	if { $id == $ll_selected } { set color $ll_select_color }
	$c create line $xmin $ymin $xmax $ymax -fill $color -dash . -width 1 -tags ll
	if { $show_ll_indices } {
		set x [expr $x + 10]
		$c create text $x $y -text $id -fill $color -anchor w -tags ll
	}
}

## @brief Draws all layer lines on the canvas.
#

proc layerLinesDrawAll { } {
	global theimg
	global fom_cutoff
	global show_ll
	set c [getImageCanvas $theimg]
	$c delete ll
	if ![winfo exists .whel] { return }
	if { $show_ll < 1 } { return }
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set ids [Bmg layerline $mg_item ids $fom_cutoff]
#	puts "drawing IDS = $ids"
	set nsel 0
	foreach id $ids {
#		puts "drawing ID = $id"
		layerLineDraw $id
		incr nsel 1
	}
	.whel.number.num config -text "$nsel"
}

## @brief Deletes the selected layer line on the canvas.
#
# @param	x
# @param	y			Layer line position in window coordinates.

proc layerLineDelete { x y } {
	global theimg ll_selected
	set ll_selected [layerLineSelect $x $y]
	if { $ll_selected < 0 } { return }
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	Bmg layerline $mg_item delete $ll_selected
	layerLinesDrawAll
}


## @brief Deletes all layer lines in list.
#

proc layerLineDeleteAll { } {
	global theimg
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	Bmg layerline $mg_item delete all
	$c delete ll
	.whel.number.num config -text "0"
}

## @brief Generates layer lines from unit cell vectors.
#

proc layerLinesGenerate { } {
	global ll_resolution
	set w .whelres
	catch {destroy $w}
	toplevel $w
	wm title $w "Set the resolution limit"
	wm iconname $w "Resolution"
	tkwait visibility $w
	grab $w
	frame $w.buttons
	button $w.buttons.ok -text OK -command "layerLinesGenerate_and_DestroyWindow $w"
	button $w.buttons.cancel -text Cancel -command "destroy $w"
	frame $w.res
	label $w.res.tag -text "Resolution limit" -anchor e
	entry $w.res.entry -width 10 -validate key -vcmd { string is double %P }
	$w.res.entry insert 0 "$ll_resolution"
	pack $w.res.tag $w.res.entry -side left -padx 5 -pady 5
	pack $w.res -side top -fill x -pady 2
	pack $w.buttons -side bottom -fill x -pady 2m
	pack $w.buttons.ok $w.buttons.cancel -side left -expand 1
	bind $w <Return> "layerLinesGenerate_and_DestroyWindow $w"
	tkwait window $w
}

## @brief Generates layer lines from unit cell vectors.
#
# @param	w 			the window to destroy.

proc layerLinesGenerate_and_DestroyWindow { w } {
	global theimg
	set wc [getControlWindow $theimg]
	set resolution [$w.res.entry get]
	set pixel_size [$wc.pixel_size.x get]
	set nx [Bimage get $theimg width]
	set ny [Bimage get $theimg height]
	if { $nx > $ny } {
		set ll_lim [expr int($nx*$pixel_size/$resolution)]
	} else {
		set ll_lim [expr int($ny*$pixel_size/$resolution)]
	}
#	puts "Radial limit = $ll_lim"
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set nll [Bmg layerline $mg_item generate $ll_lim]
#	puts "$nll layer lines generated"
	destroy $w
	updateLayerLineParam
	layerLinesDrawAll
}

## @brief Masks layer lines.
#

proc layerLinesMask { } {
	global theimg ll_width
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	Bmg layerline $mg_item mask $ll_width
	Update 1
}

## @brief Plots the selected layer line in the plot canvas
#
# @param	id			Layer line index.
# @param	order		Layer line order.
# @param	rad			Layer line length.

proc layerLinePlot { id order rad } {
	global theimg project_item
	global ll_plot_width ll_plot_height
	if ![winfo exists .whel] { return }
	Bmg set $project_item helix_radius $rad
	if { $id < 0 } { return }
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	if { [Bmg layerline $mg_item count] < 1 } { return }
	set ll_len [.whel.ll_length.entry get]
	set data [Bmg layerline $mg_item plot $id $ll_len]
	set bessel [Bmg layerline $mg_item bessel $order $rad $ll_len]
	set len [lindex $data 0]
	if { $len < 1 } { return }
	set max 0
	set bmax 0
	for {set i 1} {$i <= $len} {incr i 1} {
		set v [lindex $data $i]
		if { $max < $v } { set max $v }
		set b [lindex $bessel $i]
		if { $bmax < $b } { set bmax $b }
	}
#	set max [expr $max/$bmax]
#	set max [Bimage get $theimg max]
	set xoff 10
	set yoff 10
	set xscale [expr $ll_plot_width*1.0/$len]
	set yscale [expr $ll_plot_height*1.0/$max]
	set ybscale [expr $ll_plot_height*1.0/$bmax]
	set g .whel.plot
	$g delete llplot
	set ll ""
	set bl ""
	for {set i 1} {$i <= $len} {incr i 1} {
		set x [expr ($i-1)*$xscale + $xoff]
		set y [expr ($max - [lindex $data $i]) * $yscale + $yoff]
#		set b [expr (1 - [lindex $bessel $i]) * $ll_plot_height + $yoff]
		set b [expr ($bmax - [lindex $bessel $i]) * $ybscale + $yoff]
		append ll " " $x " " $y
		append bl " " $x " " $b
	}
	$g create line $ll -fill blue -smooth yes -width 1 -tags llplot
	$g create line $bl -fill red -smooth yes -width 1 -tags llplot
	set orix [expr $ll_plot_width/2 + $xoff]
	set oriy [expr $ll_plot_height + 2*$yoff]
	$g create line $orix 0 $orix $oriy \
		-fill black -smooth yes -width 1 -tags llplot
}

## @brief Updates layer line parameters from the micrograph parameters in memory.
#

proc updateLayerLineParam { } {
	global project_item imgtype
	global ll_width ll_selected
	set w .whel
	if ![winfo exists $w] { return }
	if ![Bmg exists] { createMgParam }
	set filename [Bmg get $project_item filename $imgtype]
	wm title $w [list "Helix:" $filename]
#	puts "Updating layer line parameters"
#	set pixel_size [Bmg get $project_item pixel_size]
	set axis [Bmg get $project_item helix_axis]
	if !$axis { set axis 90 }
	$w.axis.entry delete 0 end
	$w.axis.entry insert 0 $axis
#	set ll_width [Bmg get $project_item ll_width]
	updatePixelsize
	$w.ll_width.entry delete 0 end
	$w.ll_width.entry insert 0 $ll_width
	set ll_selected -1
}

