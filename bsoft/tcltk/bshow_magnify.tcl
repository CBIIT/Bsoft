##
# @file		bshow_magnify.tcl
#
# @brief	Procedures for rendering the magnifying window in Bshow.
#
# @author	Bernard Heymann
# @date		Created: 20010516
# @date		Modified: 20161102

set magnification 2
set magsize 75
set maglock 0
set aniproc2 0
set aniaxis z
set mag3D 0

## @brief Window to show a region magnified. 
#
# @param	x
# @param	y			Coordinates on the canvas.

proc Magnify { x y } {
	global helv12 theimg
	global magnification magsize maglock show_cross

	if ![Bimage exists $theimg] {
		tk_dialog .dialog "No image in memory!" \
				{Please read in an appropriate image first.} \
				info 0 {OK}
		return
	}
	
	set wc [getControlWindow $theimg]
	set w .wmag

	if ![winfo exists $w] {
		toplevel $w
		
		set nslices [Bimage get $theimg nslices]
		
		frame $w.frame -relief ridge -borderwidth 2 -bg black

		if { $nslices < 2 } {
			wm minsize $w 200 200
			set m $w.frame.c
			canvas $m -width 200 -height 200
			
			grid $m -in $w.frame -sticky nw -padx 2 -pady 2
			grid rowconfig    $w.frame 0 -weight 1 -minsize 0
			grid columnconfig $w.frame 0 -weight 1 -minsize 0

			image create photo theimg0
			$m create image 0 0 -image theimg0 -anchor nw -tags magimage
		} else {
			wm minsize $w 300 300
			set m $w.frame.c
			set m1 $w.frame.c1
			set m2 $w.frame.c2
			canvas $m -width 150 -height 150
			canvas $m1 -width 150 -height 150
			canvas $m2 -width 150 -height 150

			grid $m $m1 -in $w.frame -sticky nw -padx 2 -pady 2
			grid $m2 -in $w.frame -sticky nw -padx 2 -pady 2
			grid rowconfig    $w.frame 0 -weight 1 -minsize 0
			grid columnconfig $w.frame 0 -weight 1 -minsize 0

			image create photo theimg0
			$m create image 0 0 -image theimg0 -anchor nw -tags magimage
			image create photo theimg1
			$m1 create image 0 0 -image theimg1 -anchor nw -tags magimage
			image create photo theimg2
			$m2 create image 0 0 -image theimg2 -anchor nw -tags magimage
		}
		
		frame $w.n
		label $w.n.x_tag -font $helv12 -text X
		entry $w.n.x -width 6 -validate key -vcmd { string is double %P }
		label $w.n.y_tag -font $helv12 -text Y
		entry $w.n.y -width 6 -validate key -vcmd { string is double %P }
		checkbutton $w.n.cross -text "Crosses" -variable show_cross
		set show_cross 1
		
		frame $w.p
		label $w.p.size_tag -font $helv12 -text "Size"
		entry $w.p.size -width 6 -validate key -vcmd { string is double %P }
		$w.p.size insert 0 $magsize
		label $w.p.mag_tag -font $helv12 -text "Magnification"
		entry $w.p.mag -width 6 -validate key -vcmd { string is double %P }
		$w.p.mag insert 0 $magnification
		
		frame $w.q 
		label $w.q.slice_tag -font $helv12 -text "Slice" 
		scale $w.q.slice -from 0 -to [expr [Bimage get $theimg nslices] - 1] -orient horizontal -length 175 \
				-command changeSlice
		$w.q.slice set [$wc.slice.scale get]
		
		frame $w.r
		label $w.r.range_tag -font $helv12 -text "Range"
		entry $w.r.range -width 3 -validate key -validatecommand {check_length %P 3}
		$w.r.range insert 0 5
		label $w.r.speed_tag -font $helv12 -text "Speed"
		entry $w.r.speed -width 3 -validate key -validatecommand {check_length %P 3} 
		$w.r.speed insert 0 15
		button $w.r.animate -text Animate \
			-relief raised -command animateSlices -font $helv12
		tk_optionMenu $w.r.axis aniaxis x y z
		
		frame $w.instr
		label $w.instr.lock -font $helv12 -text "Press l to lock into xy position"
		label $w.instr.mag -font $helv12 -text "Press m to decrease magnification"
		
		pack $w.n.x_tag $w.n.x $w.n.y_tag $w.n.y $w.n.cross \
				-side left -ipadx 3 -ipady 2
		pack $w.p.size_tag $w.p.size $w.p.mag_tag $w.p.mag \
				-side left -ipadx 3 -ipady 2
		pack $w.q.slice_tag $w.q.slice \
				-side left -ipadx 3 -ipady 2
		pack $w.r.range_tag $w.r.range \
				$w.r.speed_tag $w.r.speed $w.r.animate \
				$w.r.axis -side left -ipadx 2 -ipady 2
		pack $w.q $w.r $w.p $w.n -side bottom -fill x
		pack $w.instr.lock $w.instr.mag -side left -ipadx 10
		pack $w.instr $w.frame -side top -fill both -expand yes
		
		
		wm title $w "Magnify"
		wm iconname $w "Magnify"
	} else {
		wm deiconify $w
		raise $w
    }
	bind .wmag <Control-w> { destroy .wmag }
	bind .wmag <Command-w> { destroy .wmag }
	bind . <Key-l> { set maglock [expr 1 - $maglock] }
	bind . <Key-m> {
		if { $magnification > 0.1 } { set magnification [expr $magnification / 2.0] }
		.wmag.p.mag delete 0 end
		.wmag.p.mag insert 0 $magnification
	}
}

proc drawCross { c lx ly } {
	$c create line [expr $lx - 5] $ly [expr $lx + 5] $ly -fill yellow -width 1 -tags cross
	$c create line $lx [expr $ly - 5] $lx [expr $ly + 5] -fill yellow -width 1 -tags cross
}

## @brief Updates the magnifying window.
#
# @param	x
# @param	y				Coordinates within the image.
# @param	currSlice		Current slice (defaults to main window's slice unless animating)
# @param	updatesource	Procedure calling update (mouse or slicer)

proc updateMagnify { x y currSlice updatesource } {
	global theimg fom_cutoff
	global mag3D show_cross project_item
	global magnification magsize maglock
	global show_boxes show_bad show_filaments
	global fom_cutoff show_markers
	global show_comp show_comp_labels comp_selected line_width
	global show_box_labels box_color box_select_color bad_color bad_select_color
	global filament_color filament_first_color filament_select_color
	
	set w .wmag
	if ![winfo exists $w] { return }
	if { $maglock && $updatesource=="mouse" } { return }
	if { $x > [Bimage get $theimg width] } { return }
	if { $y < 0 } { return }

	set x [expr int($x)]
	set y [expr int($y)]
	
	set wc [getControlWindow $theimg]
	set scale [$wc.scale.scale get]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	
	if { $currSlice=="default" } {
		set z [$wc.slice.scale get]
	} else {
		set z $currSlice
	}

	set magnification [$w.p.mag get]
	set magsize [$w.p.size get]
	set sx $magsize
	set sy $magsize
	set sz $magsize
	set width [expr int($magnification*$sx)]
	set height [expr int($magnification*$sy)]
	set nslices [Bimage get $theimg nslices]
	
	$w.frame.c delete cross box bad node marker component link
	$w.frame.c configure -width $width -height $height -scrollregion "0 0 $width $height"

	set mag3D 0
	if { $nslices > 1 } {
		set mag3D 1
		$w.frame.c1 delete cross box bad node marker component link
		$w.frame.c2 delete cross box bad node marker component link
		$w.frame.c1 configure -width $width -height $height -scrollregion "0 0 $width $height"
		$w.frame.c2 configure -width $width -height $height -scrollregion "0 0 $width $height"
	}
	wm minsize $w [expr 2*$width+60] [expr 2*$height+210]
	Bimage magnify $theimg $n $x $y $z $magnification $sx $sy $sz
#	puts "Mag image extracted"

	set lx [expr $width/2]
	set ly [expr $height/2]
	if { $show_cross } {
		drawCross $w.frame.c $lx $ly
		if { $mag3D } {
			drawCross $w.frame.c1 $lx $ly
			drawCross $w.frame.c2 $lx $ly
		}
	}
	if { $maglock==0 || $updatesource=="slicer" } {
		$w.n.x delete 0 end
		$w.n.x insert 0 $x
		$w.n.y delete 0 end
		$w.n.y insert 0 $y
	}

	if [winfo exists .wbox] {
		if { $show_boxes } {
			set ids [Bmg box $mg_item ids $fom_cutoff]
			foreach id $ids {
				drawMagBox $id $x $y $z
			}
		}
		if { $show_bad } {
			set ids [Bmg box $mg_item ids bad]
			foreach id $ids {
				drawMagBox $id $x $y $z
			}
		}
	}

	if [winfo exists .wfil] {
		if { $show_filaments } {
			set rad [.wfil.node_radius.entry get]
			set ids [Bmg node $mg_item ids]
			set pfid 0
			set pnid 0
			foreach {fid nid} $ids {
				set loc [Bmg node $mg_item location $fid $nid]
				if { $pfid == $fid } {
					set color $filament_color
				} else {
					set color $filament_first_color
				}
				drawMagSphere [lindex $loc 0] [lindex $loc 1] [lindex $loc 2] \
					 $x $y $z $rad $color $line_width $nid node
				if { $pfid == $fid } {
					drawMagLine [lindex $loc 0] [lindex $loc 1] [lindex $loc 2] \
						[lindex $ploc 0] [lindex $ploc 1] [lindex $ploc 2] \
						 $x $y $z 1 $filament_color link
				}
				set pfid $fid
				set pnid $nid
				set ploc $loc
			}
		}
	}
	
	if [winfo exists .wtomo] {
		if { $show_markers } {
			set ids [Bmg marker $mg_item ids $fom_cutoff]
			foreach id $ids {
				drawMagMarker $id $x $y $z
			}
		}
	}

	if [winfo exists .wmod] {
		if { $show_comp } {
			set comps [Bmodel component $project_item img_coords $fom_cutoff]
			if { $show_comp_labels } {
				foreach {id cx cy cz rad color} $comps {
					drawMagSphere $cx $cy $cz $x $y $z $rad $color $line_width $id component
				}
			} else {
				foreach {id cx cy cz rad color} $comps {
					drawMagSphere $cx $cy $cz $x $y $z $rad $color $line_width "" component
				}
			}
			set links [Bmodel link $project_item img_coords $fom_cutoff]
			foreach {x1 y1 z1 x2 y2 z2 rad color} $links {
				drawMagLine $x1 $y1 $z1 $x2 $y2 $z2 $x $y $z $rad $color link
			}
		}
	}
}


## @brief Draws a box in the magnification window.
#
# @param	id		Box identifier.
# @param	ptx
# @param	pty		Coordinates within the image.
# @param	ptz		Current slice (defaults to main window's slice unless animating)

proc drawMagBox { id ptx pty ptz } {
	global theimg

	set w .wmag
	if ![winfo exists $w] { return }

	set wb .wbox
	if ![winfo exists $wb] { return }
	
	global project_item
	global show_box_labels box_color box_select_color bad_color bad_select_color

	set wc [getControlWindow $theimg]

	set magnification [$w.p.mag get]
	set magsize [$w.p.size get]
	set slice_num [$wc.slice.scale get]
	set h [Bimage get $theimg height]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set loc [Bmg box $mg_item location $id]
	set x [lindex $loc 0]
#	set y [expr ($h - [lindex $loc 1] - 1) ]
	set y [yFlip [lindex $loc 1]]
	set z [lindex $loc 2]
	
	if { $id > 0 } {
		set xrad [expr [$wb.box_size.x get]*0.5]
		set yrad [expr [$wb.box_size.y get]*0.5]
		set zrad [expr [$wb.box_size.z get]*0.5]
	} else {
		set xrad [expr [$wb.bad_radius.e get]*1.0]
		set yrad $xrad
		set zrad $xrad
	}
	
	set pty [expr $h - $pty]
	set nslices [Bimage get $theimg nslices]
	
	set xr [expr abs($ptx - $x)]
	set yr [expr abs($pty - $y)]
	set zr [expr abs($ptz - $z)]
	
	set xdiff [expr ($ptx - $x)*$magnification]
	set ydiff [expr ($pty - $y)*$magnification]
	set zdiff [expr ($ptz - $z)*$magnification]

	set lx [expr $magnification*$magsize/2]
	set ly [expr $magnification*$magsize/2]
	
	set idtag [format "box%d" $id]

	if { $zr <= $zrad } {
		set fac 1.0
		if { $zrad > 0 } { set fac [expr sqrt(1-(($zr*$zr)/($zrad*$zrad)))] }
		set xy_x_radius [expr $magnification * $xrad * $fac]
		set xy_y_radius [expr $magnification * $yrad * $fac]
		set xmin [expr ($lx-$xdiff-$xrad*$magnification)]
		set xmax [expr ($lx-$xdiff+$xrad*$magnification) + $magnification - 1]
		set ymin [expr ($ly-$ydiff-$yrad*$magnification)]
		set ymax [expr ($ly-$ydiff+$yrad*$magnification) + $magnification - 1]
		set cxmin [expr ($lx-$xdiff-$xy_x_radius)]
		set cxmax [expr ($lx-$xdiff+$xy_x_radius) + $magnification - 1]
		set cymin [expr ($ly-$ydiff-$xy_y_radius)]
		set cymax [expr ($ly-$ydiff+$xy_y_radius) + $magnification - 1]
		if { $id > 0 } {
			$w.frame.c create oval $cxmin $cymin $cxmax $cymax -outline $box_color -width 1 -tags [list box $idtag]
			$w.frame.c create rect $xmin $ymin $xmax $ymax -outline $box_color -width 1 -tags [list box $idtag]
			if { $show_box_labels } {
				$w.frame.c create text [expr $lx-$xdiff] [expr $ly-$ydiff] -text $id -fill $box_color -tags [list box $idtag]
			}
		} else {
			$w.frame.c create oval $cxmin $cymin $cxmax $cymax -outline $bad_color -width 1 -tags bad
		}
	}
	
	if { $nslices > 1 } {
		if { [expr $xrad*$xrad-$xr*$xr] >= 0 } {
			set yz_z_radius [expr $magnification * $zrad * sqrt(1-(($xr*$xr)/($xrad*$xrad)))]
			set yz_y_radius [expr $magnification * $yrad * sqrt(1-(($xr*$xr)/($xrad*$xrad)))]
			set xmin [expr ($lx-$zdiff-$zrad*$magnification)]
			set xmax [expr ($lx-$zdiff+$zrad*$magnification) + $magnification - 1]
			set ymin [expr ($ly-$ydiff-$yrad*$magnification)]
			set ymax [expr ($ly-$ydiff+$yrad*$magnification) + $magnification - 1]
			set cxmin [expr ($lx-$zdiff-$yz_z_radius)]
			set cxmax [expr ($lx-$zdiff+$yz_z_radius) + $magnification - 1]
			set cymin [expr ($ly-$ydiff-$yz_y_radius)]
			set cymax [expr ($ly-$ydiff+$yz_y_radius) + $magnification - 1]
			if { $id > 0 } {
				$w.frame.c1 create oval $cxmin $cymin $cxmax $cymax -outline $box_color -width 1 -tags [list box $idtag]
				$w.frame.c1 create rect $xmin $ymin $xmax $ymax -outline $box_color -width 1 -tags [list box $idtag]
				if { $show_box_labels } {
					$w.frame.c1 create text [expr $lx-$zdiff] [expr $ly-$ydiff] -text $id -fill $box_color -tags [list box $idtag]
				}
			} else {
				$w.frame.c1 create oval $cxmin $cymin $cxmax $cymax -outline $bad_color -width 1 -tags bad
			}
		}
	
		if { [expr $yrad*$yrad-$yr*$yr] >= 0 } {
			set xz_x_radius [expr $magnification * $xrad * sqrt(1-(($yr*$yr)/($yrad*$yrad)))]
			set xz_z_radius [expr $magnification * $zrad * sqrt(1-(($yr*$yr)/($yrad*$yrad)))]
			set xmin [expr ($lx-$xdiff-$xrad*$magnification)]
			set xmax [expr ($lx-$xdiff+$xrad*$magnification) + $magnification - 1]
			set ymin [expr ($ly+$zdiff-$zrad*$magnification)]
			set ymax [expr ($ly+$zdiff+$zrad*$magnification) + $magnification - 1]
			set cxmin [expr ($lx-$xdiff-$xz_x_radius)]
			set cxmax [expr ($lx-$xdiff+$xz_x_radius) + $magnification - 1]
			set cymin [expr ($ly+$zdiff-$xz_z_radius)]
			set cymax [expr ($ly+$zdiff+$xz_z_radius) + $magnification - 1]
			if { $id > 0 } {
				$w.frame.c2 create oval $cxmin $cymin $cxmax $cymax -outline $box_color -width 1 -tags [list box $idtag]
				$w.frame.c2 create rect $xmin $ymin $xmax $ymax -outline $box_color -width 1 -tags [list box $idtag]
				if { $show_box_labels } {
					$w.frame.c2 create text [expr $lx-$xdiff] [expr $ly+$zdiff] -text $id -fill $box_color -tags [list box $idtag]
				}
			} else {
				$w.frame.c2 create oval $cxmin $cymin $cxmax $cymax -outline $bad_color -width 1 -tags bad
			}
		}
	}
}

## @brief Draws a marker in the magnification window.
#
# @param	id		Marker identifier.
# @param	ptx
# @param	pty		Coordinates within the image.
# @param	ptz		Current slice (defaults to main window's slice unless animating)

proc drawMagMarker { id ptx pty ptz } {
	global theimg
	
	set w .wmag
	if ![winfo exists $w] { return }
	global project_item markers_selected marker_radius 
	global show_marker_errors show_marker_labels
	global marker_color marker_select_color error_color error_select_color

	set wc [getControlWindow $theimg]

	set magnification [$w.p.mag get]
	set magsize [$w.p.size get]
	set slice_num [$wc.slice.scale get]	
	set h [Bimage get $theimg height]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set locerr [Bmg marker $mg_item location $id]
	set x [lindex $locerr 0]
	set y [expr ($h - [lindex $locerr 1] - 1) ]
	set z [lindex $locerr 2] 
	set ex [expr $x + [lindex $locerr 3] ]
	set ey [expr $y - [lindex $locerr 4] ]
	set ez [expr $z + [lindex $locerr 5]]

	set idtag [format "m%d" $id]
	
	set pty [expr $h - $pty]
	set nslices [Bimage get $theimg nslices]
	
	set xr [expr $ptx - $x]
	set yr [expr $pty - $y]
	set zr [expr $ptz - $z]
	
	set xdiff [expr $xr*$magnification]
	set ydiff [expr $yr*$magnification]
	set zdiff [expr $zr*$magnification]
	set exdiff [expr ($ptx - $ex)*$magnification]
	set eydiff [expr ($pty - $ey)*$magnification]
	
	if { $ez > 0 } {
		set ezdiff [expr ($ptz - $ez)*$magnification]
	} else {
		set ezdiff 0
	}
	
	set lx [expr $magnification*$magsize/2]
	set ly [expr $magnification*$magsize/2]
	
#	set rad [expr $marker_radius*$magnification]
	set rad [expr $marker_radius]
	set rad2 [expr $rad*$rad]
	
	if { [lsearch $markers_selected $id] >= 0 } {
		set markcolorshow $marker_select_color;
		set errcolorshow $error_select_color;
	} else {
		set markcolorshow $marker_color;
		set errcolorshow $error_color;
	}

	if { [expr $rad2-$zr*$zr] >= 0 } { 
		set xy_radius [expr $magnification*sqrt($rad2-$zr*$zr)]
		set xmin [expr ($lx-$xdiff-$xy_radius)]
		set xmax [expr ($lx-$xdiff+$xy_radius) + $magnification - 1]
		set ymin [expr ($ly-$ydiff-$xy_radius)]
		set ymax [expr ($ly-$ydiff+$xy_radius) + $magnification - 1]
		$w.frame.c create oval $xmin $ymin $xmax $ymax -outline $markcolorshow -width 1 -tags marker 
		if { $show_marker_labels } {
			$w.frame.c create text [expr $lx-$xdiff+$rad*$magnification+10] [expr $ly-$ydiff] -text $id -fill $markcolorshow -tags [list marker $idtag]
		}
		if { $show_marker_errors } {
			$w.frame.c create line [expr $lx-$xdiff] [expr $ly-$ydiff] [expr $lx-$exdiff] [expr $ly-$eydiff] -fill $errcolorshow -width 2 -tags [list marker $idtag]
		}
	}
	
	if { $nslices > 1 } {
		if { [expr $rad2-$xr*$xr] >= 0 } {
			set yz_radius [expr $magnification*sqrt($rad2-$xr*$xr)]
			set xmin [expr ($lx-$zdiff-$yz_radius)]
			set xmax [expr ($lx-$zdiff+$yz_radius) + $magnification - 1]
			set ymin [expr ($ly-$ydiff-$yz_radius)]
			set ymax [expr ($ly-$ydiff+$yz_radius) + $magnification - 1]
			$w.frame.c1 create oval $xmin $ymin $xmax $ymax -outline $markcolorshow -width 1 -tags marker 
			if { $show_marker_labels } {
				$w.frame.c1 create text [expr $lx-$zdiff+$rad*$magnification+10] [expr $ly-$ydiff] -text $id -fill $markcolorshow -tags [list marker $idtag]
			}
			if { $show_marker_errors } {
				$w.frame.c1 create line [expr $lx-$zdiff] [expr $ly-$ydiff] [expr $lx-$ezdiff] [expr $ly-$eydiff] -fill $errcolorshow -width 2 -tags [list marker $idtag]
			}
		}
	
		if { [expr $rad2-$yr*$yr] >= 0 } {
			set xz_radius [expr $magnification*sqrt($rad2-$yr*$yr)]
			set xmin [expr ($lx-$xdiff-$xz_radius)]
			set xmax [expr ($lx-$xdiff+$xz_radius) + $magnification - 1]
			set ymin [expr ($ly+$zdiff-$xz_radius)]
			set ymax [expr ($ly+$zdiff+$xz_radius) + $magnification - 1]
			$w.frame.c2 create oval $xmin $ymin $xmax $ymax -outline $markcolorshow -width 1 -tags marker 
			if { $show_marker_labels } {
				$w.frame.c2 create text [expr $lx-$xdiff+$rad*$magnification+10] [expr $ly+$zdiff] -text $id -fill $markcolorshow -tags [list marker $idtag]
			}
			if { $show_marker_errors } {
				$w.frame.c2 create line [expr $lx-$xdiff] [expr $ly+$zdiff] [expr $lx-$exdiff] [expr $ly+$ezdiff] -fill $errcolorshow -width 2 -tags [list marker $idtag]
			}
		}
	}
}

proc drawMagSphere { x y z ptx pty ptz rad col width lbl tag } {
	global mag3D
	
	set w .wmag
	if ![winfo exists $w] { return }
	
	set c $w.frame.c
	set c1 $w.frame.c1
	set c2 $w.frame.c2
	
	set magnification [$w.p.mag get]
	set magsize [$w.p.size get]
	
	set rm [expr $rad * $magnification]
	if { $rm < 1 } { set rm 1 }
	
	set lx [expr 0.5 * $magnification * $magsize]
	set ly $lx
	
	set max [expr $lx + $rm]
	set m2 [expr $lx * $lx]
	
	set xr [expr ($ptx - $x)*$magnification]
	set yr [expr ($y - $pty)*$magnification]
	set zr [expr ($ptz - $z)*$magnification]
	
	# Test if within the magnification window
	set xr2 [expr $xr*$xr]
	if { $xr2 > $m2 } { return }

	set yr2 [expr $yr*$yr]
	if { $yr2 > $m2 } { return }

	set zr2 [expr $zr*$zr]	
	if { $zr2 > $m2 } { return }
	
	set xxd [expr $lx - $xr]
	set yyd [expr $ly - $yr]
	
	set mag_1 [expr $magnification - 1]
	
	set rm2 [expr $rm*$rm]

	set dz2 [expr $rm2 - $zr2]
	
	if { $dz2 >= 0 } {
		set xy_radius [expr sqrt($dz2)]
		set xmin [expr $xxd - $xy_radius]
		set xmax [expr $xxd + $xy_radius + $mag_1]
		set ymin [expr $yyd - $xy_radius]
		set ymax [expr $yyd + $xy_radius + $mag_1]
		$c create oval $xmin $ymin $xmax $ymax -outline $col -width $width -tags $tag 
		if { [string length $lbl] } {
			$c create text $xmax $ymax -text $lbl -fill $col -anchor w -tags $tag
		}
	}
	
	if { $mag3D } {
		set xzd [expr $lx - $zr]
		set yzd [expr $ly + $zr]
		set dx2 [expr $rm2 - $xr2]
		if { $dx2 >= 0 } {
			set yz_radius [expr sqrt($dx2)]
			set xmin [expr $xzd - $yz_radius]
			set xmax [expr $xzd + $yz_radius + $mag_1]
			set ymin [expr $yyd - $yz_radius]
			set ymax [expr $yyd + $yz_radius + $mag_1]
			$c1 create oval $xmin $ymin $xmax $ymax -outline $col -width $width -tags $tag
			if { [string length $lbl] } {
				$c1 create text $xmax $ymax -text $lbl -fill $col -anchor w -tags $tag
			}
		}
	
		set dy2 [expr $rm2 - $yr2]
		if { $dy2 >= 0 } {
			set xz_radius [expr sqrt($dy2)]
			set xmin [expr $xxd - $xz_radius]
			set xmax [expr $xxd + $xz_radius + $mag_1]
			set ymin [expr $yzd - $xz_radius]
			set ymax [expr $yzd + $xz_radius + $mag_1]
			$c2 create oval $xmin $ymin $xmax $ymax -outline $col -width $width -tags $tag
			if { [string length $lbl] } {
				$c2 create text $xmax $ymax -text $lbl -fill $col -anchor w -tags $tag
			}
		}
	}
}

proc drawMagLine { xs ys zs xe ye ze ptx pty ptz width col tag } {
	global mag3D

	set w .wmag
	if ![winfo exists $w] { return }
	
	set c $w.frame.c
	set c1 $w.frame.c1
	set c2 $w.frame.c2
	
	set magnification [$w.p.mag get]
	set magsize [$w.p.size get]

	set width [expr int($width * $magnification)]
	if { $width < 1 } { set width 1 }
	
	set xsr [expr $ptx - $xs]
	set ysr [expr $ys - $pty]
	set zsr [expr $ptz - $zs]
	set xer [expr $ptx - $xe]
	set yer [expr $ye - $pty]
	set zer [expr $ptz - $ze]
	
	set xsdiff [expr $xsr*$magnification]
	set ysdiff [expr $ysr*$magnification]
	set zsdiff [expr $zsr*$magnification]
	set xediff [expr $xer*$magnification]
	set yediff [expr $yer*$magnification]
	set zediff [expr $zer*$magnification]
	
	set magdim [expr $magnification*$magsize]
	set lx [expr $magdim * 0.5]
	set ly [expr $magdim * 0.5]

	set z -1.0
	if { $zer != $zsr } {
		set z [expr $zsr / (1.0*($ze - $zs))]
	} else {
		if { abs($zsr) < 0.5 } { set z 0.5 }
	}
	
	set xmin [expr $lx - $xsdiff]
	set xmax [expr $lx - $xediff]
	set ymin [expr $ly - $ysdiff]
	set ymax [expr $ly - $yediff]
	if { $xmin == $xmax && $ymin == $ymax } { set ymax [expr $ymax + 1] }
	
	if { $z >= 0.0 && $z <= 1.0 } { 
		$c create line $xmin $ymin $xmax $ymax -fill $col -smooth yes -width $width -tags $tag
	}
	
	if { $mag3D } {
		set x -1.0
		if { $xe != $xs } {
			set x [expr $xsr / (1.0*($xe - $xs))]
		} else {
			if { abs($xsr) < 0.5 } { set x 0.5 }
		}
	
		set zmin [expr $ly - $zsdiff]
		set zmax [expr $ly - $zediff]
		if { $x >= 0.0 && $x <= 1.0 } {
			if { $zmin == $zmax && $ymin == $ymax } { set ymax [expr $ymax + 1] }
			$c1 create line $zmin $ymin $zmax $ymax -fill $col -smooth yes -width $width -tags $tag
		}
	
		set y -1.0
		if { $ye != $ys } {
			set y [expr $ysr / (1.0*($ys - $ye))]
		} else {
			if { abs($ysr) < 0.5 } { set y 0.5 }
		}

		set zmin [expr $magdim - $zmin]
		set zmax [expr $magdim - $zmax]
		if { $y >= 0.0 && $y <= 1.0 } {
			if { $xmin == $xmax && $zmin == $zmax } { set zmax [expr $zmax + 1] }
			$c2 create line $xmin $zmin $xmax $zmax -fill $col -smooth yes -width $width -tags $tag
		}
	}
}

proc changeSlice { deltaslice } {
	global maglock theimg
	
	set w .wmag
	
	if { $maglock == 1 } {
		updateMagnify [$w.n.x get] [$w.n.y get] $deltaslice slicer
	} else {
		set wc [getControlWindow $theimg]
		$wc.slice.scale set [$w.q.slice get]
	}
}

proc updateSlicer { aniaxis anidirection } {
	set w .wmag
	
	set x [$w.n.x get]
	set y [$w.n.y get]
	set z [$w.q.slice get]
	if { $aniaxis == "x" } {
		if { $anidirection == 1 } { updateMagnify [expr $x + 1] $y $z slicer }
		if { $anidirection == 0 } { updateMagnify [expr $x - 1] $y $z slicer }
	}
	if { $aniaxis == "y" } {
		if { $anidirection == 1 } { updateMagnify $x [expr $y + 1] $z slicer }
		if { $anidirection == 0 } { updateMagnify $x [expr $y - 1] $z slicer }
	}
	if { $aniaxis == "z" } {
		if { $anidirection == 1 } {$w.q.slice set [expr [$w.q.slice get] + 1] }
		if { $anidirection == 0 } {$w.q.slice set [expr [$w.q.slice get] - 1] }
	}
}

proc animateSlices {} {
	global theimg
	global aniproc2
	global aniaxis
	
	set w .wmag
	
	if ![winfo exists $w] { return }
	if { $aniaxis=="x" } {
		set anisource "$w.n.x"
		set lim [expr [Bimage get $theimg width] - 1]
	}
	if { $aniaxis=="y" } {
		set anisource "$w.n.y"
		set lim [expr [Bimage get $theimg height] - 1]
	}
	if { $aniaxis=="z" } {
		set anisource "$w.q.slice"
		set lim [expr [Bimage get $theimg nslices] - 1]
	}
	set aniceiling [expr [$anisource get] + [$w.r.range get]]
	set anifloor [expr [$anisource get] - [$w.r.range get]]
	if { $aniceiling > $lim } { set aniceiling $lim }
	if { $anifloor < 0 } { set anifloor 0 }
	set anidirection 1
	set speed [$w.r.speed get]
	if { $speed == "" || $speed < 0 } {
		set speed 15
	}
	if { $speed > 200 } {
		set speed 200
	}
	$w.r.speed delete 0 end
	$w.r.speed insert 0 $speed

	if { $aniproc2 == 0 } { 
		$w.r.animate configure -text "Stop"
		$w.r.axis configure -state disabled
		$w.r.range configure -state disabled
		$w.r.speed configure -state disabled
		set aniproc2 1
	} else {
		$w.r.animate configure -text "Animate"
		$w.r.axis configure -state normal
		$w.r.range configure -state normal
		$w.r.speed configure -state normal
		set aniproc2 0
	}
	
	while { [winfo exists $w] && $aniproc2 == 1 } {
		set slice [$anisource get]
		if { $slice == $aniceiling } {
			set anidirection 0
		}
		if { $slice == $anifloor } {
			set anidirection 1
		}
		after [expr 1000/[$w.r.speed get]]
		updateSlicer $aniaxis $anidirection
		update
	}
}

proc check_length {newstring max_length} {
	return [expr {[string length $newstring] <= $max_length}]
}

