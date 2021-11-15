##
# @file		bshow_ctf.tcl
#
# @brief	Procedures to do CTF fitting
#
# @author	Bernard Heymann
# @date		Created: 20020729
# @date		Modified: 20170712

set plot_width 500
set plot_height 300
set ctf_win_x 768
set ctf_win_y 768
set ctf_scale_x 1
set ctf_scale_y 0
set cont_update 0
set show_rings 1
set base_type 1
set env_type 4
set sub_base 0
set rps ""
set base ""
set envl ""
set ctf ""
set ctfdir [pwd]
set volt 120.0
set Cs 2.0
set amp_fac 0.07
set ctf_hires 0.5
set ctf_lores 50

## @brief Dialog box to manipulate CTF parameters
#
# @param	currimg		Image identifier.

proc CTF { currimg } {
	global project_item theimg imgtype
	global helv12 step
	global plot_width plot_height
	global ctf_win_x ctf_win_y
	global ctf_scale_x ctf_scale_y
    global volt Cs amp_fac base_type env_type ctf_hires ctf_lores

	if ![Bimage exists $theimg] {
		tk_dialog .dialog "No image in memory!" \
				{Please read in an appropriate image first.} \
				info 0 {OK}
		return
	}
	
	if { [string first "ps" $imgtype] < 0 } {
		set i [tk_dialog .dialog "No power spectrum!" \
			{No power spectrum was found. Please transform the image first to a power spectrum.} \
			info 0 {Cancel} {Proceed anyway}]
		if { $i < 1 } { return }
		set imgtype "ps"
	}
	set theimg $currimg
#	puts "canvas = $c  image = $theimg"
	set filename [Bimage get $theimg filename]
	set width [Bimage get $theimg width]
	set height [Bimage get $theimg height]
	setStep $theimg
	
	if ![Bmg exists] { createMgParam }
	
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]

#	puts "project_item: $mg_item  file: $filename  imgtype: $imgtype"
#	puts [Bmg get $mg_item filename ps]
	if ![string length [Bmg get $mg_item filename ps]] {
		Bmg set $mg_item filename $filename ps
		puts "Power spectrum: $filename"
	}

	setCTFPixelSize $theimg
	
	set w .wctf
	if ![winfo exists $w] {
		toplevel $w
		wm geometry $w [format "%dx%d" $ctf_win_x $ctf_win_y]

		if { $ctf_scale_y < 0.1 } {
#			set range [expr [.max.scale get] - [.min.scale get]]
			set range [expr [Bimage get $theimg max] - [Bimage get $theimg min]]
			if { $range <= 0 } { set range $plot_height }
			set ctf_scale_y [expr int($plot_height/$range + 1)]
		}
		
#		set hires [expr int(200 * [.pixel_size.entry get]) / 100.0]
		set ps [Bimage get $theimg pixel_size $n]
		set hires [expr int(200 * [lindex $ps 0]) / 100.0]
		if { $hires > $ctf_hires } { set ctf_hires $hires }
		if { $ctf_lores < [expr 2 * $ctf_hires] } {
			set ctf_lores [expr 2 * $ctf_hires]
		}

		menu $w.menuBar -tearoff 0
		$w.menuBar add cascade -menu $w.menuBar.box -label "CTF" -underline 0
		menu $w.menuBar.box -tearoff 0
		$w configure -menu $w.menuBar
		$w.menuBar.box add command -label "Help" -underline 0 \
				-command { showURL "file://$Bsoft/doc/bshow/bshow_ctf.html" }
		$w.menuBar.box add command -label "Write curves" \
				-command { writeCTF } -underline 0
		$w.menuBar.box add command -label "Save window" \
				-command [list capture $w png t.png] -underline 0
#				-command { canvas2Photo $w } -underline 0
		$w.menuBar.box add command -label "Close" \
				-command "destroy $w" -underline 0
		
		set g $w.plot
		canvas $g -width $plot_width -height $plot_height -background white \
				-relief sunken -borderwidth 2

		setupItemID $w.itemid

		# Plot parameters
		frame $w.pp -bg cyan
		label $w.pp.origin_tag -font $helv12 -text "Origin" -bg cyan
		entry $w.pp.ox_entry -width 6 -validate key -vcmd { string is double %P }
		$w.pp.ox_entry insert 0 20
		entry $w.pp.oy_entry -width 6 -validate key -vcmd { string is double %P }
		$w.pp.oy_entry insert 0 20
		label $w.pp.scale_tag -font $helv12 -text "Scale" -bg cyan
		entry $w.pp.sx_entry -width 6 -validate key -vcmd { string is double %P }
		$w.pp.sx_entry insert 0 $ctf_scale_x
		entry $w.pp.sy_entry -width 6 -validate key -vcmd { string is double %P }
		$w.pp.sy_entry insert 0 $ctf_scale_y
		label $w.pp.res_tag -font $helv12 -text "Resolution" -bg cyan
		entry $w.pp.lores_entry -width 6 -validate key -vcmd { string is double %P }
		$w.pp.lores_entry insert 0 $ctf_lores
		entry $w.pp.hires_entry -width 6 -validate key -vcmd { string is double %P }
		$w.pp.hires_entry insert 0 $ctf_hires
		label $w.pp.val_tag -font $helv12 -text "Value" -bg cyan
		label $w.pp.value -width 6 -text "" -relief sunken -bd 1 -anchor w

		# Placement in frames
		pack $w.pp.origin_tag $w.pp.ox_entry $w.pp.oy_entry \
				$w.pp.scale_tag $w.pp.sx_entry $w.pp.sy_entry \
				$w.pp.res_tag $w.pp.lores_entry $w.pp.hires_entry \
				$w.pp.value \
				-side left -ipadx 3 -ipady 2
		
		# These scales are used for CTF parameters
		scale $w.def_avg_scale -orient horizontal -from 0.02 -to 20 \
				-command { setDefocusAverage } -label "Defocus average (um)"\
				-resolution 0.001
		$w.def_avg_scale set 2

		scale $w.def_dev_scale -orient horizontal -from 0 -to 5 \
				-command { setDefocusDeviation } -label "Defocus deviation (um)"\
				-resolution 0.001
		$w.def_dev_scale set 0

		scale $w.ast_ang_scale -orient horizontal -from -90 -to 90 \
				-command { setAstigmatismAngle } -label "Astigmatism angle (degrees)"
		$w.ast_ang_scale set 0

		# Imaging parameters
		frame $w.vca -bg cyan
		label $w.vca.volt_tag -font $helv12 -text "Voltage (kV)" -bg cyan
		entry $w.vca.volt_entry -width 6 -validate key -vcmd { string is double %P }
		$w.vca.volt_entry insert 0 $volt
		label $w.vca.cs_tag -font $helv12 -text "Cs (mm)" -bg cyan
		entry $w.vca.cs_entry -width 6 -validate key -vcmd { string is double %P }
		$w.vca.cs_entry insert 0 $Cs
		label $w.vca.amp_tag -font $helv12 -text "Amp" -bg cyan
		entry $w.vca.amp_entry -width 6 -validate key -vcmd { string is double %P }
		$w.vca.amp_entry insert 0 $amp_fac
		pack $w.vca.volt_tag $w.vca.volt_entry $w.vca.cs_tag $w.vca.cs_entry \
				$w.vca.amp_tag $w.vca.amp_entry -side left -ipadx 3

		frame $w.base_eq -bg green
		label $w.base_eq.tag -font $helv12 -text "Baseline" -bg green
		entry $w.base_eq.entry
		$w.base_eq.entry insert 0 {0.2*exp(-5000*$s2) + 0.14*exp(-300*$s2) + 0.4}

		frame $w.base_buttons -bg green
		label $w.base_buttons.tag -font $helv12 -text "        Type" -bg green
		tk_optionMenu $w.base_buttons.type base_type "1" "2" "3" "4" "5" "6"
		checkbutton $w.base_buttons.subtract -text "Subtract" \
				-variable sub_base -relief raised -command "drawCTF" -bg green

		# Placement in frames
		pack $w.base_eq.tag -side left
		pack $w.base_eq.entry -side left -fill x -expand yes -ipady 2
		pack $w.base_buttons.tag $w.base_buttons.type -side left -ipadx 2
		pack $w.base_buttons.subtract -side left

		frame $w.env_eq -bg orange
		label $w.env_eq.tag -font $helv12 -text "Envelope" -bg orange
		entry $w.env_eq.entry
		$w.env_eq.entry insert 0 {0.2*exp(-1000*$s2)}
		
		frame $w.env_buttons -bg orange
		label $w.env_buttons.tag -font $helv12 -text "        Type" -bg orange
		tk_optionMenu $w.env_buttons.type env_type "1" "2" "3" "4"

		# Placement in frames
		pack $w.env_eq.tag -side left
		pack $w.env_eq.entry -side left -fill x -expand yes -ipady 2
		pack $w.env_buttons.tag $w.env_buttons.type -side left -ipadx 2

		# Buttons to do CTF fits
		frame $w.fit -bg cyan
		label $w.fit.tag -text "Fit " -bg cyan
		button $w.fit.quick -text "Quick" \
				-relief raised -command "autoCTFfit 0"
		button $w.fit.baseline -text "Baseline" \
				-relief raised -command "autoCTFfit 1"
		button $w.fit.envelope -text "Envelope" \
				-relief raised -command "autoCTFfit 2"
		button $w.fit.defocus -text "Defocus" \
				-relief raised -command "autoCTFfit 3"
		button $w.fit.astig -text "Astigmatism" \
				-relief raised -command "autoCTFfit 4"
		pack $w.fit.tag $w.fit.quick $w.fit.baseline $w.fit.envelope \
				$w.fit.defocus $w.fit.astig -side left -ipadx 5
		
		# Defocus range for quick fit
		frame $w.def -bg cyan
		label $w.def.tag -text "Defocus:" -bg cyan
		label $w.def.start_tag -text "min" -bg cyan
		entry $w.def.start -width 4
		label $w.def.end_tag -text "max" -bg cyan
		entry $w.def.end -width 4
		label $w.def.inc_tag -text "inc" -bg cyan
		entry $w.def.inc -width 4
		$w.def.start insert 0 0.1
		$w.def.end insert 0 20
		$w.def.inc insert 0 0.1
		pack $w.def.tag $w.def.start_tag $w.def.start $w.def.end_tag \
				$w.def.end $w.def.inc_tag $w.def.inc -side left -ipadx 5
		
		# Checkbutton to update the plot and show rings
		frame $w.checks
		checkbutton $w.checks.update -text "updateCTF" -variable cont_update \
				-command "updateCTF"
		checkbutton $w.checks.rings -text "Show rings" -variable show_rings \
				-command [list drawEllipse $theimg]
		pack $w.checks.update $w.checks.rings -side left -padx 2
		setupMicrographButtons $w.checks "micrograph"
				
		# Placement in window
		pack $w.itemid -side top -fill x -ipady 1
		pack $w.plot -side top -fill both -expand yes
		pack $w.pp $w.def_avg_scale $w.def_dev_scale $w.ast_ang_scale -side top -fill x -ipady 1
		pack $w.vca $w.base_eq $w.base_buttons $w.env_eq $w.env_buttons -side top -fill x -ipady 1
		pack $w.checks -side bottom -fill x -ipady 1
		pack $w.fit $w.def -side left -fill x -ipady 1
	} else {
		wm deiconify $w
		raise $w
    }
	
	getCTFparam
	wm title $w [list "CTF:" $filename "(" $n ")"]
	wm iconname $w "CTF"

	bind $g <Motion> "showFrequency %W %x %y"
	bind $w <Return> "updateCTF"

	bind $w <Control-w> { destroy $w }
	bind $w <Command-w> { destroy $w }
	.menuBar.window entryconfigure "CTF" -state normal
}

## @brief Dialog box to set the pixel size for CTF fitting
#
# @param	theimg		Image identifier.

proc setCTFPixelSize { theimg } {
	global helv12
#	set wc [getControlWindow $theimg]
	set ps [getPixelSizeEntry]
#	puts "pixel size = $pixel_size"
	if { [lindex $ps 0] < 100 } { return }
	set w .wpx
	if ![winfo exists $w] {
		toplevel $w
		
#		frame $w.pxs
#		label $w.pxs.tag -font $helv12 -text "Set a reasonable pixel size"
#		entry $w.pxs.entry -width 6 -validate key -vcmd { string is double %P }
#		$w.pxs.entry insert 0 $pixel_size

		setupVectorEntry $w.pxs "Pixel size (Å)" double $ps "Image pixel size in Å/pixel"

		frame $w.buttons
		button $w.buttons.ok -text OK -command "setPixelSize_and_DestroyWindow $w $theimg"
		button $w.buttons.cancel -text Cancel -command "destroy $w"
		
		pack $w.pxs.tag $w.pxs.entry -side left -pady 5
		pack $w.buttons.ok $w.buttons.cancel -side left -expand 1
		pack $w.pxs $w.buttons -side top -fill x -pady 2m
	
		bind $w <Return> "setPixelSize_and_DestroyWindow $w $theimg"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Pixel Size"
    wm iconname $w "Pixel Size"
	tkwait window $w
}

## @brief Sets the pixel size.
#
# @param	w 			the window to destroy.
# @param	theimg		Image identifier.

proc setPixelSize_and_DestroyWindow { w theimg } {
	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set ps [$w.pxs.entry get]
	destroy $w
	Bimage set $theimg pixel_size $img_num [lindex $ps 0] [lindex $ps 1] [lindex $ps 2]
	updatePixelSizeEntry
	setStep $theimg
}

## @brief Sets the reciprocal space step size.
#
# @param	theimg		Image identifier.

proc setStep { theimg } {
	global step
	set pixel_size [Bimage get $theimg pixel_size]
	set width [Bimage get $theimg width]
	set step [expr 1.0 / ([lindex $pixel_size 0] * $width)]
}

## @brief Updates all the CTF curves
#

proc updateCTF { } {
	global theimg project_item
	global step
	global ctf_win_x ctf_win_y
	set w .wctf
	if [winfo exists $w] {
		updateItemID $w.itemid
		setStep $theimg
		drawEllipse $theimg
		drawCTF
		set ctf_win_x [winfo width $w]
		set ctf_win_y [winfo height $w]
	} else {
		set c [getImageCanvas $theimg]
		$c delete ellipse
	}
}

## @brief Sets the defocus average
#
# @param	def_avg 	 	Defocus average.

proc setDefocusAverage { def_avg } {
	global project_item theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set w .wctf
	Bmg set $mg_item defocus [expr 1e4 * [$w.def_avg_scale get]]
	updateCTF
}

## @brief Sets the defocus deviation
#
# @param	def_dev 		Defocus deviation.

proc setDefocusDeviation { def_dev } {
	global project_item theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set w .wctf
	Bmg set $mg_item defocus_deviation [expr 1e4 * [$w.def_dev_scale get]]
	updateCTF
}

## @brief Sets the astigmatism angle
#
# @param	ast_ang	 	Astigmatism angle in degrees.

proc setAstigmatismAngle { ast_ang } {
	global project_item theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set w .wctf
	Bmg set $mg_item astigmatism_angle [$w.ast_ang_scale get]
	updateCTF
}

## @brief Calculates the wavelength
#

proc getWavelength {} {
	set volt [expr 1000 * [.wctf.vca.volt_entry get]]
	return [expr 12.26 / sqrt($volt * (1 + $volt * 0.9788e-6))]
}

## @brief Returns the first zero at the average defocus

proc getFirstZero { } {
	global project_item theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
#	puts "Getting the first zero"
	set firstzero [Bmg get $mg_item zero]
	if !$firstzero { set firstzero 0 }
	return $firstzero
}

## @brief Returns all the zeroes at the average defocus

proc getZeroes { } {
	global project_item theimg
	set wc [getControlWindow $theimg]
	set Cs [expr 1e7 * [.wctf.vca.cs_entry get]]
	set lambda [getWavelength]
#	puts "Wavelength $lambda angstrom"
	set invl2Cs [expr 1.0 / ($lambda * $lambda * $Cs)]
	set inv2l3Cs [expr 2.0 / ($lambda * $lambda * $lambda * $Cs)]
	set def_avg [expr 1e4 * [.wctf.def_avg_scale get]]
	set ctf_fz [expr $def_avg * $invl2Cs]
#	set max_s [expr 0.5 / [.pixel_size.entry get]]
	set max_s [expr 0.5 / [$wc.pixel_size.x get]]
	set nz [expr lambda*def_avg*max_s*max_s - $lambda*$lambda*$lambda*$Cs*max_s*max_s*max_s*max_s/2]
	if { $nz < 1 } { set nz 1 }
	set zeroes ""
	for {set i 0} {$i < $nz} {incr i 1} {
		set nt [expr ($i+1)*$inv2l3Cs]
		set ctf_fz2 [expr $ctf_fz * $ctf_fz]
		if { $ctf_fz2 > $nt } {
			append zeroes " " [expr sqrt(ctf_fz - sqrt(ctf_fz2 - $nt))]
		}
	}
	return $zeroes
}

## @brief Calculates the radial power spectrum
#

proc calcRPS {  } {
	global project_item theimg
	global rps wri
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set rps [Bmg rps $mg_item]
	if { [lsearch $rps "Radial"] > -1 } {
		tk_messageBox -message $rps \
			-icon question -type ok \
			-detail "The defocus or Cs values may be invalid!"
		set rps ""
	}
	if { [llength $rps] < 1 } {
		tk_messageBox -message "The radial average was not calculated!" \
			-icon question -type ok \
			-detail "Project item $mg_item"
	}
	set wri [lindex $rps end]
	set rps [lreplace $rps end end]
}

## @brief Calculates the contrast transfer function
#

proc calcCTF { } {
	global PI step
	global base envl
	global rps ctf
	set Cs [expr 1e7 * [.wctf.vca.cs_entry get]]
	set ampfac [.wctf.vca.amp_entry get]
	set phifac [expr sqrt(1 - $ampfac * $ampfac)]
	set base_eq [.wctf.base_eq.entry get]
	set env_eq [.wctf.env_eq.entry get]
	set def_avg [expr 1e4 * [.wctf.def_avg_scale get]]
	set lambda [getWavelength]
	set hpil3Cs [expr 0.5 * $PI * $lambda * $lambda * $lambda * $Cs];
	set pil [expr $PI * $lambda];
	set n [llength $rps]
	set base ""
	set envl ""
	set ctf ""
#	puts $base_eq
#	puts $env_eq
	for {set i 0} {$i < $n} {incr i 1} {
		set s [expr $i * $step]
		set s2 [expr $s * $s]
		set s3 [expr $s * $s2]
		set s4 [expr $s2 * $s2]
	    set dphi [expr $hpil3Cs*$s4 - $pil*$def_avg*$s2]
		set amp [expr $ampfac*cos($dphi) - $phifac*sin($dphi)]
		set e [expr $env_eq]
		set b [expr $base_eq]
		set c [expr $amp * $amp]
		append base " " $b
		append envl " " $e
		append ctf " " $c
	}
}

## @brief Returns the RMSD between the experimental and calculated RPS curves.

proc calcR { } {
	global theimg
	global step base envl
	global ctf rps
	global ctf_hires ctf_lores
	set wc [getControlWindow $theimg]
	set pixel_size [$wc.pixel_size.x get]
	set lores [.wctf.pp.lores_entry get]
	set hires [.wctf.pp.hires_entry get]
	if { $hires < [expr 2.1 * $pixel_size] } {
		set hires [expr 2.1 * $pixel_size]
		.wctf.pp.hires_entry delete 0 end
		.wctf.pp.hires_entry insert 0 $hires
		set ctf_hires $hires
	}
	if { $lores < $hires } {
		set lores [expr 2 * $hires]
		.wctf.pp.lores_entry delete 0 end
		.wctf.pp.lores_entry insert 0 $lores
		set ctf_lores $lores
	}
	set nr [llength $rps]
	set ilo [expr int(1.0 / ($lores*$step))]
	set ihi [expr int(1.0 / ($hires*$step))]
	if { $ihi <= $ilo } { set ihi [expr $ilo + 1] }
	if { $ilo < 1 } { set ilo 1 }
	if { $ihi >= $nr } { set ihi [expr $nr - 1] }
	set n 0.0
	set R 0.0
	for {set i $ilo} {$i <= $ihi} {incr i 1} {
		set c [expr [lindex $base $i] + [lindex $envl $i] * [lindex $ctf $i]]
		set d [expr [lindex $rps $i] - $c]
		set w 1.0
#		set w [expr 1.0 - [lindex $ctf $i]]
		if { $w > 0 } {
			set R [expr $R + $d * $d * $w]
			set n [expr $n + $w]
		}
	}
	if { $n > 0 } {
		set R [expr sqrt($R/$n)]
	} else {
		set R 1000.0
	}
	return $R
}

## @brief Attempts an automatic fit of the baseline, envelope and CTF functions

proc autoCTFfit { level } {
	global project_item theimg
	global cont_update base_type env_type
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set w .wctf
	setCTFparam
	set lo_res [$w.pp.lores_entry get]
	set hi_res [$w.pp.hires_entry get]
	set def_start [expr [$w.def.start get] * 1e4]
	set def_end [expr [$w.def.end get] * 1e4]
	set def_inc [expr [$w.def.inc get] * 1e4]
#	puts "Fitting the CTF"
#	puts "base_type=$base_type"
	set ps [Bmg get $mg_item pixel_size]
	if { [lindex $ps 0] > 1.5 } {	# Make sure we are not trying to fit the water ring
		if { $base_type > 3 } {
			set base_type [expr $base_type - 3]
		}
	}
	set old_env_type [Bmg get $mg_item envelope_type]
	Bmg set $mg_item baseline_type $base_type
	Bmg set $mg_item envelope_type $env_type
	set result [Bmg ctf_fit $mg_item $level $lo_res $hi_res $def_start $def_end $def_inc]
	if { $level == 1 && $env_type != $old_env_type } {
		set result [Bmg ctf_fit $mg_item 2 $lo_res $hi_res]
	}
#	puts $result
	if { [string length $result] < 2 } { return }
	$w.base_eq.entry delete 0 end
	$w.base_eq.entry insert 0 [lindex $result 0]
	$w.env_eq.entry delete 0 end
	$w.env_eq.entry insert 0 [lindex $result 1]
	$w.def_avg_scale set [expr 1e-4 * [lindex $result 2]]
	$w.def_dev_scale set [expr 1e-4 * [lindex $result 3]]
	$w.ast_ang_scale set [lindex $result 4]
	if { $level == 4 } { set cont_update 1 }
	updateCTF
}

## @brief Draws a number of ellipses at the CTF zeroes.
#
# @param	theimg		Image identifier.

proc drawEllipse { theimg } {
	global PI show_rings step
	set c [getImageCanvas $theimg]
	$c delete ellipse
	if { $show_rings < 1 } { return }
	set wc [getControlWindow $theimg]
	set pixel_size [$wc.pixel_size.x get]
	set Cs [expr 1e7 * [.wctf.vca.cs_entry get]]
	set phi [expr asin([.wctf.vca.amp_entry get]) / $PI]
	set lambda [getWavelength]
#	puts "Wavelength $lambda angstrom"
	set invl2Cs [expr 1.0 / ($lambda * $lambda * $Cs)]
	set inv2l3Cs [expr 2.0 / ($lambda * $lambda * $lambda * $Cs)]
	set width [Bimage get $theimg width]
	set height [Bimage get $theimg height]
	set scale [$wc.scale.scale get]
#	set originx [expr $scale * ($width + 1) / 2]
#	set originy [expr $scale * ($height + 2) / 2]
	set originx [expr $scale * $width / 2]
	set originy [expr $scale * $height / 2]
	set def_avg [expr 1e4 * [.wctf.def_avg_scale get]]
	set def_dev [expr 1e4 * [.wctf.def_dev_scale get]]
	set ast_ang [expr [.wctf.ast_ang_scale get] * $PI / 180]
	if { $def_dev > [expr 0.99 * $def_avg] } {
		set def_dev [expr 0.99 * $def_avg]
		.wctf.def_dev_scale set [expr $def_dev * 1e-4]
	}
	set ellip2 [expr ( $def_avg + $def_dev ) / ( $def_avg - $def_dev ) ]
	set fz [expr $def_avg * $invl2Cs]
	set fz2 [expr $fz * $fz]
	set nmax [expr int($lambda * $def_avg / (4 * $pixel_size * $pixel_size))]
#	puts $nmax
	if { $nmax > 10 } { set nmax 10 }
	set astart [expr $PI/2]
	set aend [expr $astart + $show_rings * $PI + 0.05]
	set astep [expr $PI/36]
	for { set n 1 } { $n < $nmax } { incr n 1 } {
		set r $width
		set yscale [expr $height*1.0/$width]
		set t [expr ($n - $phi) * $inv2l3Cs]
		if { $fz2 > $t } {
			set r [expr ($scale / $step) * sqrt($fz - sqrt($fz2 - $t))]
		}
#		puts "Pixel distance $r"
		set a $astart
		set da [expr $a - $ast_ang]
		set re [expr $r * sqrt( ($ellip2 + 1) / \
			(2 * ($ellip2 * cos($da) * cos($da) + sin($da) * sin($da))) ) ]
		set xl [expr $re*cos($a) + $originx]
		set yl [expr $re*$yscale*sin($a) + $originy]
		set yl [expr $scale * $height - $yl]
		for {set a $astart} {$a <= $aend} {set a [expr $a + $astep]} {
			set da [expr $a - $ast_ang]
			set re [expr $r * sqrt( ($ellip2 + 1) / \
				(2 * ($ellip2 * cos($da) * cos($da) + sin($da) * sin($da))) ) ]
			set x [expr $re*cos($a) + $originx]
			set y [expr $re*$yscale*sin($a) + $originy]
			# y is flipped because the canvas origin is upper left
			set y [expr $scale * $height - $y]
			$c create line $xl $yl $x $y -dash 2 -fill yellow -smooth yes -width 1 -tags ellipse
			set xl $x
			set yl $y
		}
	}
}

## @brief Draws four curves: baseline, envelope, CTF and RPS
#

proc drawCTF { } {
	global cont_update sub_base
	global bl step
	global base envl ctf rps wri
	global plot_height
	global ctf_scale_x ctf_scale_y
	if { $cont_update || [string length $rps] < 2 } {
		calcRPS
	}
	set c .wctf.plot
	set R [calcCTF]
	set base_line ""
	set env_line ""
	set ctf_line ""
	set rps_line ""
	$c delete xscale
	$c delete ctf
	$c delete rps
	set orix [.wctf.pp.ox_entry get]
	set oriy [expr $plot_height - [.wctf.pp.oy_entry get]]
	set ctf_scale_x [.wctf.pp.sx_entry get]
	set ctf_scale_y [.wctf.pp.sy_entry get]
#	set wri [lindex $rps end]
#	set rps [lreplace $rps end end]
	set n [llength $rps]
	for {set i 1} {$i < $n} {incr i 1} {
		set x [expr $orix + $ctf_scale_x * $i]
		if { [string length $rps] > 2 } {
			set y [lindex $rps $i]
		} else {
			set y 0
		}
		set b [lindex $base $i]
		set e [lindex $envl $i]
		set cv [lindex $ctf $i]
		if { $sub_base } { 
			set y [expr $y - $b]
			set cv [expr $e * $cv]
			set b 0
		} else { 
			set cv [expr $b + $e * $cv]
			set e [expr $b + $e]
		}
		if { [string length $rps] > 2 } {
			set y [expr $oriy - $ctf_scale_y * $y]
			append rps_line " " $x " " $y
		}
		set b [expr $oriy - $ctf_scale_y * $b]
		set e [expr $oriy - $ctf_scale_y * $e]
		set cv [expr $oriy - $ctf_scale_y * $cv]
		append base_line " " $x " " $b
		append env_line " " $x " " $e
		append ctf_line " " $x " " $cv
	}
	$c create line $base_line -fill green -smooth yes -width 1 -tags ctf
	$c create line $env_line -fill orange -smooth yes -width 1 -tags ctf
	$c create line $ctf_line -fill red -smooth yes -width 1 -tags ctf
	if { [string length $rps_line] > 2 } {
		$c create line $rps_line -fill blue -smooth yes -width 1.5 -tags rps
	}
	$c create line $orix 0 $orix [expr 2 * $oriy] \
		-fill black -smooth yes -width 1 -tags xscale
	$c create line $orix $oriy [expr $orix + $ctf_scale_x * $n] $oriy \
		-fill black -smooth yes -width 1 -tags xscale
	set res 0
	for {set i 0} {$i < $n} {incr i 50} {
		set x [expr $orix + $ctf_scale_x * $i]
		if { $i } { set res [format "%5.1f" [expr 1.0 /($i*$step)]] }
		$c create text $x [expr $oriy + 10] -text $res -fill black -tags xscale
		$c create line $x $oriy $x [expr $oriy - 5] -fill black -smooth yes -width 1 -tags xscale
	}
	set lores [.wctf.pp.lores_entry get]
	set x [expr $orix + $ctf_scale_x/($lores*$step)]
#	puts "$oriy $plot_height"
	$c create line $x 0 $x $plot_height \
		-fill green -smooth yes -width 1 -tags xscale
	set hires [.wctf.pp.hires_entry get]
	set x [expr $orix + $ctf_scale_x/($hires*$step)]
	$c create line $x 0 $x $plot_height \
		-fill green -smooth yes -width 1 -tags xscale
	set R [calcR]
#	set Z [expr 1.0 / [getFirstZero]]
	set Z [getFirstZero]
	$c create text [expr $orix + 300] 30 -text [format "R = %10.6f" $R] -anchor w -fill black -tags rps
	$c create text [expr $orix + 300] 50 -text [format "First zero = %10.2f A" $Z] -anchor w -fill black -tags rps
	$c create text [expr $orix + 300] 70 -text [format "Water ring index = %8.4f" $wri] -anchor w -fill black -tags rps
}

## @brief Shows the spatial frequency at that point where the mouse is moved.
#
# @param	c			Canvas object.
# @param	x
# @param	y			Coordinates on the canvas.

proc showFrequency { c x y } {
	global step
	global ctf_scale_x
	set orix [.wctf.pp.ox_entry get]
	set r [expr $ctf_scale_x * 1.0 / (($x - $orix)*$step)]
	set w .wctf
	$w.pp.value config -text [format "%6.2f" $r]
}

## @brief Saves CTF curves to a text file
#

proc writeCTF { } {
	global theimg
	global step base envl
	global ctf rps
	if { [string length $rps] < 2 } { return }
	set filename [Bimage get $theimg filename]
	set ctffile "[file rootname [file tail $filename]].txt"
	set ctffile [tk_getSaveFile -initialfile $ctffile -defaultextension .txt]
	if { [string length $ctffile] < 1 } { return }
	set fctf [open $ctffile w]
	set wc [getControlWindow $theimg]
	set pixel_size [$wc.pixel_size.x get]
	set n [llength $rps]
	set w .wctf
	puts $fctf "File name:              $filename"
	puts $fctf "Pixel size:             $pixel_size angstrom/pixel"
	puts $fctf "Acceleration voltage:   [$w.vca.volt_entry get] kV"
	puts $fctf "Cs:                     [$w.vca.cs_entry get] mm"
	puts $fctf "Amplitude contribution: [$w.vca.amp_entry get]"
	puts $fctf "Defocus average:        [$w.def_avg_scale get] um"
	puts $fctf "Defocus deviation:      [$w.def_dev_scale get] um"
	puts $fctf "Astigmatism angle:      [$w.ast_ang_scale get] degrees"
	puts $fctf "Baseline:               [$w.base_eq.entry get]"
	puts $fctf "Envelope:               [$w.env_eq.entry get]"
	puts $fctf "First zero:             [getFirstZero] angstrom"
	puts $fctf " "
	puts $fctf "s Baseline Envelope CTF RPS"
	for {set i 0} {$i < $n} {incr i 1} {
		set s [expr $step * $i]
		set b [lindex $base $i]
		set e [lindex $envl $i]
		set c [lindex $ctf $i]
		set r [lindex $rps $i]
		puts $fctf "$s $b $e $c $r"
	}
	puts $fctf " "
	close $fctf
}


## @brief Updates CTF parameters from the micrograph parameters in memory.
#

proc getCTFparam { } {
	global project_item theimg imgtype
	global filename
    global volt Cs amp_fac base_type env_type
	set w .wctf
	if ![winfo exists $w] { return }
	if ![Bmg exists] { return }
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	if { $imgtype != "ps" } { return }
#	if [regexp "Field" $project_item] { return }
#	puts "Updating CTF parameters: $project_item"
	set volt [expr 0.001 * [Bmg get $mg_item volt]]
	set Cs [expr 1e-7 * [Bmg get $mg_item Cs]]
	set amp_fac [Bmg get $mg_item amp_fac]
	set defocus [expr 1e-4 * [Bmg get $mg_item defocus]]
#	puts "Updating CTF parameters: defocus = $defocus"
	set def_dev [expr 1e-4 * [Bmg get $mg_item defocus_deviation]]
	set ast_ang [Bmg get $mg_item astigmatism_angle]
	set base_eq [Bmg get $mg_item baseline]
	set base_type [Bmg get $mg_item baseline_type]
	set env_eq [Bmg get $mg_item envelope]
	set env_type [Bmg get $mg_item envelope_type]
#	puts "Baseline type = $base_type"
	if { $volt < 1 } { set volt 120 }
	if { $Cs < 0.000001 } { set Cs 2 }
	if { $amp_fac < 0 } { set amp_fac 0 }
	if { $amp_fac > 1 } { set amp_fac 1 }
	if { $defocus < 0.001 } { set defocus 2 }
	if { [string length $base_eq] < 1 } { set base_eq {0.2*exp(-5000*$s2) + 0.14*exp(-300*$s2) + 0.4} }
	if { [string length $env_eq] < 1 } { set env_eq {0.2*exp(-1000*$s2)} }
	set filename [Bmg get $mg_item filename ps]
#	puts "Updating CTF parameters: $mg_item"
#	puts "Power spectrum filename: -$filename-"
	wm title $w [list "CTF:" $filename "(" $n ")"]
	updateItemID $w.itemid
	$w.vca.volt_entry delete 0 end
	$w.vca.volt_entry insert 0 $volt
	$w.vca.cs_entry delete 0 end
	$w.vca.cs_entry insert 0 $Cs
	$w.vca.amp_entry delete 0 end
	$w.vca.amp_entry insert 0 $amp_fac
	$w.def_avg_scale set $defocus
	$w.def_dev_scale set $def_dev
	$w.ast_ang_scale set $ast_ang
	$w.base_eq.entry delete 0 end
	$w.base_eq.entry insert 0 $base_eq
	$w.env_eq.entry delete 0 end
	$w.env_eq.entry insert 0 $env_eq
	calcRPS
	updateCTF
#	puts "Updated CTF parameters"
}

## @brief Saves CTF parameters to the micrograph parameters in memory.
#

proc setCTFparam { } {
	set w .wctf
	if ![winfo exists $w] { return }
	global project_item theimg
    global volt Cs amp_fac base_type env_type
	global ctf_hires ctf_lores
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set volt [$w.vca.volt_entry get]
	set Cs [$w.vca.cs_entry get]
	set amp_fac [$w.vca.amp_entry get]
	set ctf_lores [$w.pp.lores_entry get]
	set ctf_hires [$w.pp.hires_entry get]
	if ![Bmg exists] { createMgParam }
	set filename [Bmg get $mg_item filename ps]
	wm title $w [list "CTF:" $filename "(" $n ")"]
	Bmg set $mg_item volt [expr 1000 * $volt]
	Bmg set $mg_item Cs [expr 1e7 * $Cs]
	Bmg set $mg_item amp_fac $amp_fac
	Bmg set $mg_item defocus [expr 1e4 * [$w.def_avg_scale get]]
	Bmg set $mg_item defocus_deviation [expr 1e4 * [$w.def_dev_scale get]]
	Bmg set $mg_item astigmatism_angle [$w.ast_ang_scale get]
	Bmg set $mg_item baseline_type $base_type
	Bmg set $mg_item baseline [$w.base_eq.entry get]
	Bmg set $mg_item envelope_type $env_type
	Bmg set $mg_item envelope [$w.env_eq.entry get]
	Bmg get $mg_item zero
#	puts "Setting CTF parameters: base_type=$base_type"
}
