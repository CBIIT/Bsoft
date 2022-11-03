##
# @file		bshow_param.tcl
#
# @brief	Procedures to read and write parameter files
#
# @author	Bernard Heymann
# @date		Created: 20010516
# @date		Modified: 20210311

set paramfile ""
set paramdir [pwd]
set project_root ""
set project_item ""
set project_filter "*"
set fom_cutoff 0
set tilt_axis 0
set tilt_angle 0
set tilt_offset 0
set show_tilt_axis 1
set mg_sort ""
set focal_length 3.4
set obj_apert 200.0
set slit 0.0
set direct_beam 0.0
set thickness 500.0
set emfp 2200.0
set material "Vitreous ice"

## @brief Menu for opening and closing parameter files
#

proc setupParameterMenu {} {
	global theimg
	
	.menuBar add cascade -menu .menuBar.param -label "Micrograph" -underline 0
	menu .menuBar.param -tearoff 0
	.menuBar.param add command -label "Read parameters" \
		-command { openParam } \
		-underline 0 -accelerator "Ctrl-r"
	.menuBar.param add command -label "Append parameters" \
		-command { appendParam } \
		-underline 0 -accelerator "Ctrl-a"
	.menuBar.param add command -label "Write parameters" \
		-command { writeParam } \
		-underline 0 -accelerator "Ctrl-w"
	.menuBar.param add command -label "Edit parameters" \
		-command { editParam } \
		-underline 0 -accelerator "Ctrl-e"
	.menuBar.param add command -label "Sorted micrographs" \
		-command { sortedMgTable } \
		-underline 0
	.menuBar.param add command -label "Change image" \
		-command { selectAFile } \
		-underline 0 -accelerator "Ctrl-a"
	.menuBar.param add command -label "Calculate effective MFP" \
		-command { calculateEMFP } \
		-underline 0 -accelerator "Ctrl-l"
	.menuBar.param add command -label "Fit CTF" \
		-command { CTF $theimg } \
		-underline 0 -accelerator "Ctrl-c"
	.menuBar.param add command -label "Particles" \
		-command { Particles } \
		-underline 0 -accelerator "Ctrl-p"
	.menuBar.param add command -label "Filaments" \
		-command { Filaments } \
		-underline 0 -accelerator "Ctrl-f"
	.menuBar.param add command -label "Helix" \
		-command { Helix } \
		-underline 0 -accelerator "Ctrl-h"
	.menuBar.param add command -label "Crystallography" \
		-command { Crystallography } \
		-underline 0 -accelerator "Ctrl-x"
	.menuBar.param add command -label "Tomography" \
		-command { Tomography } \
		-underline 0 -accelerator "Ctrl-t"
}


## @brief Opens a micrograph parameter file for reading.
#

proc openParam { } {
	global paramfile paramdir
	set currdir $paramdir
	if { [string length $currdir] < 2 } { set currdir [pwd] }
    set types {
		{"STAR files"		{.star .STAR}	}
		{"XML files"		{.xml .XML}	}
		{"EMX files"		{.emx .EMX}	}
	}
	set paramfile [tk_getOpenFile -title "Open parameter file" -filetypes $types \
			 -initialdir $currdir -defaultextension .star]
	if { [string length $paramfile] < 1 } { return }
	loadParam $paramfile
	Update 0
}

## @brief Opens a micrograph parameter file for reading.
#
# @param	file			parameter file name.

proc loadParam { file } {
	global paramfile paramdir
	global project_item
	if ![file exists $file] {
		tk_messageBox -message "Error: File $file not found" -type ok
		return
	}
	set project_item "all"
	set paramfile $file
	Bmg read $paramfile
	set paramdir [file dirname $paramfile]
	set paramfile [file tail $paramfile]
#	puts "$paramdir $paramfile"
	selectAFile
}

## @brief Appends one or more micrograph parameter files.
#

proc appendParam { } {
	global paramfile paramdir
	set currdir $paramdir
	if { [string length $currdir] < 2 } { set currdir [pwd] }
    set types {
		{"STAR files"		{.star .STAR}	}
		{"XML files"		{.xml .XML}	}
		{"EMX files"		{.emx .EMX}	}
	}
	set addfiles [tk_getOpenFile -title "Append parameter files" -filetypes $types \
			 -multiple true -initialdir $currdir -defaultextension .star]
	if { [string length $addfiles] < 1 } { return }
	puts $addfiles
#	loadParam $paramfile
	Bmg add $addfiles
	Update 0
}

proc tag_line { w } {
	set cur [$w index current]
	$w tag delete seline
	$w tag add seline "$cur linestart" "$cur lineend"
	$w tag configure seline -background yellow
}

## @brief Checks if the parameter file exist, otherwise save it.
#

proc checkParamFile { } {
	global paramfile paramdir
	puts "checkParamFile: $paramfile $paramdir"
	if { ![file exists $paramfile] } {
		set pf "$paramdir/$paramfile"
		if { ![file exists $pf] } {
			writeParam
		}
	}
}

## @brief Checks parameter file name and constructs new file name.
#

proc newParamFile { suffix } {
	global paramfile project_root
	if [string length $project_root] {
		set newfile $project_root
	} else {
		set newfile [file rootname $paramfile]
	}
	if ![string length $newfile] {
		tk_messageBox -message "Please write a parameter file first" -type ok
	} else {
		append newfile $suffix
	}
	return $newfile
}

proc setupFileList { w img_files } {
	global filename project_filter
	set project_filter [$w.filter.e get]
	$w.file.list delete 1.0 end
	foreach name $img_files {
		if [string match "*$project_filter*" $name] {
			$w.file.list insert end "$name"
#			puts "$filename $name"
			if [string length $filename] {
				if [regexp $filename $name] {
					$w.file.list see end
					tag_line $w.file.list
#					puts $filename
				}
			}
			$w.file.list insert end "\n"
		}
	}
	$w.file.list delete [$w.file.list index current] end
	if ![string length $filename] {
		$w.file.list mark set current 1.0
		tag_line $w.file.list
	}
#	puts [$w.file.list index current]
}

## @brief Dialog box to select an image file.
#

proc selectAFile { } {
	global filename imgtype
	global project_item project_filter
	if { [Bmg exists] } {
		set img_files [Bmg get all image_filenames]
	} elseif { [Bmodel exists] } {
		set img_files [Bmodel get all image_filenames]
	} else {
		return
	}
#	puts "$img_files"
	set oldfilename $filename
	if { [llength $img_files] < 1 } {
		puts "Error: No image files found!"
		if { [Bmodel exists] } {
			Bmodel set all map $filename
		}
	} elseif { [llength $img_files] < 2 } {
		set list [split $img_files ":"]
		set filename [lindex $list 0]
		set imgtype [lindex $list 1]
		set project_item "[lindex $list 2]:[lindex $list 3]"
#		puts "$filename ($imgtype) $project_item"
	} else {
		set w .wsel
		catch {destroy $w}
		toplevel $w
		wm title $w "Select an image"
		wm iconname $w "Select"

		label $w.msg -justify left \
			-text "Click on the desired file name and then on Select"
		pack $w.msg -side top -fill x -pady 2m
		
		frame $w.file
		scrollbar $w.file.xscroll -command "$w.file.list xview" -orient horizontal
		scrollbar $w.file.yscroll -command "$w.file.list yview"
		text $w.file.list -width 80 -height 20 -relief sunken -borderwidth 3 \
			-xscrollcommand "$w.file.xscroll set" \
			-yscrollcommand "$w.file.yscroll set" -wrap none
		pack $w.file.xscroll -side bottom -expand 0 -fill x
		pack $w.file.yscroll -side right -expand 0 -fill y
		pack $w.file.list -side left -expand 1 -fill both
		pack $w.file -side top -fill both -pady 2m -expand 1

		## Filter
		labelframe $w.filter -text "File filter" -labelanchor nw
		entry $w.filter.e
		$w.filter.e insert 0 $project_filter
		pack $w.filter.e -side left -expand yes -fill x

		frame $w.buttons
		button $w.buttons.ok -text Select -command "setSelectionAndDestroy $w"
		button $w.buttons.cancel -text Cancel -command { destroy .wsel }
		pack $w.filter -side bottom -fill x -padx 2 -pady 2
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.filter -side top -expand 1
		pack $w.buttons.ok $w.buttons.cancel -side left -expand 1

		setupFileList $w $img_files
		
		bind $w.filter.e <Return> "setupFileList $w {$img_files}"
		bind $w.file.list <1> "tag_line $w.file.list"
#		bind $w.file.list <Double-1> "setSelectionAndDestroy $w"
#		bind $w.file.list <Double-1> {
#			tag_line $w.file.list
#			setSelectionAndDestroy $w
#		}
	
		lower . .wsel
		tkwait window $w
	}
#	puts "$oldfilename --- $filename"
#	puts "$filename ($imgtype) $project_item"
	if { [string compare $oldfilename $filename] != 0 } {
		loadImage $filename -1
	}
	updateParametersAndImage
}

proc setSelectionAndDestroy { w } {
	global filename
	global imgtype
	global project_item
	set item [$w.file.list get -- seline.first seline.last]
	set list [split $item ":"]
	set filename [lindex $list 0]
	set imgtype [lindex $list 1]
	set project_item "[lindex $list 2]:[lindex $list 3]"
#	puts "$filename ($imgtype) $project_item"
	destroy $w
}

proc updateClassParam { } {
	global project_item theimg
	if [regexp "Class" $project_item] {
		set n [Bimage get $theimg nimages]
		for { set i 0 } { $i < $n } { incr i } {
			Bimage set $theimg select $i [Bmg get $project_item select [expr $i + 1]]
		}
	}
}

proc updateParametersAndImage { } {
	global filename paramfile modelfile
	global theimg imgtype
	global project_item mg_sel
	
#	puts "Updating parameters and image: $project_item"
	
	if { [Bmg exists] } {
		Bmg set all active 0
		if [regexp "Reconstruction" $project_item] { Bmg set all active 1 }
		if { $imgtype == "frame" } {
			if { [Bimage get $theimg nslices] > 1 } { switchImage }
		}
		set wc [getControlWindow $theimg]
		set n [$wc.image.scale get]
		set mg_item [micrographItem $n]
		set mg_sel [Bmg get $mg_item select]
		wm title . [list $filename ":" $paramfile]
	} elseif { [Bmodel exists] } {
		wm title . [list $filename ":" $modelfile]
	}
	
#	Bimage set $theimg label $project_item
	updatePixelsize
	updateOrigins

	if { [Bmg exists] } {
		getCTFparam
		updateClassParam
		updateBoxParam
		updateNodeParam
		updateTomoTable
		updateMarkerParam
		updateSpotParam
		updateLayerLineParam
	} elseif { [Bmodel exists] } {
		objectsDrawAll
	}
}

proc getNeighborImageOfType { dir } {
	global filename
	global theimg imgtype
	global project_item

	set wc [getControlWindow $theimg]
	set img_num [expr [$wc.image.scale get] + $dir]
	set img_max [$wc.image.scale cget -to]
	
	if { $img_num >= 0 && $img_num <= $img_max } {
		$wc.image.scale set $img_num
		updateParametersAndImage
		return
	}
	
	if { [Bmg exists] } {
		set img_files [Bmg get all image_filenames $imgtype]
	} elseif { [Bmodel exists] } {
		set img_files [Bmodel get all image_filenames]
	} else {
		return
	}
	
	set f ""
	set prevfile ""
	set prevlist ""
	set nxt 0
	foreach name $img_files {
		set list [split $name ":"]
		if [string equal $imgtype [lindex $list 1]] {
			set f [lindex $list 0]
			if { $nxt } {
				break
			}
			if [string equal $filename $f] {
				if { $dir > 0 } {
					set nxt 1
					set f ""
				} else {
					set f $prevfile
					set list $prevlist
					break
				}
			}
			set prevfile $f
			set prevlist $list
		}
	}

	if [string length $f] {
		if { [lindex $list 2] == "Micrograph" } {
			loadMicrograph [lindex $list 3]
		} else {
			set filename $f
			set project_item "[lindex $list 2]:[lindex $list 3]"
			loadImage $filename -1
			updateParametersAndImage
		}
	}
}

proc getNeighborImageOfType2 { dir } {
	global filename
	global theimg imgtype
	global project_item

	set wc [getControlWindow $theimg]
	set img_num [expr [$wc.image.scale get] + $dir]
	
	if { [Bmg exists] } {
		set img_files [Bmg get all image_filenames]
	} elseif { [Bmodel exists] } {
		set img_files [Bmodel get all image_filenames]
	} else {
		return
	}
	
	set f ""
	set curr_item ""
	set prevfile ""
	set prevlist ""
  	set nxt 0
#	puts $project_item
	foreach name $img_files {
		set list [split $name ":"]
#		puts $name
		if [string equal $imgtype [lindex $list 1]] {
			set f [lindex $list 0]
			set curr_item "[lindex $list 2]:[lindex $list 3]"
			if { $nxt } {
				break
			}
#			puts $curr_item
			if [string equal $project_item $curr_item] {
#				puts $curr_item
				if { $dir > 0 } {
					set nxt 1
				} else {
					if [string length $prevfile] { set f $prevfile }
					if [string length $prevlist] { set list $prevlist }
					break
				}
			}
			set prevfile $f
			set prevlist $list
		}
	}

	if [string length $f] {
		set project_item "[lindex $list 2]:[lindex $list 3]"
		if [string equal $filename $f] {
		} else {
			set filename $f
			loadImage $filename -1
		}
		$wc.image.scale set $img_num
		updateParametersAndImage
	}
}

## @brief Sets up micrograph or reconstruction id and total.
#
# @param	w				window in which to set the parameters.

proc setupItemID { w } {
	global project_item
	
	set idx [Bmg get $project_item index_of_item]
	set nmg [Bmg get all number_of_item]

	if ![winfo exists $w] {
		frame $w
	}

	if [regexp "Reconstruction" $project_item] {
		label $w.tag -text "Reconstruction" -anchor e
	} else {
		label $w.tag -text "Micrograph" -anchor e
	}

	label $w.str -width 50 -text $project_item -relief sunken -anchor w
	entry $w.idx -width 5 -validate key -vcmd { string is integer %P }
	$w.idx insert 0 $idx
	label $w.nmg -width 5 -text " / $nmg" -anchor w
	
	pack $w.tag $w.str $w.idx $w.nmg -side left -padx 2

	bind $w.idx <Return> { setMicrograph %W }
}

## @brief Updates micrograph or reconstruction id.
#
# @param	w				window in which to set the parameters.

proc updateItemID { w } {
	global project_item theimg
	if ![winfo exists $w] { return }
	if [regexp "Micrograph" $project_item] {
		set item $project_item
	} elseif [regexp "Reconstruction" $project_item] {
		set item $project_item
	} else {
		set wc [getControlWindow $theimg]
		set n [$wc.image.scale get]
		set item [micrographItem $n]
	}
	set idx [Bmg get $item index_of_item]
#	puts "$project_item $item $idx"
	$w.str config -text "$item"
	$w.idx delete 0 end
	$w.idx insert 0 $idx
}

## @brief Loads a micrograph based on id.
#
# @param	id 			micrograph id .

proc loadMicrograph { id } {
	global project_item theimg filename imgtype
	set oldfilename $filename
	set mg_item "Micrograph:$id"
	if [regexp "Micrograph" $project_item] {
		set project_item $mg_item
	}
	set filename [Bmg get $mg_item filename $imgtype]
#	puts "$mg_item $filename"
	if { [string compare $oldfilename $filename] != 0 } {
		loadImage $filename -1
	}
	if { $imgtype == "mg" } {
		set n [Bmg get $mg_item img_num]
		set wc [getControlWindow $theimg]
		$wc.image.scale set $n
	}
	updateParametersAndImage
}


## @brief Sets up micrograph navigation and selection buttons.
#
# @param	w				window in which to set up the buttons.
# @param	typelabel		type of navigation (micrograph).

proc setupMicrographButtons { w typelabel } {
	global mg_sel
	
	if ![winfo exists $w] {
		frame $w
	}
	
	button $w.prev -text "Previous $typelabel" \
				-relief raised -command "getNeighborImageOfType -1"
	button $w.next -text "Next $typelabel" \
				-relief raised -command "getNeighborImageOfType 1"
	checkbutton $w.select -text "Select" -variable mg_sel \
				-command { selectMicrograph }
	pack $w.prev $w.next $w.select \
			-side left -pady 2 -padx 2 -ipadx 2 -ipady 2
}

## @brief Sets up micrograph thickness fields.
#
# @param	w				window in which to set the parameters.

proc setupMicrographThickness { w } {
	global project_item theimg direct_beam emfp thickness

	if ![regexp "Micrograph" $project_item] { return }

	if ![winfo exists $w] {
		frame $w
	}
	
	label $w.ttag -text "Thickness:" -width 10 -anchor e
	label $w.thick -width 6 -text $thickness -anchor w
	label $w.dbtag -text "Direct beam/Dose:" -width 16 -anchor e
	label $w.direct -width 6 -text $direct_beam -anchor w
	button $w.dbset -text Set -command "directBeam" -relief raised
	label $w.ltag -text "EMFP:" -width 8 -anchor e
	label $w.emfp -width 6 -text $emfp -anchor w
	button $w.lamcalc -text Calculate -command "calculateEMFP" -relief raised

	pack $w.ttag $w.thick $w.dbtag $w.direct $w.dbset $w.ltag \
		$w.emfp $w.lamcalc -side left -padx 1m
	
	updateMicrographThickness $w
}

## @brief Updates micrograph thickness fields.
#
# @param	w				window in which to set the parameters.

proc updateMicrographThickness { w } {
	global project_item theimg direct_beam emfp thickness
	if ![winfo exists $w] { return }
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]

	if { $emfp < 100 } {
		set emfp [getEMFP]
	}

	set db [Bmg get $mg_item dose]
	if { $db > 0.01 } {
		set direct_beam $db
	} else {
		set db "nd"
	}
	
	set intensity [Bmg get $mg_item intensity]
	
	if { $intensity > 0.001 && $direct_beam > 0.01 } {
		set thickness [expr $emfp * log($direct_beam/$intensity)]
		set td "$thickness"
	} else {
		set td "nd"
	}
	
	$w.thick config -text $td
	$w.direct config -text $db
	$w.emfp config -text $emfp
}

## @brief Creates micrograph parameters from the image.
#

proc createMgParam { } {
	global filename theimg imgtype
	global paramfile project_item
	if [Bmg exists] { return }
	set id_name [file rootname [file tail $filename]]
#	puts "New ID: $id_name"
	set nz [Bimage get $theimg nslices]
	if { $imgtype == "" } { setImageType }
	if { $imgtype == "mg" || $imgtype == "frame" } {
		if { $nz > 1 } { switchImage }
		set nz [Bimage get $theimg nslices]
	}
	set nimg [Bimage get $theimg nimages]
	Bmg create $filename $imgtype
	if { $nz > 1 } {
		Bmg set all active 1
		if { $nimg > 1 } {
			set imgtype "part"
		} else {
			set imgtype "rec"
		}
		set project_item "Reconstruction:$id_name"
		Bmg set $project_item id $id_name
	} else {
		Bmg set all active 0
		set project_item "Micrograph:$id_name"
		Bmg set $project_item id $id_name
		if { $nimg > 1 && $imgtype != "frame" } {
			set project_item "Field:$id_name"
			Bmg set $project_item id $id_name
		}
	}
#	puts "New item: $project_item"
	set ori [Bimage get $theimg origin]
	set ps [Bimage get $theimg pixel_size 0]
	Bmg set $project_item origin [lindex $ori 0] [lindex $ori 1] [lindex $ori 2]
	Bmg set $project_item pixel_size [lindex $ps 0] [lindex $ps 1] [lindex $ps 2]
#	Bmg set $project_item defocus 2e4
#	puts "Pixel size: $ps"
#	Bimage set $theimg label $project_item
}

## @brief Returns the current project item.
#
# @param	n			micrograph index (first = 0)

proc micrographItem { n } {
	global project_item

	if { ![Bmg exists] && ![Bmodel exists] } { return "" }
	
	if [regexp "Model" $project_item] {
		return $project_item
	} elseif [regexp "Reconstruction" $project_item] {
		return $project_item
	} elseif [regexp "Class" $project_item] {
		return $project_item
	} elseif [regexp "Micrograph" $project_item] {
		set id [Bmg get $project_item id $n]
		return "Micrograph:$id"
#		return $project_item
	} elseif [regexp "Field" $project_item] {
		set id [Bmg get $project_item id $n]
		return "Micrograph:$id"
	} else {
		puts "Error: project_item not defined! ($project_item)"
	}
}

proc selectMicrograph { } {
	global project_item theimg mg_sel
	set mg_item $project_item
	if ![regexp "Micrograph" $mg_item] {
		set wc [getControlWindow $theimg]
		set n [$wc.image.scale get]
		set mg_item [micrographItem $n]
	}
#	puts "selectMicrograph: $mg_item"
	Bmg set $mg_item select $mg_sel
}

## @brief Updates an existing project for micrograph parameters in memory.
#

proc setMgParam { } {
	global filament_width filament_node_radius
	global marker_radius fom_cutoff
	global theimg imgtype
	global filename project_item
	if ![Bmg exists] { createMgParam }
#	puts "setMgParam: $project_item"
#	Bimage set $theimg label $project_item
	if { [Bimage get $theimg nslices] > 1 } {	# Converting from mg to rec
		Bmg set all active 1
		if { [string first "mg" $imgtype] > -1 } {
			set imgtype "rec"
			puts "Trying to convert to reconstruction"
		}
#		puts "Reconstruction: $imgtype $filename"
	} else {
		if { [string first "rec" $imgtype] > -1 } { set imgtype "mg" }
	}
#	puts "Bmg set $project_item filename $filename $imgtype"
	Bmg set $project_item filename $filename $imgtype
	set ps [getPixelSizeEntry]
	Bmg set $project_item pixel_size [lindex $ps 0] [lindex $ps 1] [lindex $ps 2]
#	if { [Bimage get $theimg nslices] < 2 } {	# Micrograph parameters only
#		if [winfo exists .wctf] { setCTFparam }
#	}
#	puts "setMgParam: $filename"
	setBoxSize
	if [winfo exists .wfil] {
		set filament_width [.wfil.fil_width.entry get]
		set filament_node_radius [.wfil.node_radius.entry get]
		Bmg set $project_item filament_width $filament_width
		Bmg set $project_item filament_node_radius $filament_node_radius
	}
	if [winfo exists .wtomo] {
		set marker_radius [.wtomo.marker_radius.e get]
		Bmg set $project_item marker_radius $marker_radius
	}
}

## @brief Saves the micrograph data base to a parameter file.
#

proc writeParam { } {
	global filename
	global paramfile paramdir
	set currdir $paramdir
	if { [string length $currdir] < 2 } { set currdir [pwd] }
	if { [string first $currdir ".."] > -1 } { set currdir [pwd] }
	puts "Current directory: $currdir"
	createMgParam
	set thename [file rootname [file tail $filename]]
	if { [string length $paramfile] < 1 } {
		set paramfile $thename
		if { ![string match "*.star" $paramfile] } {
			append paramfile ".star"
		}
	}
	set paramfile [file tail $paramfile]
    set types {
		{"STAR files"		{.star .STAR}	}
		{"XML files"		{.xml .XML}	}
		{"EMX files"		{.emx .EMX}	}
	}
	set paramfile [tk_getSaveFile -filetypes $types \
			-initialfile $paramfile -initialdir $currdir -defaultextension .star]
	if { [string length $paramfile] < 1 } { return }
	set paramfile [relativePath $paramdir $paramfile]
	setMgParam
	wm title . [list $filename ":" $paramfile]
#	puts "Writing $paramfile"
	Bmg write $paramfile
#	puts "Writing .bshow"
	writeSettingsFile
}

## @brief Edits a text parameter file.
#

proc editParam { } {
	global paramfile paramdir
	editTextFile $paramfile
}



## @brief Sets up the window for the sorted micrograph table.
#

proc sortedMgTable { } {
	if ![Bmg exists] { return }
	global cour12 helv12
	global mg_sort
	set w .wsort
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

#		label $w.header -font $helv12 -anchor w \
#			-text "     Mg       $mg_sort"
		frame $w.header
		label $w.header.tag -text "     Mg       "
		tk_optionMenu $w.header.menu mg_sort \
			"Intensity" "Defocus" "Water ring" "Particles"
		$w.header.menu configure -width 10

		frame $w.buttons
		button $w.buttons.update -text Update -command "updateSortedMgTable"
		button $w.buttons.close -text Close -command "destroy $w"
		
		pack $w.header -side top -fill x -padx 5
		pack $w.header.tag $w.header.menu -side left -pady 2 -anchor w
		pack $w.frame -expand yes -fill both -padx 5 -pady 5
		
		pack $w.buttons.update $w.buttons.close -side left -expand 1

		pack $w.buttons -side bottom -fill x -pady 2m
	
		$w.table tag configure selmg -background #a0b7ce
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Sorted Micrographs"
    wm iconname $w "Sorted"
	updateSortedMgTable
}

## @brief Updates the sorted micrograph table.
#

proc updateSortedMgTable {} {
	global project_item mg_sort
	if ![winfo exists .wsort] { return }
#	puts "Updating sorted micrograph table parameters"
	set m [Bmg sort $project_item $mg_sort]
	set table ""
	foreach {id val} $m {
		append table [format " %s  %7.3f\n" $id $val]
	}
	.wsort.table delete 1.0 end
	.wsort.table insert 1.0 $table
	bind .wsort.table <Double-Button-1> { setMicrographFromTable %W }
}


## @brief Loads a micrograph based on id.
#
# @param	w 			dialog window with micrograph id.

proc setMicrographFromTable { w } {
	global project_item theimg filename imgtype
	set oldfilename $filename
	tag_line $w;
	set item [$w get -- seline.first seline.last]
	loadMicrograph [lindex $item 0]
}

## @brief Loads a micrograph based on id.
#
# @param	w 			dialog window with micrograph id or number.

proc setMicrograph { w } {
	global project_item
	set idx [$w get]
	set n [Bmg get $project_item number_of_item]
	if { $idx > $n } { set idx $n }
	set id [Bmg get $project_item id_from_index $idx]
#	puts "$idx $id"
	loadMicrograph $id
}

## @brief Dialog box to set the direct beam/dose

proc directBeam { } {
	global helv12
	global project_item theimg
	set wc [getControlWindow $theimg]
	set n [$wc.image.scale get]
	set mg_item [micrographItem $n]
	set db [Bmg get $mg_item dose]
	set w .wdb
	if ![winfo exists $w] {
		toplevel $w
		
		setupEntry $w.db "Direct beam/Dose" double $db \
			"Direct beam/Applied dose in electrons/angstrom square"

		checkbutton $w.all -text "Set dose for all micrographs" \
			-variable dball -anchor w

		pack $w.db $w.all -side top -fill x -pady 2 -padx 2
		
		frame $w.buttons
		button $w.buttons.ok -text OK -command "setDirectBeam $w"
		button $w.buttons.cancel -text Cancel -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.ok $w.buttons.cancel -side left -expand 1
	
		bind $w <Return> "setTilts .wtilt"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Dose"
    wm iconname $w "Dose"
}

## @brief Sets the direct beam/dose.
#
# @param	w 			the window to destroy.

proc setDirectBeam { w } {
	global project_item theimg dball
	set db [expr [$w.db.e get]]
	destroy $w
	if { $dball } {
		Bmg set "all" dose $db
	} else {
		set wc [getControlWindow $theimg]
		set n [$wc.image.scale get]
		set mg_item [micrographItem $n]
		Bmg set $mg_item dose $db
	}
	updateMicrographThickness .wbox.thick
}

## @brief Returns the effective MFP calculated from current parameters.
#
proc getEMFP {} {
	global emfp volt focal_length obj_apert slit material
	if { $obj_apert < 1 } {
		set obj_apert 200
	}
	set emfp_material [Bmg emfp [expr 1000 * $volt] \
		[expr 1e7 * $focal_length] [expr 1e4 * $obj_apert] $slit \
		$material]
	set emfp [lindex $emfp_material 0]
	set material [lindex $emfp_material 1]
	return $emfp
}

## @brief Returns thickness if direct beam intensity is specified.
#
proc getThickness { intensity } {
	global imgtype thickness emfp direct_beam
	if { $imgtype != "mg" } { return thickness }
	if { $direct_beam < 0.1 } { return thickness }
	if { $intensity < 0.1 } { return thickness }
	set emfp [getEMFP]
	set thickness [expr $emfp * log($direct_beam/$intensity)]
	return thickness
}

proc panelEMFP { w } {
	global volt focal_length obj_apert slit emfp material

	set emfp [getEMFP]
	
	set matlist [Bmg material]

	setupEntry $w.volt "Voltage           " double $volt \
		"The microscope acceleration voltage in kV."

	setupEntry $w.focal "Focal length      " double $focal_length \
		"The microscope focal length in mm."

	setupEntry $w.objap "Objective aperture" double $obj_apert \
		"The objective aperture diameter in µm."

	setupEntry $w.slit "Slit width        " double $slit \
		"The energy filter slit width in electron volts."

	frame $w.material
	label $w.material.tag -text "Material" -width 20 -anchor w
	tk_optionMenu $w.material.set material {*}$matlist
	pack $w.material.tag $w.material.set -side left -padx 1m

	setupEntry $w.emfp "Effective MFP     " double $emfp \
		"The calculated or entered effective mean free path in angstrom."

	pack $w.volt $w.focal $w.objap $w.slit $w.material $w.emfp \
		-side top -fill x -pady 2 -padx 5
	
	bind $w.volt.e <Return> "doCalculateEMFP $w"
	bind $w.focal.e <Return> "doCalculateEMFP $w"
	bind $w.objap.e <Return> "doCalculateEMFP $w"
	bind $w.slit.e <Return> "doCalculateEMFP $w"
	bind $w.material.set <Return> "doCalculateEMFP $w"
}

## @brief Calculates the effective MFP from microscope and material parameters.
#
proc calculateEMFP { } {
	global paramfile volt focal_length obj_apert slit emfp material
 	set w .wlam
	if ![winfo exists $w] {
		toplevel $w

		panelEMFP $w
		
		frame $w.buttons
		button $w.buttons.do -text Calculate \
			-command "doCalculateEMFP $w"
		button $w.buttons.close -text Close -command "destroy $w"
		pack $w.buttons -side bottom -fill x -pady 2m
		pack $w.buttons.do $w.buttons.close -side left -expand 1
	} else {
		wm deiconify $w
		raise $w
	}
	wm title $w "Effective MFP"
	wm iconname $w "EMFP"
}

## @brief Calculates lamda from microscope parameters.
#
# @param	w 			dialog window for parameters.

proc doCalculateEMFP { w } {
	global paramfile volt focal_length obj_apert slit emfp material
	set volt [$w.volt.e get]
	set focal_length [$w.focal.e get]
	set obj_apert [$w.objap.e get]
	set slit [$w.slit.e get]
#	set material [$w.material.e get]
	set emfp [getEMFP]
#	puts "EMFP = $emfp"
	updateEMFP $w
#	destroy $w
}

## @brief Updates the lamda value from the global variable.
#
# @param	w 			dialog window for update.

proc updateEMFP { w } {
	global emfp
	if ![winfo exists $w] { return }
	$w.emfp.e delete 0 end
	$w.emfp.e insert 0 $emfp
}

## @brief Updates the thickness from the global variable.
#
# @param	w 			dialog window for update.

proc updateThickness { w } {
	global thickness
	if ![winfo exists $w] { return }
	$w.thick.e delete 0 end
	$w.thick.e insert 0 $thickness
}

