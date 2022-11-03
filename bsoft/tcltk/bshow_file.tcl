##
# @file		bshow_file.tcl
#
# @brief	Procedures to manage file I/O
#
# @author	Bernard Heymann
# @date		Created: 20010516
# @date		Modified: 20220224

set img_load_phase 0
set img_filter "*"
#set browser_dir "./"
set browser_dir ""

## @brief Menu for opening and closing image files and quiting
#

proc setupFileMenu {} {
	global filename
	
	.menuBar add cascade -menu .menuBar.file -label "File" -underline 0
	menu .menuBar.file -tearoff 0
	.menuBar.file add command -label "Open ..." -command { fileDialog "open" } \
		-underline 0 -accelerator "Ctrl-o"
	.menuBar.file add command -label "Save" -command { fileDialog "save" } \
		-underline 0
	.menuBar.file add command -label "Save as ..." -command { fileDialog "saveAs" } \
		-underline 0 -accelerator "Ctrl-S"
	.menuBar.file add command -label "Save region ..." -command { extractImage } \
		-underline 0
	.menuBar.file add command -label "File browser" -command { fileBrowserWindow .fb } \
		-underline 0
	.menuBar.file add command -label "Revert" -command { loadImage $filename -1 } \
		-underline 0 -accelerator "Ctrl-u"
	.menuBar.file add command -label "Get info" -command { getInfo $filename } \
		-underline 0 -accelerator "Ctrl-i"
	.menuBar.file add command -label "Write PostScript" -command { writePostScript } \
		-underline 0
	.menuBar.file add command -label "Quit" -command { Quit } -underline 0 \
		-accelerator "Ctrl-q"
}

#   Type names		Extension(s)	Mac File Type(s)
#
#---------------------------------------------------------
set filetypes {
	{"All files"		*}
	{"CCP4"		{.map .MAP}		}
	{"MRC"		{.mrc .MRC .mrcs .stk .STK .st .ST .ali .ALI .rec .REC}		}
	{"PIF"		{.pif .PIF}		}
	{"GRD"		{.grd .GRD}		}
	{"EER"		{.eer .EER}		}
	{"EM"		{.em  .EM }		}
	{"IMAGIC"	{.img .IMG .hed .HED}	}
	{"SPIDER"	{.spi .SPI}		}
	{"SUPRIM"	{.sup .SUP .spm .SPM}	}
	{"DM3"		{.dm3 .DM3}		}
	{"DM4"		{.dm4 .DM4}		}
	{"SER"		{.ser .SER}		}
	{"TIFF"		{.tif .tiff .TIF .TIFF}	}
	{"PNG"		{.png .PNG}		}
	{"PNM"		{.pbm .pgm .ppm .PBM .PGM .PPM}		}
	{"JPEG"		{.jpg .jpeg .JPG .JPEG}	}
	{"ASCII"	{.txt .asc .ascii .TXT .ASC .ASCII}	}
}

## @brief Selects a file for reading or writing using a standard dialog
#
# @param	operation 	"open", "save" or "saveAs".

proc fileDialog { operation } {
	global theimg
	global filename imgdir
	global filetypes
	global paramfile
	set currdir $imgdir
	if { [string length $currdir] < 2 } { set currdir [pwd] }
#	puts "Image directory: $imgdir"
	if { [string compare $operation "open"] == 0 } {
		set new_filename [tk_getOpenFile -title "Open image file" \
			-filetypes $filetypes -initialdir $currdir]
		if { [string length $new_filename] < 1 } { return }
		set filename [relativePath [pwd] $new_filename]
		set imgdir [file dirname $filename]
#		tk_messageBox -icon info -type ok -title "Directory and file name" -message \
#			"imgdir = $imgdir\nfilename = $filename"
		loadImage $filename -1
    } else {
		if { [string compare $operation "saveAs"] == 0 } {
			set new_filename [tk_getSaveFile -filetypes $filetypes \
				-initialfile [file tail $filename] -initialdir $currdir \
				-defaultextension .mrc]
			if { [string length $new_filename] < 1 } { return }
        	set filename [relativePath [pwd] $new_filename]
		}
		Bimage save $theimg $filename
		writeSettingsFile
    }
#	puts "Image directory: $imgdir"
#	tk_messageBox -icon info -type ok -title "Directory and file name" -message \
#		"$imgdir $filename"
	if [Bmg exists] { setMgParam }
	wm title . [list $filename ":" $paramfile]
}


## @brief	Populates a file browser tree node
#
# @param	tree		File browser tree.
# @param	node		File browser node.

proc openNode {tree node} {
	global img_filter browser_dir

    if {[$tree set $node type] ne "directory"} {
		return
    }

    set path [$tree set $node fullpath]
    $tree delete [$tree children $node]

    foreach f [lsort -dictionary [glob -nocomplain -dir $path *]] {
#    	set r [string match "./$img_filter" $f]
#   	puts "$f $r"
		if [string match "$browser_dir/$img_filter" $f] {
			set type [file type $f]
			set id [$tree insert $node end -text [file tail $f] \
				-values [list $f $type]]

			if {$type eq "directory"} {
	    	## Make it so that this node is openable
	    		$tree insert $id 0 -text dummy ;# a dummy
	    		$tree item $id -text [file tail $f]/
	
			} elseif {$type eq "file"} {
	    		set size [file size $f]
	    		## Format the file size nicely
	    		if {$size >= 1024*1024*1024} {
					set size [format %.1f\ GB [expr {$size/1024/1024/1024.}]]
	    		} elseif {$size >= 1024*1024} {
					set size [format %.1f\ MB [expr {$size/1024/1024.}]]
	    		} elseif {$size >= 1024} {
					set size [format %.1f\ kB [expr {$size/1024.}]]
	    		} else {
					append size " bytes"
	    		}
				$tree set $id size $size
			}
		}
    }

    # Stop this code from rerunning on the current node
    $tree set $node type processedDirectory
}

## @brief	Points the file browser to the current directory
#
# @param	w			File browser tree.

proc setupFileBrowser { w } {
	global img_filter browser_dir
	
	$w.tree set {} fullpath $browser_dir
	$w.tree set {} type "directory"
	
	set img_filter [$w.filter.e get]
	if { $img_filter eq "" } {
		set img_filter "*"
		$w.filter.e insert 0 $img_filter
	}

	set browser_dir [$w.dir.e get]
	if { $browser_dir eq "" } {
		set browser_dir "./"
		$w.dir.e insert 0 $browser_dir
	}

	openNode $w.tree {}
}

## @brief	Opens a file indicated by a double-click in the file browser
#
# @param	wtree		File browser tree.

proc openFileFromTree { wtree } {
	global filename imgdir
	set type [$wtree set [$wtree focus] type]
	if { $type eq "file" || $type eq "link" } { 
		set fn [$wtree set [$wtree focus] fullpath]
		set tbt [exec file $fn]
		if [string match "*text*" $tbt] {
			editTextFile $fn
		} else {
			set ft [fileType $fn]
			if { $ft == 1 } {
				set filename [relativePath [pwd] $fn]
				set imgdir [file dirname $filename]
				loadImage $fn -1
			} else {
				tk_messageBox -icon info -type ok \
					-title "Error" -message "File type not supported!"
			}
		}
	}
}

proc changeBrowserDirectory { w } {
	global browser_dir
	
#	file normalize $browser_dir

	set browser_dir [tk_chooseDirectory \
        -initialdir $browser_dir -title "Choose a directory"]
 
 	$w.dir.e delete 0 end
	$w.dir.e insert 0 $browser_dir

	setupFileBrowser $w
}

## @brief	Sets up the file browser window
#
# @param	w		Window name.

proc fileBrowserWindow { w } {
	global img_filter browser_dir
	
	if { $browser_dir eq "" } {
		set browser_dir [pwd]
	}

	catch {destroy $w}
	toplevel $w
	wm title $w "File browser"
	wm iconname $w "Browser"
	wm geometry $w 500x800+500+50

	## Buttons
	ttk::frame $w.buttons
	ttk::button $w.buttons.dir -text "Change directory" \
		-command [list changeBrowserDirectory $w]
	ttk::button $w.buttons.close -text "Close" \
		-command [list destroy [winfo toplevel $w]]
	ttk::button $w.buttons.code -text "Refresh" \
		-command [list setupFileBrowser $w]
	pack $w.buttons.dir $w.buttons.code $w.buttons.close -side left
	pack $w.buttons -side bottom -padx 2 -pady 2

	## Filter
	labelframe $w.filter -text "File filter" -labelanchor nw
	ttk::entry $w.filter.e
	$w.filter.e insert 0 $img_filter
	pack $w.filter.e -side left -expand yes -fill x
	pack $w.filter -side bottom -fill x -padx 2 -pady 2

	## Directory
	labelframe $w.dir -text "Directory" -labelanchor nw
	ttk::entry $w.dir.e
	$w.dir.e insert 0 $browser_dir
	pack $w.dir.e -side left -expand yes -fill x
	pack $w.dir -side bottom -fill x -padx 2 -pady 2

	## Create the tree and set it up
	ttk::treeview $w.tree -columns {fullpath type size} -displaycolumns {size} \
		-yscroll "$w.vsb set" -xscroll "$w.hsb set"
	if {[tk windowingsystem] ne "aqua"} {
		ttk::scrollbar $w.vsb -orient vertical -command "$w.tree yview"
		ttk::scrollbar $w.hsb -orient horizontal -command "$w.tree xview"
	} else {
	    scrollbar $w.vsb -orient vertical -command "$w.tree yview"
		scrollbar $w.hsb -orient horizontal -command "$w.tree xview"
	}
	$w.tree heading \#0 -text "File"
	$w.tree heading size -text "Size"
	$w.tree column size -stretch 0 -width 70
	
	setupFileBrowser $w

	## Arrange the tree and its scrollbars in the toplevel
	lower [ttk::frame $w.gframe]
	pack $w.gframe -fill both -expand 1
	grid $w.tree $w.vsb -sticky nsew -in $w.gframe
	grid $w.hsb -sticky nsew -in $w.gframe
	grid columnconfigure $w.gframe 0 -weight 1
	grid rowconfigure $w.gframe 0 -weight 1

	bind $w.tree <<TreeviewOpen>> {openNode %W [%W focus]}
	bind $w.tree <Double-Button-1> {openFileFromTree %W }
	bind $w <Return> [list setupFileBrowser $w]
}

## @brief	Loads the selected file into a photo image for display.
#
# @param	filename	Image file name.
# @param	img_num		Sub-image number.

proc loadImage { filename img_num } {
	global theimg
	global imgdir
	global paramfile modelfile
	global mode
	global window_width window_height window_scale rem_scale
	global img_load_phase
	global avg_check kernel
	if { [string length $filename] < 1 } { return }
	if { ![winfo exists .] } { return }
	set imgdir [file dirname $filename]
	set filename [file tail $filename]
#	tk_messageBox -icon info -type ok -title "Directory and file name" -message \
#		"$imgdir $filename"
	if { [Bmg exists] } {
		wm title . [list $filename ":" $paramfile]
	} elseif { [Bmodel exists] } {
		wm title . [list $filename ":" $modelfile]
	} else {
		wm title . $filename
	}
#	puts "Opening image $imgdir/$filename ($img_num)"

	if { [catch { set nimages [Bimage open $theimg $imgdir/$filename $img_num] } err] } {
#		bgerror $err
		fileDialog "open"
		return
	}
#	puts "Image opened with $nimages sub-images"
	if { $nimages < 1 } { return } 
	set nslices [Bimage get $theimg nslices]
#	puts "Image opened with $nslices slices"
	if { $nslices < 1 } { return }

	if { $avg_check } {
#		puts "Average filtering kernel: $kernel"
		Bimage average $theimg $kernel $kernel $kernel
	}
	
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
#	puts "Loading: $window_width $window_height $window_scale"
	set img_load_phase 1
	set width [Bimage get $theimg width]
	set height [Bimage get $theimg height]
	set origin [Bimage get $theimg origin]
#	puts $origin
	if { $img_num < 0 } { set img_num 0 }
#	puts "Image opened with origin $origin"
	set oz [lindex $origin 2]
	if { $oz < 0 } { set oz 0 }
	if { $oz >= $nslices } { set oz 0 }
	set scale $window_scale
#	$wc.scale.scale set $window_scale
	if { $window_width < 2 || $rem_scale == 0 } {
		set scale 1.0
		set window_width 512
		set window_height 512
		if { $width > $window_width } { set scale [expr 0.1*(($window_width*10)/$width) ] }
	}
	if { $scale < 0.1 } { set scale 0.1 }
	set window_scale $scale
#	puts "Image scale = $scale  mode = $mode"
	set width [expr $width * $scale ]
	set height [expr $height * $scale ]
	set size "0 0 $width $height"
	if { $width < 200 } { set width 200 }
	if { $height < 200 } { set height 200 }
	if { $width > $window_width } { set width $window_width }
	if { $height > $window_height } { set height $window_height }
	set window_width $width
	set window_height $height
#	puts "Width $width  height $height  scrollregion $size"
	$c configure -width $width -height $height -scrollregion $size
#	$c configure -scrollregion $size
#	puts "Loading2: $window_width $window_height $window_scale"
	reconfigureRangeStep
#	puts "Range step reconfigured"
	reconfigureScales
#	puts "Scales reconfigured"
	if { ![ regexp {^([0-9.]+)$} $oz ] } { set oz 0 }
	$wc.image.scale set $img_num
	$wc.slice.scale set $oz
	$wc.scale.scale set $scale
	updatePixelsize
#	puts "Pixel size updated"
	update
	set img_load_phase 0
	Update 1
}

## @brief Derives the relative path from an absolute one relative to a given directory.
#
# @param	cdir			Current directory.
# @param	path			Absolute path to convert.

proc relativePath { cdir path } {
	if { $cdir == "/" } { return $path }
	set dlist [file split $cdir]
	set plist [file split $path]
	set len [llength $dlist]
#	tk_messageBox -icon info -type ok -title "Current directory and path in relativePath" -message \
#		"Directory = $cdir ($len) Path = $path [lindex $dlist 0]"
#	puts "Current directory: $cdir  Path: $path"
	if { [lindex $dlist 0] == [lindex $plist 0] } {
		for { set i 0 } { ( $i < $len ) && ( [lindex $dlist $i] == [lindex $plist $i] ) } { incr i 1 } {
		}
		set path ""
		for { set j $i } { $j < $len } { incr j 1 } {
			append path "../"
		}
		set len [llength $plist]
		for { set j $i } { $j < $len } { incr j 1 } {
			set path [file join $path [lindex $plist $j]]
		}
	} elseif { [lindex $plist 0] == ".." } {
	}
#	puts "Relative path: $path"
#	tk_messageBox -icon info -type ok -title "Modified path in relativePath" -message "$path"
	return $path
}

## @brief Loads the information for the current image file into a text window for display.
#
# @param	filename			Image file name.

proc getInfo { filename } {
	global Bsoft
	global imgdir
#	puts "$Bsoft $imgdir $filename"
	if ![winfo exists .info] {
		toplevel .info
		frame .info.buttons
		pack .info.buttons -side bottom -fill x
		button .info.buttons.dismiss -text Dismiss \
            -default active -command "destroy .info"
		pack .info.buttons.dismiss -side left -expand 1 -pady 2
		frame .info.frame
		pack  .info.frame -expand yes -fill both -padx 1 -pady 1
		text .info.text -height 40 -width 100 -wrap word -background white\
			-xscrollcommand ".info.xscroll set" \
			-yscrollcommand ".info.yscroll set" \
			-setgrid 1 -highlightthickness 0 -pady 2 -padx 3
		scrollbar .info.xscroll -command ".info.text xview" \
			-highlightthickness 0 -orient horizontal
		scrollbar .info.yscroll -command ".info.text yview" \
			-highlightthickness 0 -orient vertical

		grid .info.text -in .info.frame -padx 1 -pady 1 \
			-row 0 -column 0 -rowspan 1 -columnspan 1 -sticky news
		grid .info.yscroll -in .info.frame -padx 1 -pady 1 \
			-row 0 -column 1 -rowspan 1 -columnspan 1 -sticky news
		grid rowconfig    .info.frame 0 -weight 1 -minsize 0
		grid columnconfig .info.frame 0 -weight 1 -minsize 0
	} else {
		wm deiconify .info
		raise .info
    }
	set filename [file tail $filename]
    wm title .info $filename
    wm iconname .info $filename
#    puts "$Bsoft/bin/bhead -verbose 7 -info $imgdir/$filename"
	catch { exec $Bsoft/bin/bhead -verbose 7 -info $imgdir/$filename } theinfo
    .info.text delete 1.0 end
    .info.text insert 1.0 $theinfo
    .info.text mark set insert 1.0
}

## @brief Writes the current canvas to a postscript file.
#

proc writePostScript { } {
	global theimg
	global imgdir
	global filename
	set currdir $imgdir
	if { [string length $currdir] < 2 } { set currdir [pwd] }
#	puts "Image directory: $imgdir"
	set init_filename [file rootname $filename]
	append init_filename .ps
    set ps_filename [tk_getSaveFile -filetypes {{"PostScript"	{.ps	.PS}		}} \
		-initialfile $init_filename -initialdir $currdir -defaultextension .ps]
	if { [string length $ps_filename] < 1 } { return }
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set scale [$wc.scale.scale get]
	set w [expr [Bimage get $theimg width] * $scale]
	set h [expr [Bimage get $theimg height] * $scale]
	$c postscript -file $ps_filename -pageanchor sw -x 0 -y 0 -width $w -height $h
#	tk_messageBox -icon info -type ok -title "Directory and file name" -message \
#		"$imgdir $filename"
}


