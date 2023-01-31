##
# @file		brun_panel.tcl
#
# @brief	Panel to construct and run Bsoft command lines
#
# @author	Bernard Heymann
# @date		Created: 19990822
# @date		Modified: 20191023

set font {Helvetica 14}
set bigfont {Helvetica 24}
set smallfont {Helvetica 10}

set bshow_script $Bsoft/tcltk
set bshow_lib $Bsoft/lib
set bsoft_bin $Bsoft/bin

## @brief Writes a command to a file.
#
# @param	filename	file name.
# @param	mode		access mode.
proc saveCommand { filename mode } {
	global command

	if { [catch {open $filename $mode} f] } { return }
	
	puts $f "$command\n"
	
	close $f
}

## @brief Selects a file for reading or writing using a standard dialog
#
# @param	operation 	"save" or "append".

proc fileDialog { operation } {
	set filename [tk_getSaveFile]

	if { [string length $filename] < 1 } { return }

	if { [string compare $operation "save"] == 0 } {
		saveCommand $filename "w"
	} elseif { [string compare $operation "append"] == 0 } {
		saveCommand $filename "a"
	}
}

## @brief Constructs the command line from a list containing the program name, the options, and the file names
#
# @param c 			the entry for the command
# @param f 			the frame for files
# @param t 			the text with options

proc constructCommand { c f t } {
	set cmd [lindex [$c get] 0]

	foreach w [winfo children $t] {
		if [winfo exists $w.label] {
			upvar [$w.label cget -variable] on
			set opt [$w.label cget -text]
			if { $on == 1 } {
				append cmd " " $opt " " [$w.e1 get]
			}
		}
	}

	foreach w [winfo children $f] {
		upvar [$w.label cget -variable] on
		if { $on == 1 } {
			append cmd " " [$w.name get]
		}
	}

	$c delete 0 end
	$c insert 0 $cmd
}

## @brief Parses the usage line
#
# @param c 			the entry for the command
# @param f 			the frame for files
# @param t 			the text with options
# @param line 		the line

proc parseUsageLine { c f t line } {
	global font

	destroy .usage
	label .usage -font $font -justify left -text $line -background yellow

	set words [split $line]
	for { set i 3 } { $i < [llength $words] } { incr i } {
		set w $f.$i
		frame $w
		set fn [lindex $words $i]
		checkbutton $w.label -text $fn -variable $fn -relief flat \
			-background cyan -command "constructCommand $c $f $t"
		$w.label deselect
		entry $w.name -width 50
    	button $w.button -text "Browse ..." -command "getFilename $w.name"
		pack $w.label $w.name $w.button -side left -pady 1 -padx 4
		pack $w -side top -fill x -expand no
		bind $w.name <KeyRelease> "constructCommand $c $f $t"
	}

	if { $i > 3 } { pack $f -side bottom -pady 2 -fill x -expand no }
	pack .usage -side bottom -pady 2 -anchor w
}

## @brief Parses one line starting with an option
#
# @param c 			the entry for the command
# @param f 			the frame for files
# @param t 			the text with options
# @param line 		the line

proc parseOptionLine { c f t line } {
	global font smallfont
	if [regexp "^--" $line] { return }

	set words [split $line]
	set b [lindex $words 0]
	set v [lindex $words 1]

	set d ""
	foreach word [lrange $words 2 end] {
		if { $word != "" } {
			lappend d $word
		}
	}

	set w $t.$b
	frame $w
	entry $w.e1 -width [string length $v]
	$w.e1 insert 0 $v
	
	checkbutton $w.label -text "$b" -variable $b -relief flat \
			-background cyan -command "constructCommand $c $f $t"
	$w.label deselect

	if { $v > " " } {
		pack $w.label $w.e1 -side left -pady 1 -padx 4
	} else {
		pack $w.label -side left -pady 1 -padx 4
	}

	$t window create end -window $w -padx 4 -pady 1
	$t insert end [join $d " "]
	$t insert end "\n"
	
	bind $w.e1 <KeyRelease> "constructCommand $c $f $t"
}

## @brief Parses a command line to fill in options
#
# @param t 			text window name
# @param cmd 		command line

proc parseCommand { t cmd } {
	set vf 0
	foreach word $cmd {
		if ![string first "-" $word] {
			set vf 0
#			set word [string trimleft $word "-"]
			foreach w [winfo children $t] {
				if [regexp -- $word $w] {
					$w.label select
					if { [$w.e1 get] > " " } {
						set vf 1
					}
					break
				}
			}
		} elseif { $vf } {
			$w.e1 delete 0 end
			$w.e1 insert 0 $word
			set vf 0
		}
	}
}

## @brief Selects the program and get the options from the program
#
# @param w 			window name
# @param cmd 		command line

proc selectProgram { w cmd } {
	global font
	global smallfont

	set wn $w
	if [string equal "." $w] { set wn "" }

	set c $wn.cmd
	set f $wn.file
	set t $wn.txt
	set g $wn.grid

	destroy $f
	destroy $t
	destroy $g
	destroy $wn.hscroll
	destroy $wn.vscroll

	set prog [lindex $cmd 0]
#	puts "Program: $prog"
	if { [ string first "brun" $cmd ] > -1 } { return }
	if ![ file executable $prog ] { return }

	catch { exec $prog } result

	$c.e delete 0 end
	$c.e insert 0 $cmd

	frame $f

	frame $g
	pack $g -side bottom -expand yes -fill both -padx 1 -pady 1

	grid rowconfig    $g 0 -weight 1 -minsize 0
	grid columnconfig $g 0 -weight 1 -minsize 0

	text $t -wrap none -relief sunken -bd 2 \
			-xscrollcommand "$wn.hscroll set" \
			-yscrollcommand "$wn.vscroll set" \
			-setgrid 1 -highlightthickness 0 -pady 2 -padx 3

	scrollbar $wn.hscroll -orient horizontal -command "$t xview"
	scrollbar $wn.vscroll -orient vertical -command "$t yview"

	grid $t -padx 1 -in $g -pady 1 \
		-row 0 -column 0 -rowspan 1 -columnspan 1 -sticky news
	grid $wn.vscroll -in $g -padx 1 -pady 1 \
		-row 0 -column 1 -rowspan 1 -columnspan 1 -sticky news
	grid $wn.hscroll -in $g -padx 1 -pady 1 \
		-row 1 -column 0 -rowspan 1 -columnspan 1 -sticky news

	$t insert end "Description:\n"
	foreach line [split $result "\n"] {
		if ![string first "-" $line] {
			parseOptionLine $c.e $f $t $line
		} elseif ![string first "Usage" $line] {
			parseUsageLine $c.e $f $t $line
		} else {
			$t insert end "$line\n"
		}
	}

	parseCommand $t $cmd
}

## @brief Generates the program interface
#
# @param w 			window name
# @param cl			command line

proc programInterface { w cl } {
	global bsoft_bin

	wm title $w "Bsoft program interface"

	set wn $w	
	if [string equal "." $w] { set wn "" }

	$w configure -menu $wn.menu

	menu $wn.menu -tearoff 0
	menu $wn.menu.file -tearoff 0
	$wn.menu add cascade -label "File" -menu $wn.menu.file -underline 0
	$wn.menu.file add command -label "Save ..." -command { fileDialog "save" } \
		-underline 0 -accelerator "Ctrl-s"
	$wn.menu.file add command -label "Append ..." -command { fileDialog "append" } \
		-underline 0 -accelerator "Ctrl-a"
	$wn.menu.file add command -label "Quit" -command { exit } -underline 0 \
		-accelerator "Ctrl-q"

	bind all <Control-q> { exit }
	bind all <Command-q> { exit }

#	puts "cl: $cl"

	set command "Select a program and fill in the options"
	set program ""
	if [string length $cl] {
		set command "$bsoft_bin/$cl"
		set p "$bsoft_bin/"
		append p [lindex $cl 0]
		if [ file executable $p ] { set program $p }
	}
#	puts "Command: $command"
#	puts "Program: $program"
	
	set c $wn.cmd
	frame $c
	entry $c.e -width 100 -xscrollcommand "$c.s set" -textvariable command
	scrollbar $c.s -relief sunken -orient horiz -command "$c.e xview"
	pack $c.s $c.e -side bottom -fill x -expand 1 -padx 1m -pady 1

	frame $wn.buttons
	menubutton $wn.buttons.program -text "Program" -underline 0 -direction below \
			-menu $wn.buttons.program.m -relief raised -width 20 \
			-textvariable program
	menu $wn.buttons.program.m -tearoff 0
	set file_list [lsort [glob [file join $bsoft_bin/*]]]
	set i 0
	foreach p [split $file_list] {
		set prog [file tail $p] 
		if { [file executable $p] && [file isfile $p] } {
			$wn.buttons.program.m add command -label $prog \
				-command "selectProgram $w $p" -columnbreak [expr ($i % 20) == 0]
			incr i
		}
	}

	if [string length $command] { selectProgram $w $command }

#	set logfile "b.log"

	button $wn.buttons.run -text Execute -command { executeCommand $command "b.log" }
	button $wn.buttons.dismiss -text Dismiss -command "destroy $w"
	
	pack $wn.buttons.program $wn.buttons.run $wn.buttons.dismiss -side left -expand yes
	pack $wn.buttons $c -side top -fill x -pady 2m -expand no

	bind . <Return> { selectProgram . [.cmd.e get] }
}


