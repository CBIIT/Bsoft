##
# @file		btext.tcl
#
# @brief	Window to show text information
#
# @author	Bernard Heymann
# @date		Created: 20010516
# @date		Modified: 20130123

set textfile ""

## @brief Open a text editing window
#
# @param	w				window name to contain ".text" window.

proc openTextWindow { w } {
	global textfile

	set wn $w
	if [string equal $w "."] { set wn "" }
	set twin $wn.text
#	puts "Window name: $wn"
	
	if { ![winfo exists $w] || $w == "." } {
		if { $w != "." } { toplevel $w }

		set m $wn.menu
#		puts "Text window menu: $m"
		menu $m -tearoff 0
		$m add cascade -menu $m.file -label "File" -underline 0
		menu $m.file -tearoff 0
		$w configure -menu $m
		$m.file add command -label "Open" \
				-command "openTextFile $twin" -underline 0
		$m.file add command -label "Save" \
				-command [ list "saveTextFile" $twin $textfile ] -underline 0
		$m.file add command -label "Save as ..." \
				-command "saveAsTextFile $twin" -underline 0
		$m.file add command -label "Close" \
				-command "$wn.text delete 1.0 end" -underline 0

		frame $wn.frame
		pack  $wn.frame -expand yes -fill both -padx 1 -pady 1
		text $wn.text -width 80 -height 60 -wrap none -undo 1 -background white\
			-xscrollcommand "$wn.xscroll set" \
			-yscrollcommand "$wn.yscroll set" \
			-setgrid 1 -highlightthickness 0 -pady 2 -padx 3
		scrollbar $wn.xscroll -command "$wn.text xview" \
			-highlightthickness 0 -orient horizontal
		scrollbar $wn.yscroll -command "$wn.text yview" \
			-highlightthickness 0 -orient vertical

		grid $wn.text -in $wn.frame -padx 1 -pady 1 \
			-row 0 -column 0 -rowspan 1 -columnspan 1 -sticky news
		grid $wn.xscroll -in $wn.frame -padx 1 -pady 1 \
			-row 1 -column 0 -rowspan 1 -columnspan 1 -sticky news
		grid $wn.yscroll -in $wn.frame -padx 1 -pady 1 \
			-row 0 -column 1 -rowspan 1 -columnspan 1 -sticky news
		grid rowconfig    $wn.frame 0 -weight 1 -minsize 0
		grid columnconfig $wn.frame 0 -weight 1 -minsize 0
		
	} else {
		wm deiconify $wn
		raise $wn
    }
    wm title $w "Text Editor"
    wm iconname $w "Ted"
	bind $w <Control-z> { $wn.text edit undo }
	bind $w <Control-Z> { $wn.text edit redo }
	bind $w <Command-Z> { $wn.text edit redo }
	bind $w <Control-w> { destroy $w }
	bind $w <Command-w> { destroy $w }
	return $wn.text
}

## @brief Opens a text window with the supplied information
#
# @param	theinfo		text information.

proc showInfo { theinfo } {
	set w .info
	set wtxt [openTextWindow $w]
    wm title $w "Bsoft information window"
    wm iconname $w "Bsoft info"
    $wtxt delete 1.0 end
    $wtxt insert 1.0 $theinfo
    $wtxt mark set insert 1.0
	return $wtxt
}

## @brief Reads a text file for editing
#
# @param	filename		file name (including path).

proc editTextFile { filename } {
	if { [string length $filename] < 1 } { return }
	set w .wtxt
	readTextFile [openTextWindow $w] $filename
	set thename [file tail $filename]
    wm title $w $thename
    wm iconname $w $thename
}

## @brief Opens a text file for editing
#
# @param	wtxt				text window.

proc openTextFile { wtxt } {
	set filename [tk_getOpenFile -initialdir [pwd]]
#	set filename [file tail $filename]
#	puts "Open in $wtxt: $filename"
	readTextFile $wtxt $filename
}

## @brief Reads a text file
#
# @param	wtxt 			a text window
# @param	filename		file name (including path).

proc readTextFile { wtxt filename } {
	global textfile
	if { [string length $filename] < 1 } { return }
	if { [catch {open $filename r} ftxt] } {
		tk_messageBox -message "File $filename not opened!" -icon error -type ok
		return
	}
	set thetext [read $ftxt]
	close $ftxt
    $wtxt delete 1.0 end
    $wtxt insert 1.0 $thetext
	set textfile $filename
}

## @brief Saves a text file using its current name
#
# @param	wtxt 			a text window
# @param	filename		file name (including path).

proc saveTextFile { wtxt filename } {
	global textfile
	if { [string length $filename] < 1 } { return }
	if { [catch {open $filename w} ftxt] } { return }
	set thetext [$wtxt get 1.0 end]
	puts $ftxt $thetext
	close $ftxt
	set textfile $filename
	return $filename
}

## @brief Saves a text file with a new name
#
# @param	wtxt 			a text window

proc saveAsTextFile { wtxt } {
	global textfile
	set filename [tk_getSaveFile -initialdir [pwd] -initialfile [file tail $textfile]]
	saveTextFile $wtxt $filename
	return $filename
}

## @brief Updates text originating from an open channel
#
# @param	wtxt 			a text window
# @param	chan 			the channel

proc updateText { wtxt chan } {
	$wtxt configure -state normal
	$wtxt insert end [read $chan]
	$wtxt configure -state disabled
	$wtxt see end
    if {[eof $chan]} {
		puts "Closing text output channel"
		close $chan
		$wtxt configure -state normal
    }
}

## @brief Executes a command and sends the output to a channel for updating
#
# @param	cmd 			the command (with arguments)
# @param	logfile 		a text file to write to.

proc executeCommand { cmd logfile } {
	global textfile
	set p [lindex $cmd 0]
	if [ file executable $p ] {
		set textfile $logfile
		puts "Executing: $cmd"
		set wtxt [showInfo "$cmd\n\n"]
		set chan [open |[concat $cmd 2>@1]]
		fconfigure $chan -blocking 0 -buffering line
		fileevent $chan readable [list updateText $wtxt $chan]
	} else {
		puts "$p not executable!"
	}
}

## @brief Gets a file name entry
#
# @param	w 			an entry window

proc getFilename { w } {
	set name [tk_getOpenFile]
	if { $name > " " } {
		$w delete 0 end
		$w insert 0 $name
	}
}

## @brief Gets a directory entry
#
# @param	w 			an entry window

proc getDirectoryName { w } {
	set dir [tk_chooseDirectory \
        -initialdir ~ -title "Choose a directory"]
	if { $dir > " " } {
		$w.e delete 0 end
		$w.e insert 0 $dir
	}
}

