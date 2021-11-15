##
# @file		bshow_voxels.tcl
#
# @brief	Procedures to record voxel values in Bshow.
#
# @author	Bernard Heymann
# @date		Created: 20101202
# @date		Modified: 20200117

## @brief Handles selections from the image data.
#

set line_width 1

proc Voxels { } {
	global helv12 theimg
	global line_width trace_line
	set c [getImageCanvas $theimg]

	set w .wvox
	
	if ![winfo exists $w] {
		toplevel $w
		
		frame $w.frame
		pack  $w.frame -expand yes -fill both -padx 1 -pady 1
		text $w.text -width 60 -height 24 -wrap none -undo 1 -background white\
			-xscrollcommand "$w.xscroll set" \
			-yscrollcommand "$w.yscroll set" \
			-setgrid 1 -highlightthickness 0 -pady 2 -padx 3
		scrollbar $w.xscroll -command "$w.text xview" \
			-highlightthickness 0 -orient horizontal
		scrollbar $w.yscroll -command "$w.text yview" \
			-highlightthickness 0 -orient vertical

		grid $w.text -in $w.frame -padx 1 -pady 1 \
			-row 0 -column 0 -rowspan 1 -columnspan 1 -sticky news
		grid $w.xscroll -in $w.frame -padx 1 -pady 1 \
			-row 1 -column 0 -rowspan 1 -columnspan 1 -sticky news
		grid $w.yscroll -in $w.frame -padx 1 -pady 1 \
			-row 0 -column 1 -rowspan 1 -columnspan 1 -sticky news
		grid rowconfig    $w.frame 0 -weight 1 -minsize 0
		grid columnconfig $w.frame 0 -weight 1 -minsize 0

		frame $w.type
		checkbutton $w.type.line -text "Line" -variable trace_line 
		setupEntry $w.type.width "Width" double $line_width "Width of line"

		frame $w.buttons
		button $w.buttons.save -text Save -command "saveAsTextFile $w.text"
#		button $w.buttons.clear -text Clear -command "$w.text delete 1.0 end"
		button $w.buttons.clear -text Clear -command "clearVoxelText $w"
		button $w.buttons.close -text Close -command "destroy $w"
		
		pack $w.type.line $w.type.width -side left -expand 1
#		pack $w.type.line -side left -expand 1
		pack $w.buttons.save $w.buttons.clear $w.buttons.close -side left -expand 1
		pack $w.buttons $w.type -side bottom -fill x -pady 2m

	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Record voxels"
    wm iconname $w "Voxels"
	.menuBar.window entryconfigure "Voxels" -state normal
	bind $c <1> "mousePressed %W %x %y"
	bind $w <Control-w> { destroy $w }
}

proc clearVoxelText { w } {
	$w.text delete 1.0 end
}

