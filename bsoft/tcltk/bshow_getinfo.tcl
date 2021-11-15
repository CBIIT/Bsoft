##
# @file		bshow_getinfo.tcl
#
# @brief	Procedures to show image information
#
# @author	Bernard Heymann
# @date		Created: 20010516
# @date		Modified: 030523

## @brief Loads the information for the current image file into a text window for display.
#
# @param	filename			File name.

proc getInfo {filename} {
	if ![winfo exists .info] {
		toplevel .info
		frame .info.buttons
		pack .info.buttons -side bottom -fill x
		button .info.buttons.dismiss -text Dismiss \
            -default active -command "destroy .info"
		pack .info.buttons.dismiss -side left -expand 1 -pady 2
		frame .info.frame
		pack  .info.frame -expand yes -fill both -padx 1 -pady 1
		text .info.text -height 40 -wrap word\
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
    wm title .info $filename
    wm iconname .info $filename
	catch { exec bhead -verbose 7 -info $filename } theinfo
    .info.text delete 1.0 end
    .info.text insert 1.0 $theinfo
    .info.text mark set insert 1.0
}

