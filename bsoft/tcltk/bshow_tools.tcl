##
# @file		bshow_tools.tcl
#
# @brief	The toolbox window
#
# @author	Bernard Heymann
# @date		Created: 20051130
# @date		Modified: 20130924

set tool point

## @brief Toolbox to manage overlay objects
#

proc Toolbox { } {
	global helv12
	global tool
	global bshow_script
	
	if ![winfo exists .wtool] {
		toplevel .wtool
		
		image create photo point_off -file $bshow_script/point.gif
		image create photo point_on -file $bshow_script/point_active.gif
		radiobutton .wtool.point -image point_off -selectimage point_on \
			-relief raised -indicatoron 0 -variable tool -value point
#		ttk::radiobutton .wtool.point -image { point_off selected point_on } \
#			-compound image -variable tool -value point
		
		image create photo voxels_off -file $bshow_script/voxels.gif
		image create photo voxels_on -file $bshow_script/voxels_active.gif
		radiobutton .wtool.voxels -image voxels_off -selectimage voxels_on \
			-relief raised -indicatoron 0 -variable tool -value voxels -command { Voxels }
		
		image create photo select_off -file $bshow_script/select.gif
		image create photo select_on -file $bshow_script/select_active.gif
		radiobutton .wtool.select -image select_off -selectimage select_on \
			-relief raised -indicatoron 0 -variable tool -value select -command { Selection }
		
		image create photo measure_off -file $bshow_script/measure.gif
		image create photo measure_on -file $bshow_script/measure_active.gif
		radiobutton .wtool.measure -image measure_off -selectimage measure_on \
			-relief raised -indicatoron 0 -variable tool -value measure -command { Measure }
		
		image create photo box_off -file $bshow_script/box.gif
		image create photo box_on -file $bshow_script/box_active.gif
		radiobutton .wtool.box -image box_off -selectimage box_on \
			-relief raised -indicatoron 0 -variable tool -value box -command { Particles }
		
		image create photo filament_off -file $bshow_script/filament.gif
		image create photo filament_on -file $bshow_script/filament_active.gif
		radiobutton .wtool.filament -image filament_off -selectimage filament_on \
			-relief raised -indicatoron 0 -variable tool -value filament -command { Filaments }
		
		image create photo helix_off -file $bshow_script/helix.gif
		image create photo helix_on -file $bshow_script/helix_active.gif
		radiobutton .wtool.helix -image helix_off -selectimage helix_on \
			-relief raised -indicatoron 0 -variable tool -value helix -command { Helix }
		
		image create photo xtal_off -file $bshow_script/xtal.gif
		image create photo xtal_on -file $bshow_script/xtal_active.gif
		radiobutton .wtool.xtal -image xtal_off -selectimage xtal_on \
			-relief raised -indicatoron 0 -variable tool -value xtal -command { Crystallography }
		
		image create photo marker_off -file $bshow_script/marker.gif
		image create photo marker_on -file $bshow_script/marker_active.gif
		radiobutton .wtool.marker -image marker_off -selectimage marker_on \
			-relief raised -indicatoron 0 -variable tool -value marker -command { Tomography }
		
		image create photo marker_group_off -file $bshow_script/marker_group.gif
		image create photo marker_group_on -file $bshow_script/marker_group_active.gif
		radiobutton .wtool.marker_group -image marker_group_off -selectimage marker_group_on \
			-relief raised -indicatoron 0 -variable tool -value marker_group -command { Tomography }
		
		image create photo model_off -file $bshow_script/model.gif
		image create photo model_on -file $bshow_script/model_active.gif
		radiobutton .wtool.model -image model_off -selectimage model_on \
			-relief raised -indicatoron 0 -variable tool -value model -command { Model }
		
		image create photo toolhelp_off -file $bshow_script/help.gif
		image create photo toolhelp_on -file $bshow_script/help_active.gif
		radiobutton .wtool.toolhelp -image toolhelp_off -selectimage toolhelp_on \
			-relief raised -indicatoron 0 -variable tool -value toolhelp \
				-command { showURL "file://$Bsoft/doc/bshow/bshow_tools.html" }
#				-command { showHelp bshow_tools.hlp }
		
		set tool point
		
		pack .wtool.point .wtool.voxels .wtool.select .wtool.measure \
			.wtool.box .wtool.filament .wtool.helix .wtool.xtal \
			.wtool.marker .wtool.marker_group .wtool.model \
			.wtool.toolhelp -side top -padx 4 -pady 4
		
	} else {
		wm deiconify .wtool
		raise .wtool
    }
	
    wm title .wtool "Toolbox"
    wm iconname .wtool "Toolbox"
	wm resizable .wtool no yes
	wm protocol .wtool WM_DELETE_WINDOW return
#	wm overrideredirect .wtool 1

	toolTip .wtool.point "Select objects"
	toolTip .wtool.voxels "Record voxel values"
	toolTip .wtool.select "Select a region"
	toolTip .wtool.measure "Measure distances"
	toolTip .wtool.box "Pick particles"
	toolTip .wtool.filament "Pick filaments"
	toolTip .wtool.helix "Pick layer lines"
	toolTip .wtool.xtal "Pick reflections"
	toolTip .wtool.marker "Pick fiducial markers"
	toolTip .wtool.marker_group "Select all fiducial markers in current image"
	toolTip .wtool.model "Pick components and links"
	toolTip .wtool.toolhelp "Display help"
}


