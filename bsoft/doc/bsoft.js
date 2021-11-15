// Bsoft javascript library
// Created: 20110305
// Modified: 20110305
// Bernard Heymann

function changeIframeContent(id, url) {
    var iframe = document.getElementById(id);
    if ( iframe && iframe.src ) {
        iframe.src = url;
        return false;
    }
    return true;
}

// The id refers to the id tag in Safari and the name tag in Firefox
function resizeIframe(id, tx, ty) {
	var iframe = document.getElementById(id);
	if ( window.frames[id] && window.frames[id].document) {
		var body = window.frames[id].document.body;
		if ( body ) {
			var obj = (iframe.style) ? iframe.style : iframe;
			obj.width = window.innerWidth - tx + 'px';
			obj.height = window.innerHeight - ty + 'px';
		} else {
			alert("body not found!");
		}
	} else {
		alert("iframe " + id + " not found!");
	}
}

function setIFrameSize(id) {
    f = $("ifrm2");
    f[0].setAttribute("width", f.parent().width());
    f[0].setAttribute("height", f.parent().height());
}
