
local app         = "bsoft"
local version     = myModuleVersion()
local base = pathJoin("/usr/local/apps/",app,version)

load("gcc/10.2.0")

setenv("BSOFT", base)
setenv("BPARAM",pathJoin(base,"parameters"))
prepend_path("PATH", pathJoin(base,"bin"))
prepend_path("LD_LIBRARY_PATH", pathJoin(base,"lib"))

if (mode() == "load") then
    LmodMessage("[+] Loading  ", app, version, " ...")
end
if (mode() == "unload") then
    LmodMessage("[-] Unloading ", app, version, " ...")
end

