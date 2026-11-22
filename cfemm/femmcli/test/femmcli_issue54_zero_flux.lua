-- Regression test for GitHub issue #54.
-- A copied linear material must keep MuMax == 0 so loading a zero-flux
-- magnetostatic solution does not enter the incremental-permeability path.

newdocument(0)
mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)
mi_getmaterial("Air")

mi_addnode(0, 1)
mi_addnode(0, -1)
mi_addsegment(0, -1, 0, 1)
mi_addarc(0, -1, 0, 1, 180, 100)

mi_addblocklabel(0.5, 0)
mi_selectlabel(0.5, 0)
mi_setblockprop("Air", 1, 0, "", 0, 0, 0)
mi_clearselected()

mi_saveas("femmcli_issue54_zero_flux.fem")
mi_analyze()
mi_loadsolution()

write("SUCCESS\n")
