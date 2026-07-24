-- femmcli_moverotate_magdir.lua
-- Regression test for GitHub issue #19: moving/rotating a group must rotate
-- the magnetic block label's numeric magnetization direction as well as its
-- position.

newdocument(0)
mi_probdef(0, "millimeters", "planar", 1e-8, 1, 30)
mi_addmaterial("TestMagnet", 1, 1, 1000000, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
mi_addblocklabel(1, 0)
mi_selectlabel(1, 0)
mi_setblockprop("TestMagnet", 1, 0, "<None>", 10, 4, 1)
mi_clearselected()

mi_selectgroup(4)
mi_moverotate(0, 0, 45, 4)
mi_saveas("femmcli_moverotate_magdir.result.fem")
