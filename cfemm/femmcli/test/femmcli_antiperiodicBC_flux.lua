-- femmcli_antiperiodicBC_flux.lua
-- This checks the correct flux density calculation for a motor with antiperiodic BC.
-- (The calculation contained errors in xfemm 2.0)
-- Output:
-- SUCCE

MAX_DIFF = 0
MAX_DIFF_REL = 0
-- check variable <name>,
-- compare <value> against <expected> value
-- if the absolute or relative difference is greater than the margin, complain and return 1
-- if the expected value is 0, the relative margin is ignored
-- relative margin is in percent
function check(name, value, expected, marginAbs, marginRel)
   diff=value - expected
   diffRel=0
   if (expected~=0) then
      diffRel=100*diff/expected
   end
   if abs(diff) > marginAbs or abs(diffRel) > marginRel then
      fail=1
      result="[FAILED] "
   else
      fail=0
      result="[  ok  ] "
   end
	if abs(diff) > MAX_DIFF then
		MAX_DIFF = abs(diff)
	end
	if abs(diffRel) > MAX_DIFF_REL then
		MAX_DIFF_REL = abs(diffRel)
	end
   print(result .. name .. ": " .. value .. " (expected: " .. expected
   .. ", diff: " .. diff .. " [" .. diffRel .. "%]"
      .. ", margin: " .. marginAbs .. " [" .. marginRel .. "%])")
   return fail
end


show_console()
open("femmcli_antiperiodicBC_flux.fem")

mi_saveas("femmcli_antiperiodicBC_flux.result.fem")
mi_analyze(1)
mi_loadsolution()

tolerance = 0.021
tolerance_rel = 70

-- reference values from FEMM42:
B_abs_ref = {}
B_abs_ref[0] = 2.0172349211547e-005
B_abs_ref[1] = 9.980515629468256e-005
B_abs_ref[2] = 0.0001299999163508703
B_abs_ref[3] = 0.0001123790790807998
B_abs_ref[4] = 0.0001035638039144443
B_abs_ref[5] = 4.712390174349349e-005
B_abs_ref[6] = 1.52445645287574e-005
B_abs_ref[7] = 2.238963159664765e-005
B_abs_ref[8] = 5.255076860801661e-005
B_abs_ref[9] = 0.0001501846812706276
B_abs_ref[10] = 0.0001875769644342095
B_abs_ref[11] = 0.3434866174139751
B_abs_ref[12] = 0.3134173948888234
B_abs_ref[13] = 0.2044994410897019
B_abs_ref[14] = 0.007451582859305774
B_abs_ref[15] = 0.0008821347188708024
B_abs_ref[16] = 6.601243412719332e-005
B_abs_ref[17] = 7.025226578942951e-005
B_abs_ref[18] = 0.4604507368660987
B_abs_ref[19] = 0.7077796643063485
B_abs_ref[20] = 0.000977829137071488
B_abs_ref[21] = 0.001296664349720595
B_abs_ref[22] = 0.3067033815500693
B_abs_ref[23] = 0.001112765134440902
B_abs_ref[24] = 0.001093521065145376
B_abs_ref[25] = 0.1565090701118708
B_abs_ref[26] = 0.2693002040535217
B_abs_ref[27] = 0.000862706441337431
B_abs_ref[28] = 0.9772489557271694
B_abs_ref[29] = 0.00123522025689945
B_abs_ref[30] = 0.01587650688903416
B_abs_ref[31] = 0.6762699821164218
B_abs_ref[32] = 0.01412679651697672
B_abs_ref[33] = 0.002241449307702012
B_abs_ref[34] = 0.348186232581696
B_abs_ref[35] = 0.0007352310639869735
B_abs_ref[36] = 0.002163325834308883
B_abs_ref[37] = 0.1874576095909183
B_abs_ref[38] = 1.162091824546714
B_abs_ref[39] = 1.0038388444812
B_abs_ref[40] = 0.2763546286608281
B_abs_ref[41] = 1.001015575730773
B_abs_ref[42] = 2.847901526518666
B_abs_ref[43] = 0.7007337071437815
B_abs_ref[44] = 0.00312702704756613

failed=0
idx = 0
for x = -40, -20, 5 do
	for y = -20, 20, 5 do
		A,Bx,By,Sig,E,H1,H2,Je,Js,Mu1,Mu2,Pe,Ph = mo_getpointvalues(x,y)
		B_abs = abs(Bx) + abs(By)
		--write("B_abs_ref["..idx.."] = " .. B_abs)
		failed= failed +check("|B| @ "..x..", "..y, B_abs, B_abs_ref[idx], tolerance, tolerance_rel)
		idx = idx + 1
	end
end
print("Maximum difference (abs,rel) = " .. MAX_DIFF ..", " ..MAX_DIFF_REL)
print("Maximum tolerated difference (abs,rel) = " .. tolerance .. ", " .. tolerance_rel)

assert(failed==0)
write("SUCCESS\n")
quit()
