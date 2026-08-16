-- femmcli_TorqueBenchmark.lua
-- This checks the behaviour documented at http://www.femm.info/wiki/TorqueBenchmark/
-- Output:
-- SUCCESS
showconsole()

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
	print(result .. name .. ": " .. value .. " (expected: " .. expected
	.. ", diff: " .. diff .. " [" .. diffRel .. "%]"
		.. ", margin: " .. marginAbs .. " [" .. marginRel .. "%])")
	return fail
end

-- enable for additional output:
-- XFEMM_VERBOSE = 1

-- Reference values = analytically predicted torque
-- using absolute error margins because reference value for 0° is 0.
tq_tolerance = 0.000043 -- includes both Triangle and Tangle discretisation errors.
tq_toleranceRel = 0.006 -- worst relative error was 0.0057% at 10°
tq_ref = {}
tq_ref[0] = 0
tq_ref[10] = 0.173648
tq_ref[20] = 0.342020
tq_ref[30] = 0.5
tq_ref[40] = 0.642788
tq_ref[50] = 0.766044
tq_ref[60] = 0.866025
tq_ref[70] = 0.939693
tq_ref[80] = 0.984808
tq_ref[90] = 1

open("femmcli_TorqueBenchmark.fem")

failed=0
for deg = 0, 90, 10 do
	mi_modifyboundprop("AGE",10,deg)
	mi_modifyboundprop("AGE",11,0)
	mi_saveas("femmcli_TorqueBenchmarki_" .. deg .. ".fem")
	mi_createmesh()
	mi_analyze()
	mi_loadsolution()

	tq=mo_gapintegral("AGE", 0)
	failed= failed +check("Torque_"..deg, tq, tq_ref[deg], tq_tolerance, tq_toleranceRel)
end

assert(failed==0)
write("SUCCESS\n")
