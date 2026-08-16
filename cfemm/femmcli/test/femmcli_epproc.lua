-- femmcli_epproc.lua
-- The file femmcli_epproc.fee is the same as cfemm/esolver/test/test.fee
-- SUCCESS
showconsole()

-- check variable <name>,
-- compare <value> against <expected> value
-- if the relative difference is greater than the margin (in percent), complain and return 1
function check(name, value, expected, margin)
	diff=100*(value - expected) / expected
	if abs(diff) > margin then
		fail=1
		result="[FAILED] "
	else
		fail=0
		result="[  ok  ] "
	end
	print(result .. name .. ": " .. value .. " (expected: " .. expected .. ", diff: " .. diff .. "%, margin: " .. margin .. "%)")
	--print("check(\""..name.."\", "..name..", "..value..", "..margin..")")
	return fail
end

-- enable for additional output:
-- XFEMM_VERBOSE = 1

open("femmcli_epproc.fee")
ei_analyze(0)
ei_loadsolution()

V,Dx,Dy,Ex,Ey,ex,ey,nrg  = eo_getpointvalues(0.250, 0)

-- check result against FEMM42 output:
-- FIXME: error margin needs sane values
failed=0
if getenv("XFEMM_TEST_MESHER_BACKEND") == "Tangle" then
	V_ref = 48.90577131945231
	Dx_ref = 8.387310310792659e-010
	Dy_ref = 6.043844822898733e-011
	Ex_ref = 23.68176077680935
	Ey_ref = 1.706493285265354
	nrg_ref = 9.982882720090222e-009
else
	V_ref = 48.37056814422403
	Dx_ref = 1.157764975200258e-009
	Dy_ref = 7.559208128960357e-011
	Ex_ref = 32.68975650415626
	Ey_ref = 2.134359549590023
	nrg_ref = 1.900419790445539e-008
end
failed= failed +check("V", V, V_ref, 1)
failed= failed +check("Dx", Dx, Dx_ref, 1.5)
failed= failed +check("Dy", Dy, Dy_ref, 15)
failed= failed +check("Ex", Ex, Ex_ref, 2)
failed= failed +check("Ey", Ey, Ey_ref, 15)
failed= failed +check("ex", ex, 4, 0.1)
failed= failed +check("ey", ey, 4, 0.1)
failed= failed +check("nrg", nrg, nrg_ref, 3)

assert(failed==0)
write("SUCCESS\n")
