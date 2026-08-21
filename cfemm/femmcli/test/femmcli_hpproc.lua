-- femmcli_hpproc.lua
-- This test is basically the same as running fmesher and fsolver on the .fem file.
-- The file femmcli_femfile.feh is the same as cfemm/hsolver/test/Temp0.feh
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

open("femmcli_hpproc.feh")
hi_analyze()
hi_loadsolution()

T,Fx,Fy,Gx,Gy,kx,ky= ho_getpointvalues(1.1,1.1)

failed=0
-- check result against FEMM42 output:
-- FIXME: error margin needs sane values
print("Checks against femm42 output:")
if getenv("XFEMM_MESHER_BACKEND") == "Tangle" then
	T_ref = 304.8996961164884
	Fx_ref = 0.07573299004184463
	Fy_ref = 0.07330545755661612
	Gx_ref = 2.862950337283372
	Gy_ref = 2.771181810205393
else
	T_ref = 304.8641290114103
	Fx_ref = 0.2199070927061962
	Fy_ref = 0.1428113935654898
	Gx_ref = 8.313999477015031
	Gy_ref = 5.399252187839117
end
failed = failed + check("T", T, T_ref, 2)
failed = failed + check("Fx", Fx, Fx_ref, 4)
failed = failed + check("Fy", Fy, Fy_ref, 4)
failed = failed + check("Gx", Gx, Gx_ref, 4)
failed = failed + check("Gy", Gy, Gy_ref, 4)
failed = failed + check("kx", kx, 0.02645021728882154, 2)
failed = failed + check("ky", ky, 0.02645021728882154, 2)

assert(failed==0)
write("SUCCESS\n")
