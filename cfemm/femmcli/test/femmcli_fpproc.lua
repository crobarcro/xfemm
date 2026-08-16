-- femmcli_fpproc.lua
-- This test is basically the same as running fmesher and fsolver on the .fem file.
-- The file femmcli_femfile.fem is the same as cfemm/fsolver/test/Temp.fem
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

open("femmcli_fpproc.fem")
mi_analyze()
mi_loadsolution()

A,B1,B2,Sig,E,H1,H2,Je,Js,Mu1,Mu2,Pe,Ph = mo_getpointvalues(0.250, 0)

-- check result against FEMM42 output:
-- FIXME: error margin needs sane values
failed=0
if getenv("XFEMM_TEST_MESHER_BACKEND") == "Tangle" then
	A_ref = 1.240205148710215e-014
	B1_ref = -1.050120917482506e-013
	B2_ref = 2.635532930443618e-014
	E_ref = 4.664092424066253e-021
	H1_ref = -8.356596743076861e-008
	H2_ref = 2.097290467807851e-008
else
	A_ref = 1.245741227364988e-014
	B1_ref = -9.855007421888915e-014
	B2_ref = 3.052725906923963e-014
	E_ref = 4.235125240802008e-021
	H1_ref = -7.842365727004682e-008
	H2_ref = 2.429282089958189e-008
end
failed = failed + check("A", A, A_ref, 2)
failed = failed + check("B1", B1, B1_ref, 2)
failed = failed + check("B2", B2, B2_ref, 2)
failed = failed + check("Sig", Sig, 0, 2)
failed = failed + check("E", E, E_ref, 3)
failed = failed + check("H1", H1, H1_ref, 2)
failed = failed + check("H2", H2, H2_ref, 2)
failed = failed + check("Je", Je, 0, 2)
failed = failed + check("Js", Js, 0, 2)
failed = failed + check("Mu1", Mu1, 1, 2)
failed = failed + check("Mu2", Mu2, 1, 2)
failed = failed + check("Pe", Pe, 0, 2)
failed = failed + check("Ph ", Ph , 0, 2)

assert(failed==0)
write("SUCCESS\n")
