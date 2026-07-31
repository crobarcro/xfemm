function Test_mfemmsession
%TEST_MFEMMSESSION Exercise the high-level session and native MEX gateway.
    triangle = analysis_session_uniform_field_example('Backend', 'triangle');
    assert(max(abs(triangle.B - [0, -0.25])) < eps);

    tangle = analysis_session_uniform_field_example('Backend', 'tangle');
    assert(tangle.success && triangle.success);
    assert(max(abs(tangle.B - triangle.B)) < eps);
    fprintf('mfemmsession analytical integration test passed.\n');
end
