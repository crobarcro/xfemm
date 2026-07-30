function Test_mfemmsession
%TEST_MFEMMSESSION Exercise the high-level session and native MEX gateway.
    triangle = analysis_session_uniform_field_example('Backend', 'triangle');
    assert(max(abs(triangle.B - [0, -0.25])) < eps);

    tangle = analysis_session_uniform_field_example('Backend', 'tangle');
    assert(max(abs(tangle.A - triangle.A)) < 1e-12);
    fprintf('mfemmsession analytical integration test passed.\n');
end
