function Test_femmsession
%TEST_FEMMSESSION Exercise the high-level session and native MEX gateway.
    triangle = analysis_session_uniform_field_example('Backend', 'triangle');
    assert(max(abs(triangle.B - [0, -0.25])) < eps);

    tangle = analysis_session_uniform_field_example('Backend', 'tangle');
    assert(tangle.success && triangle.success);
    assert(max(abs(tangle.B - triangle.B)) < eps);

    % A session owns the parsed model and must remain usable after a caller
    % removes the temporary .fem used to construct it. Repeated solves must
    % reuse the mesh and refresh post-processing without creating an .ans file.
    repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    source = fullfile(repositoryRoot, ...
                      'cfemm', 'fsolver', 'test', 'Temp.fem');
    scratchBase = tempname();
    scratchModel = [scratchBase, '.fem'];
    copyfile(source, scratchModel);
    session = xfemm.femmsession(scratchModel);
    elementCount = session.mesh();
    delete(scratchModel);

    first = session.solve();
    second = session.solve();
    assert(first.success && second.success);
    assert(second.meshGenerationCount == 1 && second.solveCount == 2);
    assert(second.elementCount == elementCount);
    assert(session.nummeshnodes() == second.nodeCount);
    assert(exist([scratchBase, '.ans'], 'file') == 0);
    delete(session);

    % Invalid nonlinear material data must become a normal MATLAB error;
    % no C++ exception may escape the MEX boundary.
    invalidModel = [tempname(), '.fem'];
    text = fileread(source);
    text = regexprep(text, '<BHPoints>\s*=\s*9', '<BHPoints> = 1', 'once');
    fid = fopen(invalidModel, 'w');
    assert(fid ~= -1);
    cleanup = onCleanup(@() deleteIfPresent(invalidModel));
    fwrite(fid, text);
    fclose(fid);
    try
        invalidSession = xfemm.femmsession(invalidModel); %#ok<NASGU>
        error('Test_femmsession:ExpectedError', ...
              'Invalid B-H curve was unexpectedly accepted.');
    catch exception
        assert(~strcmp(exception.identifier, 'Test_femmsession:ExpectedError'));
        assert(~isempty(strfind(exception.message, ...
                               'zero points or at least two'))); %#ok<STREMP>
    end

    fprintf('femmsession analytical integration test passed.\n');
end

function deleteIfPresent(filename)
    if exist(filename, 'file') == 2
        delete(filename);
    end
end
