function result = analysis_session_uniform_field_example(varargin)
%ANALYSIS_SESSION_UNIFORM_FIELD_EXAMPLE Verify FEMMSESSION analytically.
%   RESULT = ANALYSIS_SESSION_UNIFORM_FIELD_EXAMPLE creates a one metre
%   square of air with vector potential A = B0*x prescribed along its
%   boundary. Because that potential is linear, the finite-element solution
%   is exactly A = B0*x on any conforming first-order mesh and the magnetic flux
%   density is the uniform field B = [0,-B0]. The function solves through
%   FEMMSESSION and asserts that the native nodal result agrees with the
%   analytical vector potential.
%
%   ...('Backend', NAME) selects 'triangle' (default) or 'tangle'.
%   ...('KeepProblem', true) retains the generated .fem file and returns its
%   path in RESULT.problemFile. This is useful for inspecting the example.
%
%   This example is also the MATLAB and GNU Octave MEX integration test.

    options.Backend = 'triangle';
    options.KeepProblem = false;
    options = mfemmdeps.parseoptions(options, varargin);

    B0 = 0.25; % tesla
    problem = newproblem_mfemm('planar', 'Frequency', 0, ...
                               'LengthUnits', 'meters');
    [problem, ~, nodeids] = addnodes_mfemm(problem, [0, 1, 1, 0], ...
                                                   [0, 0, 1, 1]);
    [problem, ~, boundaryName] = addboundaryprop_mfemm(problem, ...
                                         'UniformField', 0, 'A1', B0);
    problem = addsegments_mfemm(problem, nodeids, ...
                                nodeids([2, 3, 4, 1]), ...
                                'BoundaryMarker', boundaryName);
    problem = addblocklabel_mfemm(problem, 0.5, 0.5, ...
                                  'BlockType', 'Air', 'MaxArea', 1);

    problemFile = [tempname(), '.fem'];
    writefemmfile(problemFile, problem);
    cleanup = onCleanup(@() deleteProblem(problemFile, options.KeepProblem));

    session = xfemm.femmsession(problemFile);
    session.setBackend(options.Backend);
    elementCount = session.mesh();
    status = session.solve();
    assert(status.success && status.elementCount == elementCount);
    assert(~isfield(status, 'A') && ~isfield(status, 'x') && ~isfield(status, 'y'), ...
           'MFEMM:session:unexpectedFieldData', ...
           'SOLVE status must not copy the nodal field solution.');

    sampleX = [0.25, 0.5, 0.75];
    sampleA = session.geta(sampleX, 0.5 * ones(size(sampleX)));
    assert(max(abs(sampleA - B0 * sampleX)) < 1e-7, ...
           'MFEMM:session:analyticalMismatch', ...
           'Session result does not match A = B0*x.');

    % A solved session exposes the complete fpproc interface directly.
    assert(session.nummeshnodes() == status.nodeCount);
    B = session.getb(0.5, 0.5);
    assert(max(abs(B - [0; -B0])) < 1e-7, ...
           'MFEMM:session:postProcessorMismatch', ...
           'Session post-processing does not match the analytical field.');

    session.accept();
    stateFile = [tempname(), '.mat'];
    stateCleanup = onCleanup(@() deleteIfPresent(stateFile));
    session.saveState(stateFile);
    session.reject();
    session.loadState(stateFile);

    result = status;
    result.B = [0, -B0];
    result.elementCount = elementCount;
    if options.KeepProblem, result.problemFile = problemFile; end
    clear cleanup stateCleanup session;
end

function deleteProblem(filename, keep)
    if ~keep, deleteIfPresent(filename); end
end

function deleteIfPresent(filename)
    if exist(filename, 'file'), delete(filename); end
end
