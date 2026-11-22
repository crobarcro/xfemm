function Test_femmsession_tiled
%TEST_FEMMSESSION_TILED Smoke-test the tiled-magnetic JSON session loader.
%   Loads a tiled model, meshes each tile once, repeats it, and solves.

    repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    modelFile = fullfile(repositoryRoot, 'mfemm', 'testing', 'tiled_ring.json');
    session = xfemm.femmsession(modelFile);
    cleanup = onCleanup(@() delete(session)); %#ok<NASGU>

    elementCount = session.mesh();
    assert(elementCount > 0, 'tiled model produced an empty mesh');

    info = session.instancedInfo();
    assert(info.instanceCount == 36, 'tiled model did not create 36 instances');
    assert(info.templateTopologyIdentity ~= 0);
    assert(info.materializationCount == 1);

    session.setCircuit('phase', 'current', 1);
    status = session.solve();
    assert(status.success, 'tiled model session solve failed');

    % A second solve reuses the materialised topology.
    second = session.solve();
    assert(second.success);
    assert(second.meshGenerationCount == 1 && second.solveCount == 2);
    info = session.instancedInfo();
    assert(info.materializationCount == 1);

    fprintf('femmsession tiled model smoke test passed.\n');
end
