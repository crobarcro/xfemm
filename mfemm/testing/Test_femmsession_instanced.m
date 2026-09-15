function Test_femmsession_instanced
%TEST_FEMMSESSION_INSTANCED Smoke-test rotational instancing through the MEX gateway.
%   Meshes one tile once, repeats it eight times about the origin, solves the
%   materialised ring, and checks that a repeated solve reuses the topology.

    repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    tile = fullfile(repositoryRoot, 'cfemm', 'fpproc', 'test', 'fixtures', ...
                    'd_ring_tile.fem');
    session = xfemm.femmsession(tile);
    cleanup = onCleanup(@() delete(session)); %#ok<NASGU>

    session.setRotationalInstances(0, 0, 8, 360);
    elementCount = session.mesh();
    info = session.instancedInfo();
    assert(info.instanceCount == 8);
    assert(info.templateTopologyIdentity ~= 0);
    assert(info.instanceLayoutIdentity ~= 0);
    assert(info.materializationCount == 1);
    assert(elementCount > 0);

    status = session.solve();
    assert(status.success);
    assert(info.materializationCount == 1);
    assert(status.elementCount > 0);
    assert(session.nummeshnodes() == status.nodeCount);

    % A second solve reuses the materialised topology.
    second = session.solve();
    assert(second.success);
    assert(second.meshGenerationCount == 1 && second.solveCount == 2);
    info = session.instancedInfo();
    assert(info.materializationCount == 1);

    fprintf('femmsession rotational instancing smoke test passed.\n');
end
