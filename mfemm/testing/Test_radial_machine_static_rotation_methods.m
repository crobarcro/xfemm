function report = Test_radial_machine_static_rotation_methods (positionIndices)
%TEST_RADIAL_MACHINE_STATIC_ROTATION_METHODS Compare AGE and remeshing.
%
% report = Test_radial_machine_static_rotation_methods ()
% report = Test_radial_machine_static_rotation_methods (positionIndices)
%
% Runs the static portion of the RNFoundry radial-machine example twice
% with the current xfemm: first using the in-memory Air Gap Element session,
% then by redrawing and remeshing the magnets at every position.

    testDir = fileparts (mfilename ('fullpath'));
    addpath (fullfile (testDir, 'radial_machine'));
    if nargin < 1
        positionIndices = [];
    end

    sliding = radial_machine_fixture_case (...
        'sliding', 'PositionIndices', positionIndices);
    redraw = radial_machine_fixture_case (...
        'redraw', 'PositionIndices', positionIndices);

    % A is defined up to a constant in the periodic sector. Direct circuit
    % linkage is therefore gauge-dependent for the partial winding model;
    % fluxLinkage below is formed from opposing coil sides and is invariant.
    % The legacy redraw's detailed vector potential and torque integrals do not
    % reproduce the sliding session at non-zero positions, so those observables
    % are reported only in the tiled-machine and sliding-versus-redraw tests.
    sliding = rmfield (sliding, {'circuitFluxLinkage', 'torque', 'airgapTorque', 'randomA'});
    redraw = rmfield (redraw, {'circuitFluxLinkage', 'torque', 'airgapTorque', 'randomA'});

    report = compare_radial_machine_static_results (redraw, sliding, ...
        'FluxAbsoluteTolerance', 2e-6);
    fprintf ('Sliding-mesh versus redraw static machine test passed.\n');
end
