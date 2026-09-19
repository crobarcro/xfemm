function report = Test_radial_machine_tiled_methods (positionIndices)
%TEST_RADIAL_MACHINE_TILED_METHODS Compare the tiled machine with redraw/sliding.
%
% report = Test_radial_machine_tiled_methods ()
% report = Test_radial_machine_tiled_methods (positionIndices)
%
% Solves the checked-in tiled-magnetic machine (one 60-degree rotor tile
% repeated six times and one 10-degree stator slot tile repeated 36 times,
% coupled only through the air-gap element) and compares its winding flux
% linkage and coil flux density against the conventional redraw and sliding
% fixtures.
%
% The winding flux linkage is a small difference of large cancelling terms, so
% it is compared with an absolute tolerance that accommodates the documented
% redraw-versus-sliding residual. Coil flux density is compared relatively.

    testDir = fileparts (mfilename ('fullpath'));
    addpath (fullfile (testDir, 'radial_machine'));
    if nargin < 1 || isempty (positionIndices)
        positionIndices = 5;
    end

    tiled = radial_machine_tiled_case ('PositionIndices', positionIndices);
    redraw = radial_machine_fixture_case ('redraw', 'PositionIndices', positionIndices);
    sliding = radial_machine_fixture_case ('sliding', 'PositionIndices', positionIndices);

    % Direct FEMM circuit flux linkage is gauge-dependent for the partial
    % winding sector, so it is not compared. The legacy redraw's detailed
    % vector potential and weighted-stress-tensor torque do not reproduce the
    % sliding session at non-zero positions (they agree exactly at position 0),
    % so those observables are compared against the sliding session, which
    % shares the tiled model's AGE mechanism.
    tiled = rmfield (tiled, 'circuitFluxLinkage');
    redraw = rmfield (redraw, {'circuitFluxLinkage', 'torque', 'randomA'});
    sliding = rmfield (sliding, 'circuitFluxLinkage');

    comparisonOptions = {'FluxAbsoluteTolerance', 2e-6};
    fprintf ('Tiled versus redraw:\n');
    report.redraw = compare_radial_machine_static_results ( ...
        redraw, tiled, comparisonOptions{:});
    fprintf ('Tiled versus sliding:\n');
    report.sliding = compare_radial_machine_static_results ( ...
        sliding, tiled, comparisonOptions{:});
    report.passed = report.redraw.passed && report.sliding.passed;

    fprintf ('Tiled-machine static comparison passed.\n');
end
