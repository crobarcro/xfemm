function [results, design, simoptions] = radial_machine_static_case (rnfoundryRoot, varargin)
%RADIAL_MACHINE_STATIC_CASE Run the static part of the RNFoundry example.
%
% results = radial_machine_static_case (rnfoundryRoot)
% results = radial_machine_static_case (..., 'Parameter', value)
%
% This is the non-ODE part of
% example_radial_flux_permanent_magnet_machine_sim_using_femmsession.m.
% It deliberately stops after the magnetostatic position sweep and returns
% the quantities used by the regression comparisons.
%
% Parameters:
%   SolveMethod    - 'femmsession' (default) or 'xfemm_legacy'
%   RotationMethod - 'SlidingMesh' (default) or 'MagnetRedraw'
%   NPositions     - number of positions over one pole pitch (default 10)
%   OutputFile     - optional MAT file to which results is saved
%   Quiet          - suppress xfemm solver output (default true)

    options.SolveMethod = 'femmsession';
    options.RotationMethod = 'SlidingMesh';
    options.NPositions = 10;
    options.OutputFile = '';
    options.Quiet = true;
    options = parse_options (options, varargin{:});

    assert (ischar (rnfoundryRoot) && exist (rnfoundryRoot, 'dir') == 7, ...
            'XFEMM:TestSetup', 'RNFoundry root does not exist: %s', rnfoundryRoot);
    assert (isscalar (options.NPositions) && options.NPositions >= 2 ...
            && options.NPositions == fix (options.NPositions), ...
            'XFEMM:TestSetup', 'NPositions must be an integer greater than one.');

    % Keep the selected xfemm installation ahead of RNFoundry's helper
    % directories. This matters to the release-comparison runner, where each
    % solver is executed in a separate process with a different mfemm path.
    addpath (genpath (rnfoundryRoot), '-end');
    required = {'completedesign_RADIAL_SLOTTED', 'simfun_RADIAL_SLOTTED', 'fr'};
    for ind = 1:numel (required)
        assert (exist (required{ind}, 'file') ~= 0, 'XFEMM:TestSetup', ...
                'RNFoundry function %s is not available after adding %s.', ...
                required{ind}, rnfoundryRoot);
    end

    design = example_design ();
    simoptions = struct ();
    design = completedesign_RADIAL_SLOTTED (design, simoptions);
    design.Rgm = mean ([design.Rmo, design.Rai]);

    simoptions.GetVariableGapForce = false;
    simoptions.SkipInductanceFEA = true;
    simoptions.NMagFEAPositions = options.NPositions;
    simoptions.MagFEASim.SolveMethod = options.SolveMethod;
    simoptions.MagFEASim.RotationMethod = options.RotationMethod;
    simoptions.MagFEASim.UseParFor = false;
    simoptions.MagFEASim.UseFemm = false;
    simoptions.MagFEASim.QuietFemm = options.Quiet;

    [design, simoptions] = simfun_RADIAL_SLOTTED (design, simoptions);

    results.schemaVersion = 1;
    results.solveMethod = options.SolveMethod;
    results.rotationMethod = options.RotationMethod;
    results.positions = design.MagFEASimPositions(:);
    results.fluxLinkage = design.FemmDirectFluxLinkage;
    results.windingFluxLinkage = winding_flux_linkage (design);
    results.slotPositions = design.intAdata.slotPos(:);
    results.slotIntegralA = design.intAdata.slotIntA;
    results.coggingTorque = design.RawCoggingTorque(:);
    results.armatureToothFluxDensity = design.ArmatureToothFluxDensity(:);

    assert_all_finite (results.fluxLinkage, 'fluxLinkage');
    assert_all_finite (results.windingFluxLinkage, 'windingFluxLinkage');
    assert_all_finite (results.slotIntegralA, 'slotIntegralA');
    assert_all_finite (results.coggingTorque, 'coggingTorque');
    assert_all_finite (results.armatureToothFluxDensity, ...
                       'armatureToothFluxDensity');

    if ~isempty (options.OutputFile)
        save (options.OutputFile, 'results', '-v7');
    end
end


function fluxLinkage = winding_flux_linkage (design)
% Form the same gauge-invariant winding linkage used by finfun without
% running the rest of its ODE-preparation and material calculations.
    [slotPos, order] = sort (design.intAdata.slotPos);
    slotIntA = design.intAdata.slotIntA(order,:,:);
    [slotPos, uniqueIndices] = unique (slotPos);
    slotIntA = slotIntA(uniqueIndices,:,:);

    pos = slotPos(slotPos <= slotPos(1) + 2);
    intA = slotIntA(slotPos <= slotPos(1) + 2,1:2,1);
    if pos(end) < slotPos(1) + 2
        pos(end+1) = slotPos(1) + 2;
        intA(end+1,:) = interp1 (slotPos, slotIntA(:,1:2,1), pos(end));
    end

    intAslm(1) = slmengine (pos, intA(:,1), ...
                           'EndCon', 'periodic', ...
                           'knots', ceil (numel (pos) / 2), ...
                           'Plot', 'off');
    intAslm(2) = slmengine (pos, intA(:,2), ...
                           'EndCon', 'periodic', ...
                           'knots', ceil (numel (pos) / 2), ...
                           'Plot', 'off');

    coilPitch = design.thetas * design.yd / design.thetap;
    searchPositions = linspace (0, 1, 1000);
    searchFlux = fluxlinkagefrmintAslm (intAslm, coilPitch, ...
        searchPositions, design.CoilTurns, design.CoilArea, ...
        'Skew', design.MagnetSkew, ...
        'NSkewPositions', design.NSkewMagnetsPerPole);
    [~, peakIndex] = max (abs (searchFlux));

    fluxLinkage = fluxlinkagefrmintAslm (intAslm, coilPitch, ...
        linspace (0, 2, 200), design.CoilTurns, design.CoilArea, ...
        'Skew', design.MagnetSkew, ...
        'NSkewPositions', design.NSkewMagnetsPerPole, ...
        'Offset', searchPositions(peakIndex(1)));
end


function design = example_design ()
% Keep these values in step with the RNFoundry tutorial on which this test
% is based. Using a local constructor makes the test independent of figures,
% console exploration, SLM fitting, and the tutorial's ODE simulation.
    design.Poles = 12;
    design.Phases = 3;
    design.CoilLayers = 2;
    design.Qc = design.Phases * design.Poles;
    design.qc = fr (design.Qc, design.Poles * design.Phases);
    design.yd = 4;
    design.CoilFillFactor = 0.6;
    design.CoilTurns = 200;
    design.Branches = 1;

    design.Ryo = 95e-3;
    design.tm = 6.3e-3;
    design.tbi = 29.9523e-3;
    design.ty = 17.4e-3;
    design.tc = [16.4960e-3, 2.1995e-3];
    design.tsb = 2.604e-3;
    design.tsg = 1.7364e-3;
    design.g = 2e-3;
    design.thetam = (2*pi / design.Poles) * 0.667;
    design.thetacg = 84.7661e-3;
    design.thetacy = 93.0181e-3;
    design.thetasg = 38.6622e-3;
    design.ls = 88.9e-3;

    design.ArmatureType = 'external';
    design.MagnetPolarisation = 'radial';
    design.MagFEASimMaterials.AirGap = 'Air';
    design.MagFEASimMaterials.Magnet = 'NdFeB 40 MGOe';
    design.MagFEASimMaterials.FieldBackIron = '1117 Steel';
    design.MagFEASimMaterials.ArmatureYoke = ...
        design.MagFEASimMaterials.FieldBackIron;
    design.MagFEASimMaterials.ArmatureCoil = '36 AWG';
end


function options = parse_options (options, varargin)
    assert (mod (numel (varargin), 2) == 0, 'XFEMM:TestSetup', ...
            'Optional arguments must be parameter/value pairs.');
    names = fieldnames (options);
    for ind = 1:2:numel (varargin)
        match = find (strcmpi (varargin{ind}, names), 1);
        assert (~isempty (match), 'XFEMM:TestSetup', ...
                'Unknown option: %s', varargin{ind});
        options.(names{match}) = varargin{ind+1};
    end
end


function assert_all_finite (values, name)
    assert (all (isfinite (values(:))), 'XFEMM:StaticMachineResult', ...
            'Static result %s contains a non-finite value.', name);
end
