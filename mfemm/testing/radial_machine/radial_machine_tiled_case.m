function results = radial_machine_tiled_case (varargin)
%RADIAL_MACHINE_TILED_CASE Solve the checked-in tiled-machine fixture.
%
%   results = radial_machine_tiled_case ()
%   results = radial_machine_tiled_case ('Parameter', value)
%
% Solves the tiled-magnetic model (a 60-degree rotor tile repeated six times
% and a 10-degree stator slot tile repeated 36 times, coupled only through the
% air-gap element) at the requested rotor positions and extracts the same
% observables as radial_machine_fixture_case.
%
% Parameters:
%   PositionIndices - fixture indices into data/positions.txt (default all)
%   OutputFile      - optional MAT file to which results is saved
%   Quiet           - suppress solver output (default true)

    options.OutputFile = '';
    options.Quiet = true;
    options.PositionIndices = [];
    options = parse_options (options, varargin{:});

    dataDirectory = fullfile (fileparts (mfilename ('fullpath')), 'data');
    positions = dlmread (fullfile (dataDirectory, 'positions.txt'));
    positions = positions(:);
    if isempty (options.PositionIndices)
        fixtureIndices = (1:numel (positions))';
    else
        fixtureIndices = options.PositionIndices(:);
        assert (all (fixtureIndices >= 1 & fixtureIndices <= numel (positions) ...
                     & fixtureIndices == fix (fixtureIndices)), ...
                'PositionIndices contains an invalid fixture index.');
        positions = positions(fixtureIndices);
    end

    modelFile = fullfile (dataDirectory, 'radial_machine_tiled.json');
    [coilLabels, coilOverrides, gapName, rotorLabels] = parse_tiled_model (modelFile);

    results.schemaVersion = 3;
    results.method = 'tiled';
    results.positions = positions;
    results.fluxLinkage = [];
    results.circuitFluxLinkage = [];
    results.coilFluxDensity = [];
    results.torque = [];
    results.randomA = [];

    session = xfemm.femmsession (modelFile);
    cleanup = onCleanup (@() delete (session));
    session.mesh ();
    for ind = 1:numel (positions)
        fprintf ('Tiled static solve %d of %d (fixture %d).\n', ...
                 ind, numel (positions), fixtureIndices(ind));
        session.setAGEPosition (gapName, 30 * positions(ind), 0);
        session.solve ();
        results = extract_results (results, ind, session, ...
                                   coilLabels, coilOverrides, rotorLabels);
    end

    assert (all (isfinite (results.fluxLinkage(:))));
    assert (all (isfinite (results.coilFluxDensity(:))));
    assert (all (isfinite (results.torque(:))));
    assert (all (isfinite (results.randomA(:))));
    if ~isempty (options.OutputFile)
        save (options.OutputFile, 'results', '-v7');
    end
end


function [coilLabels, coilOverrides, gapName, rotorLabels] = parse_tiled_model (modelFile)
    model = jsondecode (fileread (modelFile));
    tiles = as_array (model.tiles);
    stator = [];
    rotor = [];
    for t = 1:numel (tiles)
        if strcmp (tiles{t}.name, 'stator')
            stator = tiles{t};
        end
        if strcmp (tiles{t}.name, 'rotor')
            rotor = tiles{t};
        end
    end
    assert (~isempty (stator), 'tiled model has no stator tile');
    assert (~isempty (rotor), 'tiled model has no rotor tile');
    labels = as_array (stator.geometry.labels);
    coilLabels = {};
    for k = 1:numel (labels)
        if strncmp (labels{k}.name, 'coil', 4)
            coilLabels{end+1} = labels{k}; %#ok<AGROW>
        end
    end
    assert (numel (coilLabels) == 2, 'expected two coil labels per slot tile');
    coilOverrides = {};
    overrides = as_array (model.overrides);
    for o = 1:numel (overrides)
        coilOverrides{end+1} = overrides{o}; %#ok<AGROW>
    end
    couplings = as_array (model.couplings);
    gapName = couplings{1}.boundary;

    rotorLabels = {};
    rotorTileLabels = as_array (rotor.geometry.labels);
    for k = 1:numel (rotorTileLabels)
        material = rotorTileLabels{k}.material;
        if strcmp (material, 'NdFeB 40 MGOe') || strcmp (material, '1117 Steel')
            rotorLabels{end+1} = rotorTileLabels{k}; %#ok<AGROW>
        end
    end
end


function results = extract_results (results, positionIndex, session, ...
                                    coilLabels, coilOverrides, rotorLabels)
    nSlots = 36;
    pitchDegrees = 10;
    circuitNames = {'1', '2', '3'};
    flux = zeros (1, numel (circuitNames));
    positiveSides = zeros (1, numel (circuitNames));

    for layer = 1:numel (coilLabels)
        label = coilLabels{layer};
        override = [];
        for o = 1:numel (coilOverrides)
            if strcmp (coilOverrides{o}.label, label.name)
                override = coilOverrides{o};
            end
        end
        circuits = as_array (override.circuit);
        scales = override.turnScale;
        for slot = 1:nSlots
            [x, y] = rotate_point (label, (slot - 1) * pitchDegrees);
            turns = label.turns * scales(slot);
            session.clearblock ();
            session.selectblock (x, y);
            area = session.blockintegral (5);
            intA = session.blockintegral (1);
            circuitIndex = find (strcmp (circuitNames, circuits{slot}), 1);
            flux(circuitIndex) = flux(circuitIndex) + turns * intA / area;
            positiveSides(circuitIndex) = positiveSides(circuitIndex) + (turns > 0);
        end
    end
    results.fluxLinkage(positionIndex, :) = flux ./ positiveSides;

    % One period of the coil flux-density pattern, in the same slot-major,
    % layer-minor order as radial_machine_fixture_case.
    coilFluxDensity = zeros (1, 6 * numel (coilLabels));
    for slot = 1:6
        for layer = 1:numel (coilLabels)
            label = coilLabels{layer};
            [x, y] = rotate_point (label, (slot - 1) * pitchDegrees);
            B = session.getb (x, y);
            coilFluxDensity((slot - 1) * numel (coilLabels) + layer) = ...
                hypot (B(1), B(2));
        end
    end
    results.coilFluxDensity(positionIndex, :) = coilFluxDensity;

    % Weighted-stress-tensor torque over the rotor regions of all six
    % instances. The tiled model is the full machine, so no pole scaling.
    session.clearblock ();
    for instance = 0:5
        angle = instance * 60 * pi / 180;
        for k = 1:numel (rotorLabels)
            label = rotorLabels{k};
            x = label.x * cos (angle) - label.y * sin (angle);
            y = label.x * sin (angle) + label.y * cos (angle);
            session.selectblock (x, y, false);
        end
    end
    results.torque(positionIndex, 1) = session.blockintegral (22);

    % Random vector-potential samples in the meshed air-gap halves.
    [x, y] = radial_machine_sample_points ();
    A = session.geta (x, y);
    results.randomA(positionIndex, :) = A(:)';
end


function [x, y] = rotate_point (label, angleDegrees)
    angle = angleDegrees * pi / 180;
    x = label.x * cos (angle) - label.y * sin (angle);
    y = label.x * sin (angle) + label.y * cos (angle);
end


function v = as_array (x)
    if iscell (x)
        v = x;
    else
        v = num2cell (x);
    end
end


function options = parse_options (options, varargin)
    assert (mod (numel (varargin), 2) == 0, ...
            'Optional arguments must be parameter/value pairs.');
    names = fieldnames (options);
    for ind = 1:2:numel (varargin)
        match = find (strcmpi (varargin{ind}, names), 1);
        assert (~isempty (match), 'Unknown option: %s', varargin{ind});
        options.(names{match}) = varargin{ind+1};
    end
end
