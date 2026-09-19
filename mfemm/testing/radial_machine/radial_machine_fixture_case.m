function results = radial_machine_fixture_case (method, varargin)
%RADIAL_MACHINE_FIXTURE_CASE Solve checked-in static machine fixtures.
%
% method is 'sliding' for one in-memory AGE session or 'redraw' for the ten
% independently drawn and remeshed FEM problems.

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

    results.schemaVersion = 3;
    results.method = method;
    results.positions = positions;
    results.fluxLinkage = [];
    results.circuitFluxLinkage = [];
    results.coilFluxDensity = [];
    results.torque = [];
    results.randomA = [];

    switch lower (method)
        case 'sliding'
            filename = fullfile (dataDirectory, 'radial_machine_sliding.fem');
            problem = loadfemmfile (filename);
            ageIndices = find ([problem.BoundaryProps.BdryType] == 6 ...
                               | [problem.BoundaryProps.BdryType] == 7);
            assert (numel (ageIndices) == 1, ...
                    'Sliding fixture must contain exactly one AGE property.');
            ageName = problem.BoundaryProps(ageIndices).Name;
            session = xfemm.femmsession (filename);
            cleanup = onCleanup (@() delete (session));
            for ind = 1:numel (positions)
                fprintf ('Sliding static solve %d of %d (fixture %d).\n', ...
                         ind, numel (positions), fixtureIndices(ind));
                session.setAGEPosition (ageName, 30 * positions(ind), 0);
                session.solve ();
                results = extract_results (results, ind, problem, session);
            end

        case 'redraw'
            for ind = 1:numel (positions)
                fprintf ('Redraw static solve %d of %d (fixture %d).\n', ...
                         ind, numel (positions), fixtureIndices(ind));
                source = fullfile (dataDirectory, ...
                    sprintf ('radial_machine_redraw_%02d.fem', fixtureIndices(ind)));
                problem = loadfemmfile (source);
                [answerFile, problemFile] = analyse_mfemm (...
                    problem, 'Quiet', options.Quiet, 'KeepMesh', false);
                cleanupFiles = onCleanup (@() delete_files (answerFile, problemFile));
                solution = fpproc (answerFile);
                results = extract_results (results, ind, problem, solution);
                clear solution cleanupFiles;
            end

        otherwise
            error ('XFEMM:TestSetup', ...
                   'Method must be ''sliding'' or ''redraw'' (got %s).', method);
    end

    assert (all (isfinite (results.fluxLinkage(:))));
    assert (all (isfinite (results.circuitFluxLinkage(:))));
    assert (all (isfinite (results.coilFluxDensity(:))));
    assert (all (isfinite (results.torque(:))));
    assert (all (isfinite (results.randomA(:))));
    if ~isempty (options.OutputFile)
        save (options.OutputFile, 'results', '-v7');
    end
end


function results = extract_results (results, positionIndex, problem, solution)
    circuitNames = {problem.Circuits.Name};
    for circuitIndex = 1:numel (circuitNames)
        props = solution.getcircuitprops (circuitNames{circuitIndex});
        results.circuitFluxLinkage(positionIndex, circuitIndex) = props(3);

        labelIndices = find (strcmp ({problem.BlockLabels.InCircuit}, ...
                                    circuitNames{circuitIndex}));
        signedFlux = 0;
        positiveSides = 0;
        for labelIndex = labelIndices
            label = problem.BlockLabels(labelIndex);
            if label.Turns == 0
                continue;
            end
            solution.clearblock ();
            solution.selectblock (label.Coords(1), label.Coords(2));
            area = solution.blockintegral (5);
            intA = solution.blockintegral (1);
            signedFlux = signedFlux + label.Turns * intA / area;
            positiveSides = positiveSides + (label.Turns > 0);
        end
        assert (positiveSides > 0, 'Circuit %s has no positive coil side.', ...
                circuitNames{circuitIndex});
        results.fluxLinkage(positionIndex, circuitIndex) = ...
            signedFlux / positiveSides;
    end

    coilIndices = [];
    for circuitIndex = 1:numel (circuitNames)
        coilIndices = [coilIndices, ...
            find(strcmp({problem.BlockLabels.InCircuit}, ...
                        circuitNames{circuitIndex}))];
    end
    coilIndices = sort (coilIndices);
    coords = reshape ([problem.BlockLabels(coilIndices).Coords], 2, [])';
    B = solution.getb (coords(:,1), coords(:,2));
    results.coilFluxDensity(positionIndex, :) = sqrt (sum (B.^2, 1));

    % Weighted-stress-tensor torque on the rotor, normalised to the full
    % machine, matching feasim_RADIAL_SLOTTED.
    rotorIndices = rotor_label_indices (problem);
    solution.clearblock ();
    for labelIndex = rotorIndices
        label = problem.BlockLabels(labelIndex);
        solution.selectblock (label.Coords(1), label.Coords(2), false);
    end
    results.torque(positionIndex, 1) = ...
        12 * solution.blockintegral (22) / 2;

    % Random vector-potential samples in the meshed air-gap halves.
    [x, y] = radial_machine_sample_points ();
    A = solution.geta (x, y);
    results.randomA(positionIndex, :) = A(:)';
end


function indices = rotor_label_indices (problem)
    indices = [];
    for labelIndex = 1:numel (problem.BlockLabels)
        label = problem.BlockLabels(labelIndex);
        radius = hypot (label.Coords(1), label.Coords(2));
        if strcmp (label.BlockType, 'NdFeB 40 MGOe') ...
                || (strcmp (label.BlockType, '1117 Steel') && radius < 0.0565)
            indices(end+1) = labelIndex; %#ok<AGROW>
        end
    end
end


function delete_files (varargin)
    for ind = 1:nargin
        if exist (varargin{ind}, 'file') == 2
            delete (varargin{ind});
        end
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
