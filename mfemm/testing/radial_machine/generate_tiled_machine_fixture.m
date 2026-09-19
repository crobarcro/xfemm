function generate_tiled_machine_fixture (rnfoundryRoot, outputFile)
%GENERATE_TILED_MACHINE_FIXTURE Write a tiled-machine JSON from RNFoundry.
%
%   generate_tiled_machine_fixture (rnfoundryRoot, outputFile)
%
% RNFoundry is a fixture-generation dependency only; ordinary CI consumes the
% checked-in JSON and never needs RNFoundry. The model is the 12-pole/36-slot
% design of radial_machine_static_case, decomposed into one 60-degree
% pole-pair rotor tile (6 instances) and one 10-degree slot tile (36
% instances), coupled only through the air-gap element. The air gap is split
% at a third and two thirds of the gap, matching the AGE arc radii used by the
% full-sector sliding model, so each tile carries one smooth AGE arc.

    assert (ischar (rnfoundryRoot) && exist (rnfoundryRoot, 'dir') == 7, ...
            'RNFoundry root does not exist: %s', rnfoundryRoot);
    addpath (genpath (rnfoundryRoot), '-end');

    % RNFoundry assigns random group numbers; seed for a reproducible fixture.
    rng (0);

    design = example_design ();
    design = completedesign_RADIAL_SLOTTED (design, struct ());
    design.Rgm = mean ([design.Rmo, design.Rai]);

    library = fullfile (fileparts (which ('matstr2matstruct_mfemm')), '..', 'matlib.mat');
    materialNames = { design.MagFEASimMaterials.AirGap, ...
                      design.MagFEASimMaterials.Magnet, ...
                      design.MagFEASimMaterials.FieldBackIron, ...
                      design.MagFEASimMaterials.ArmatureYoke, ...
                      design.MagFEASimMaterials.ArmatureCoil };

    rotor = build_rotor_tile (design, library, materialNames);
    stator = build_stator_tile (design, library, materialNames);
    seamName = rotor.BoundaryProps(1).Name;
    assert (strcmp (seamName, stator.BoundaryProps(1).Name), ...
            'tiles disagree on the seam boundary name');

    [materials, ~] = merge_materials ({rotor, stator});
    [boundaries, ~] = merge_boundaries ({rotor, stator});
    [circuits, ~] = merge_circuits ({rotor, stator});

    model = struct ();
    model.xfemm = struct ('format', 'tiled-magnetic', 'version', 1);
    model.problem = struct ('type', 'planar', 'coords', 'cartesian', ...
                            'units', 'meters', 'depth', design.ls, ...
                            'frequency', 0, 'precision', 1e-8, 'minAngle', 30);
    model.pointProps = {};
    model.boundaries = boundaries;
    model.materials = materials;
    model.circuits = circuits;

    gapName = rotor.BoundaryProps(2).Name;
    assert (strcmp (gapName, stator.BoundaryProps(2).Name), ...
            'tiles disagree on the air-gap boundary name');
    rotorTile = tile_from_problem (rotor, 'rotor', 6, seamName, gapName);
    statorTile = tile_from_problem (stator, 'stator', 36, seamName, gapName);
    model.tiles = {rotorTile, statorTile};

    gapInner = design.Rmo + design.g/3;
    gapOuter = design.Rmo + 2*design.g/3;
    coupling = struct ('inner', 'rotor', 'outer', 'stator', 'boundary', gapName, ...
                       'center', [0, 0], 'innerRadius', gapInner, ...
                       'outerRadius', gapOuter);
    model.couplings = {coupling};
    model.overrides = stator_phase_overrides (design);

    text = jsonencode (model);
    fid = fopen (outputFile, 'w');
    assert (fid > 0, 'could not open %s for writing', outputFile);
    fwrite (fid, text);
    fclose (fid);
    fprintf ('Wrote tiled machine fixture to %s\n', outputFile);
end


function problem = build_rotor_tile (design, library, materialNames)
% Draw one 60-degree pole-pair rotor and extend the air gap to the first AGE
% arc. The rotor tile naturally contains the N/S magnet pair.
    problem = newproblem_mfemm ('planar', 'Depth', design.ls, 'MinAngle', 30);
    [problem, matinds] = addmaterials_mfemm (problem, materialNames, ...
                                             'MaterialsLibrary', library);
    [problem, ~, ~] = addboundaryprop_mfemm (problem, 'seam', 4);
    [problem, ~, ~] = addboundaryprop_mfemm (problem, 'gap', 6);
    seamName = problem.BoundaryProps(end-1).Name;
    gapName = problem.BoundaryProps(end).Name;

    magnetMesh = choosemesharea_mfemm (design.tm, design.Rmm * design.thetam, 1/10);
    backIronMesh = choosemesharea_mfemm (min (design.tbi), 2 * design.Rbm * design.thetap, 1/10);
    outerMesh = [choosemesharea_mfemm(design.tm, design.Rbo * design.thetap, 1/5), -1];

    problem = radialfluxrotor2dfemmprob ( ...
        design.thetap, design.thetam, design.tm, design.tbi, [false, true], ...
        design.Rmo, 'FemmProblem', problem, 'NPolePairs', 1, ...
        'MagArrangement', 'NS', 'PolarisationType', design.MagnetPolarisation, ...
        'MagnetMaterial', matinds(2), 'BackIronMaterial', matinds(3), ...
        'OuterRegionsMaterial', matinds(1), 'MagnetSpaceMaterial', matinds(1), ...
        'MagnetRegionMeshSize', magnetMesh, ...
        'BackIronRegionMeshSize', backIronMesh, ...
        'OuterRegionsMeshSize', outerMesh);

    pitch = 2 * design.thetap;
    problem = add_gap_air_region (problem, design.Rmo, design.Rmo + design.g/3, ...
                                  pitch, seamName, gapName, materialNames{1});
    problem = consolidate_seam (problem, seamName, pitch);
end


function problem = build_stator_tile (design, library, materialNames)
% Draw one 10-degree slot (coil, teeth and yoke) and extend the air gap from
% the second AGE arc to the stator bore.
    problem = newproblem_mfemm ('planar', 'Depth', design.ls, 'MinAngle', 30);
    [problem, matinds] = addmaterials_mfemm (problem, materialNames, ...
                                             'MaterialsLibrary', library);
    [problem, ~, ~] = addboundaryprop_mfemm (problem, 'seam', 4);
    [problem, ~, ~] = addboundaryprop_mfemm (problem, 'gap', 6);
    seamName = problem.BoundaryProps(end-1).Name;
    gapName = problem.BoundaryProps(end).Name;

    Rs = design.Rmo + design.g + design.tc(1) + design.tsb + design.ty/2;
    Inputs = stator_inputs (design, Rs, matinds);

    [problem, statorinfo] = radialfluxstator2dfemmprob ( ...
        design.Qs, design.Poles, Rs, design.thetap, design.thetac, ...
        design.thetasg, design.ty, design.tc(1), design.tsb, design.tsg, [1, 0], ...
        'NSlots', 1, 'SlotMaterial', matinds(5), 'ShoeGapMaterial', matinds(1), ...
        'ShoeGapRegionMeshSize', Inputs.ShoeGapRegionMeshSize, ...
        'FemmProblem', problem, 'NWindingLayers', design.CoilLayers, ...
        'SplitSlot', false, 'Tol', Inputs.Tol, 'DrawCoilInsulation', false, ...
        'CoilBaseFraction', Inputs.CoilBaseFraction);

    [problem, statorinfo] = stator_iron_boundary (problem, design, Inputs, ...
                                                  statorinfo, Rs);
    [problem, statorinfo] = stator_outer_regions (problem, design, Inputs, ...
                                                  statorinfo, Rs, matinds(1));

    % Yoke label.
    [yokeX, yokeY] = pol2cart (2*pi/design.Qs/2, Rs);
    problem = addblocklabel_mfemm (problem, yokeX, yokeY, ...
        'BlockType', materialNames{4}, 'MaxArea', Inputs.YokeRegionMeshSize);

    % Coil labels and circuits for the first slot; per-instance overrides
    % supply the other slots' phase connections.
    problem = addcoil_labels (problem, design, Inputs, matinds(5), ...
                              statorinfo.CoilLabelLocations);

    % Air gap extension from the second AGE arc to the stator bore.
    pitch = 2*pi/design.Qs;
    problem = add_gap_air_region (problem, design.Rai, design.Rmo + 2*design.g/3, ...
                                  pitch, seamName, gapName, materialNames{1});
    problem = consolidate_seam (problem, seamName, pitch);
end


function Inputs = stator_inputs (design, Rs, matinds)
    Inputs.NPolePairs = 1/6;
    Inputs.Tol = 1e-5;
    Inputs.CoilBaseFraction = design.tc(2) / design.tc(1);
    Inputs.YokeRegionMeshSize = mean ([ ...
        choosemesharea_mfemm(design.ty, 2*(design.Rym*design.thetap), 1/10), ...
        choosemesharea_mfemm(design.tc(1), (design.Rcm*(design.thetas - mean(design.thetac))), 1/10)]);
    Inputs.CoilRegionMeshSize = choosemesharea_mfemm (design.tc(1), ...
                                                      (design.Rcm*mean(design.thetac)));
    Inputs.AirGapMeshSize = choosemesharea_mfemm (design.g, ...
                                                  (design.Rmm*design.thetap), 1/10);
    Inputs.ShoeGapRegionMeshSize = choosemesharea_mfemm (max (design.tsg, design.tsb), ...
                                                         (design.Rmo*design.thetasg), 1/20);
    Inputs.StatorOuterRegionSize = [2*design.tm, 10*design.tm];
    Inputs.StatorOuterRegionsMeshSize = [choosemesharea_mfemm(design.tm, Rs*design.thetap, 1/5), -1];
    Inputs.StatorOuterRegionMaterials = {design.MagFEASimMaterials.AirGap, ...
                                         design.MagFEASimMaterials.AirGap};
    Inputs.GapMaterialIndex = matinds(1);
end


function problem = add_gap_air_region (problem, rConnect, rNew, pitchRad, ...
                                       seamName, gapName, airMaterial)
% Add an annular air region between an existing boundary at rConnect and a new
% arc at rNew. The new arc carries the AGE marker; the radial sides are seams.
    maxSegDegrees = 1.0;
    pitchDegrees = rad2deg (pitchRad);
    coords = reshape ([problem.Nodes.Coords], 2, [])';
    bottomNode = find_node_at (coords, rConnect, 0.0);
    topNode = find_node_at (coords, rConnect, pitchDegrees);
    assert (~isempty (bottomNode) && ~isempty (topNode), ...
            'could not find existing gap boundary nodes at r=%g', rConnect);

    [bx, by] = pol2cart (0.0, rNew);
    [tx, ty] = pol2cart (pitchRad, rNew);
    [problem, ~, newIds] = addnodes_mfemm (problem, [bx; tx], [by; ty]);
    newBottom = newIds(1);
    newTop = newIds(2);

    problem = addsegments_mfemm (problem, bottomNode, newBottom, ...
                                 'BoundaryMarker', seamName);
    problem = addsegments_mfemm (problem, topNode, newTop, ...
                                 'BoundaryMarker', seamName);
    problem = addarcsegments_mfemm (problem, newBottom, newTop, pitchDegrees, ...
        'MaxSegDegrees', maxSegDegrees, 'BoundaryMarker', gapName);

    [lx, ly] = pol2cart (pitchRad/2, 0.5*(rConnect + rNew));
    problem = addblocklabel_mfemm (problem, lx, ly, 'BlockType', airMaterial, ...
                                   'MaxArea', -1);
end


function problem = consolidate_seam (problem, seamName, pitchRad)
% Replace the collinear seam segments on each tile edge with one long radial
% segment. Tangle splits it at the retained collinear vertices, so the periodic
% boundary still consists of exactly the two declared sides.
    pitchDegrees = rad2deg (pitchRad);
    coords = reshape ([problem.Nodes.Coords], 2, [])';
    keep = true (1, numel (problem.Segments));
    longSegments = zeros (0, 2);
    for edge = [0, pitchDegrees]
        endpoints = [];
        for k = 1:numel (problem.Segments)
            seg = problem.Segments(k);
            a = coords(seg.n0 + 1, :);
            b = coords(seg.n1 + 1, :);
            d = b - a;
            length = hypot (d(1), d(2));
            if length <= 1e-12
                continue;
            end
            if abs (d(1)*a(2) - d(2)*a(1)) / length > 1e-9
                continue;
            end
            mid = 0.5 * (a + b);
            angle = mod (atan2d (mid(2), mid(1)), 360);
            if angular_distance (angle, edge) < 1e-6
                endpoints = [endpoints; seg.n0; seg.n1]; %#ok<AGROW>
                keep(k) = false;
            end
        end
        endpoints = unique (endpoints);
        assert (numel (endpoints) >= 2, ...
                'seam edge at %.3g degrees has no segments', edge);
        radius = hypot (coords(endpoints + 1, 1), coords(endpoints + 1, 2));
        [~, order] = sort (radius);
        longSegments(end+1, :) = [endpoints(order(1)), endpoints(order(end))]; %#ok<AGROW>
    end
    problem.Segments = problem.Segments(keep);
    for k = 1:size (longSegments, 1)
        problem = addsegments_mfemm (problem, longSegments(k,1), longSegments(k,2), ...
                                     'BoundaryMarker', seamName);
    end
end


function node = find_node_at (coords, radius, angleDegrees)
    node = [];
    best = inf;
    for k = 1:size (coords, 1)
        r = hypot (coords(k,1), coords(k,2));
        if abs (r - radius) > 1e-7
            continue;
        end
        angle = mod (atan2d (coords(k,2), coords(k,1)), 360);
        if angular_distance (angle, angleDegrees) > 1e-6
            continue;
        end
        if abs (r - radius) < best
            best = abs (r - radius);
            node = k - 1;
        end
    end
end


function d = angular_distance (a, b)
    d = abs (mod (a - b + 180, 360) - 180);
end


function problem = addcoil_labels (problem, design, Inputs, coilMatInd, locations)
    for i = 1:design.Phases
        cname = num2str (i);
        if ~hascircuit_mfemm (problem, cname)
            problem = addcircuit_mfemm (problem, cname);
        end
        problem = setcircuitcurrent (problem, cname, 0);
    end
    phases = design.WindingLayout.Phases(1, :);
    for layer = 1:design.CoilLayers
        problem = addblocklabel_mfemm (problem, locations(layer,1), locations(layer,2), ...
            'BlockType', design.MagFEASimMaterials.ArmatureCoil, ...
            'InCircuit', num2str (abs (phases(layer))), ...
            'Turns', design.CoilTurns * sign (phases(layer)), ...
            'MaxArea', Inputs.CoilRegionMeshSize);
    end
end


function [FemmProblem, statorinfo] = stator_iron_boundary (FemmProblem, design, Inputs, statorinfo, Rs)
% Copy of the external-armature branch of slottedfemmprob_radial's local
% helper, with the fractional pole-pair count for a single slot tile.
    elcount = elementcount_mfemm (FemmProblem);
    statorirongp = getgroupnumber_mfemm (FemmProblem, 'StatorIronOutline');
    [edgenodes(:,1), edgenodes(:,2)] = pol2cart ( ...
        [0; design.thetap*2*Inputs.NPolePairs; 0; ...
         design.thetap*Inputs.NPolePairs; design.thetap*2*Inputs.NPolePairs], ...
        [design.Rmo+design.g; design.Rmo+design.g; Rs+design.ty/2; ...
         Rs+design.ty/2; Rs+design.ty/2]);
    [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, edgenodes(:,1), edgenodes(:,2));
    statorinfo.node_id_set_stator_inner = [nodeids(1), nodeids(2)];
    statorinfo.node_id_set_stator_outer = [nodeids(3), nodeids(4), nodeids(5)];
    statorinfo.node_id_pair_stator_bottom_edge = [nodeids(1), nodeids(3)];
    statorinfo.node_id_pair_stator_top_edge = [nodeids(2), nodeids(5)];
    [FemmProblem, ~, statorboundname] = addboundaryprop_mfemm (FemmProblem, ...
        'Radial Stator Back Iron Periodic', 4);
    segprop = struct ('BoundaryMarker', statorboundname, 'InGroup', statorirongp);
    [FemmProblem, statorinfo.BottomSegInds] = addsegments_mfemm (FemmProblem, ...
        statorinfo.node_id_pair_stator_bottom_edge(1), ...
        statorinfo.node_id_pair_stator_bottom_edge(2), segprop);
    [FemmProblem, statorinfo.TopSegInds] = addsegments_mfemm (FemmProblem, ...
        statorinfo.node_id_pair_stator_top_edge(1), ...
        statorinfo.node_id_pair_stator_top_edge(2), segprop);
    FemmProblem = addarcsegments_mfemm (FemmProblem, statorinfo.node_id_set_stator_inner(1), ...
        statorinfo.OuterNodes(1), rad2deg(((2*pi/design.Qs)-design.thetac(1))/2), ...
        'InGroup', statorirongp);
    FemmProblem = addarcsegments_mfemm (FemmProblem, statorinfo.OuterNodes(4), ...
        statorinfo.node_id_set_stator_inner(2), rad2deg(((2*pi/design.Qs)-design.thetac(1))/2), ...
        'InGroup', statorirongp);
    FemmProblem = addarcsegments_mfemm (FemmProblem, statorinfo.node_id_set_stator_outer(1), ...
        statorinfo.node_id_set_stator_outer(2), rad2deg(design.thetap*Inputs.NPolePairs));
    FemmProblem.ArcSegments(end).InGroup = statorirongp;
    FemmProblem = addarcsegments_mfemm (FemmProblem, statorinfo.node_id_set_stator_outer(2), ...
        statorinfo.node_id_set_stator_outer(3), rad2deg(design.thetap*Inputs.NPolePairs));
    FemmProblem.ArcSegments(end).InGroup = statorirongp;
    FemmProblem = translatenewelements_mfemm (FemmProblem, elcount, 0, 0);
end


function [FemmProblem, statorinfo] = stator_outer_regions (FemmProblem, design, Inputs, statorinfo, Rs, GapMatInd)
% Copy of the external-armature branch of slottedfemmprob_radial's local
% helper, with the fractional pole-pair count for a single slot tile.
    elcount = elementcount_mfemm (FemmProblem);
    statorinfo.touterregion = Inputs.StatorOuterRegionSize;
    statorinfo.routerregion = cumsum ([Rs + design.ty/2, statorinfo.touterregion]);
    sectionangle = design.thetap * 2 * Inputs.NPolePairs;
    statorinfo.node_id_sets_stator_outer_region = statorinfo.node_id_set_stator_outer;
    for ind = 1:numel (statorinfo.touterregion)
        [edgenodes(:,1), edgenodes(:,2)] = pol2cart ( ...
            [0; design.thetap*Inputs.NPolePairs; design.thetap*2*Inputs.NPolePairs], ...
            repmat (statorinfo.routerregion(ind+1), [3,1]));
        [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, edgenodes(:,1), edgenodes(:,2));
        statorinfo.node_id_sets_stator_outer_region = [ ...
            statorinfo.node_id_sets_stator_outer_region; nodeids];
        [FemmProblem, ~, boundname] = addboundaryprop_mfemm (FemmProblem, ...
            'Radial Stator Outer Periodic', 4);
        segprops = struct ('BoundaryMarker', boundname, 'InGroup', 0);
        [FemmProblem, statorinfo.BottomSegInds(ind)] = addsegments_mfemm (FemmProblem, ...
            statorinfo.node_id_sets_stator_outer_region(ind,1), ...
            statorinfo.node_id_sets_stator_outer_region(ind+1,1), segprops);
        [FemmProblem, statorinfo.TopSegInds(ind)] = addsegments_mfemm (FemmProblem, ...
            statorinfo.node_id_sets_stator_outer_region(ind,3), ...
            statorinfo.node_id_sets_stator_outer_region(ind+1,3), segprops);
        FemmProblem = addarcsegments_mfemm (FemmProblem, ...
            [statorinfo.node_id_sets_stator_outer_region(ind+1,1), ...
             statorinfo.node_id_sets_stator_outer_region(ind+1,2)], ...
            [statorinfo.node_id_sets_stator_outer_region(ind+1,2), ...
             statorinfo.node_id_sets_stator_outer_region(ind+1,3)], ...
            rad2deg(repmat(design.thetap*Inputs.NPolePairs, 1, 2)));
        [labelloc(1), labelloc(2)] = pol2cart (design.thetap * Inputs.NPolePairs, ...
                                               mean (statorinfo.routerregion(ind:ind+1)));
        FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (labelloc(1), labelloc(2), ...
            'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
            'MaxArea', Inputs.StatorOuterRegionsMeshSize(ind));
    end
    FemmProblem = translatenewelements_mfemm (FemmProblem, elcount, 0, 0);
end


function overrides = stator_phase_overrides (design)
    phases = design.WindingLayout.Phases;
    assert (size (phases, 1) == design.Qs, ...
            'unexpected winding layout size for the tiled overrides');
    overrides = {};
    for layer = 1:design.CoilLayers
        circuit = cell (1, design.Qs);
        turnScale = zeros (1, design.Qs);
        baseSign = sign (phases(1, layer));
        for slot = 1:design.Qs
            circuit{slot} = num2str (abs (phases(slot, layer)));
            turnScale(slot) = sign (phases(slot, layer)) * baseSign;
        end
        overrides{end+1} = struct ('tile', 'stator', ...
                                   'label', sprintf ('coil%d', layer), ...
                                   'circuit', {circuit}, 'turnScale', turnScale);
    end
end


function tile = tile_from_problem (problem, name, count, seamName, gapName)
    tile = struct ();
    tile.name = name;
    geometry = struct ();
    geometry.nodes = nodes_json (problem);
    geometry.segments = segments_json (problem);
    geometry.arcs = arcs_json (problem);
    geometry.labels = labels_json (problem);
    tile.geometry = geometry;
    tile.repeat = struct ('kind', 'rotation', 'center', [0, 0], ...
                          'count', count, 'closure', 'closed');
    tile.seamBoundary = seamName;
    tile.airGapBoundary = gapName;
end


function out = nodes_json (problem)
    out = {};
    for k = 1:numel (problem.Nodes)
        node = problem.Nodes(k);
        out{end + 1} = struct ('x', node.Coords(1), 'y', node.Coords(2), ...
                               'group', node.InGroup); %#ok<AGROW>
    end
end


function out = segments_json (problem)
    out = {};
    for k = 1:numel (problem.Segments)
        segment = problem.Segments(k);
        entry = struct ('n0', segment.n0, 'n1', segment.n1, ...
                        'maxSideLength', segment.MaxSideLength, ...
                        'hidden', double (segment.Hidden), 'group', segment.InGroup);
        if isfield (segment, 'BoundaryMarker') && ischar (segment.BoundaryMarker) ...
                && ~isempty (segment.BoundaryMarker)
            entry.boundary = segment.BoundaryMarker;
        end
        out{end + 1} = entry; %#ok<AGROW>
    end
end


function out = arcs_json (problem)
    out = {};
    for k = 1:numel (problem.ArcSegments)
        arc = problem.ArcSegments(k);
        entry = struct ('n0', arc.n0, 'n1', arc.n1, 'arcLength', arc.ArcLength, ...
                        'maxSegDegrees', arc.MaxSegDegrees, ...
                        'hidden', double (arc.Hidden), 'group', arc.InGroup);
        if isfield (arc, 'BoundaryMarker') && ischar (arc.BoundaryMarker) ...
                && ~isempty (arc.BoundaryMarker)
            entry.boundary = arc.BoundaryMarker;
        end
        out{end + 1} = entry; %#ok<AGROW>
    end
end


function out = labels_json (problem)
    out = {};
    coilCount = 0;
    for k = 1:numel (problem.BlockLabels)
        label = problem.BlockLabels(k);
        if isfield (label, 'InCircuit') && ischar (label.InCircuit) ...
                && ~isempty (label.InCircuit)
            coilCount = coilCount + 1;
            labelName = sprintf ('coil%d', coilCount);
        else
            labelName = sprintf ('label%d', k);
        end
        entry = struct ('name', labelName, ...
                        'x', label.Coords(1), 'y', label.Coords(2), ...
                        'turns', label.Turns, 'maxArea', label.MaxArea, ...
                        'group', label.InGroup);
        if isfield (label, 'BlockType') && ischar (label.BlockType) ...
                && ~isempty (label.BlockType)
            entry.material = label.BlockType;
        end
        if isfield (label, 'InCircuit') && ischar (label.InCircuit) ...
                && ~isempty (label.InCircuit)
            entry.circuit = label.InCircuit;
        end
        if isfield (label, 'MagDir') && ~isempty (label.MagDir) && ~isnan (label.MagDir)
            entry.magDir = label.MagDir;
        end
        if isfield (label, 'MagDirFctn') && ischar (label.MagDirFctn) ...
                && ~isempty (label.MagDirFctn)
            entry.magDirFctn = label.MagDirFctn;
        end
        out{end + 1} = entry; %#ok<AGROW>
    end
end


function [merged, mapping] = merge_materials (problems)
    merged = {};
    mapping = cell (1, numel (problems));
    for p = 1:numel (problems)
        mapping{p} = containers.Map ();
        for k = 1:numel (problems{p}.Materials)
            prop = problems{p}.Materials(k);
            existing = find (cellfun (@(m) strcmp (m.name, prop.Name), merged));
            if isempty (existing)
                merged{end + 1} = material_json (prop); %#ok<AGROW>
                existing = numel (merged);
            end
            mapping{p}(prop.Name) = existing;
        end
    end
end


function entry = material_json (prop)
    entry = struct ('name', prop.Name, 'muX', prop.Mu_x, 'muY', prop.Mu_y, ...
                    'Hc', prop.H_c, 'sigma', prop.Sigma, 'lamD', prop.d_lam, ...
                    'phiH', prop.Phi_h, 'phiHx', prop.Phi_hx, 'phiHy', prop.Phi_hy, ...
                    'lamType', prop.LamType, 'lamFill', prop.LamFill, ...
                    'strands', prop.NStrands, 'wireD', prop.WireD);
    if prop.J_re ~= 0 || prop.J_im ~= 0
        entry.J = [prop.J_re, prop.J_im];
    end
    if ~isempty (prop.BHPoints)
        bh = cell (1, size (prop.BHPoints, 1));
        for r = 1:size (prop.BHPoints, 1)
            bh{r} = prop.BHPoints(r, :);
        end
        entry.bh = bh;
    end
end


function [merged, mapping] = merge_boundaries (problems)
    merged = {};
    mapping = cell (1, numel (problems));
    for p = 1:numel (problems)
        mapping{p} = containers.Map ();
        for k = 1:numel (problems{p}.BoundaryProps)
            prop = problems{p}.BoundaryProps(k);
            existing = find (cellfun (@(b) strcmp (b.name, prop.Name), merged));
            if isempty (existing)
                merged{end + 1} = boundary_json (prop); %#ok<AGROW>
                existing = numel (merged);
            end
            mapping{p}(prop.Name) = existing;
        end
    end
end


function entry = boundary_json (prop)
    entry = struct ('name', prop.Name, 'type', bdrytype (prop.BdryType), ...
                    'format', prop.BdryType, 'A0', prop.A0, 'A1', prop.A1, ...
                    'A2', prop.A2, 'phi', prop.Phi, 'mu', prop.Mu_ssd, ...
                    'sigma', prop.Sigma_ssd, ...
                    'c0', [prop.c0, prop.c0i], 'c1', [prop.c1, prop.c1i], ...
                    'innerAngle', prop.InnerAngle, 'outerAngle', prop.OuterAngle);
end


function [merged, mapping] = merge_circuits (problems)
    merged = {};
    mapping = cell (1, numel (problems));
    for p = 1:numel (problems)
        mapping{p} = containers.Map ();
        for k = 1:numel (problems{p}.Circuits)
            prop = problems{p}.Circuits(k);
            existing = find (cellfun (@(c) strcmp (c.name, prop.Name), merged));
            if isempty (existing)
                merged{end + 1} = struct ('name', prop.Name, 'type', prop.CircType, ...
                                          'amps', [prop.TotalAmps_re, prop.TotalAmps_im]); %#ok<AGROW>
                existing = numel (merged);
            end
            mapping{p}(prop.Name) = existing;
        end
    end
end


function name = bdrytype (format)
    switch format
        case 0, name = 'dirichlet';
        case 1, name = 'smallskin';
        case 2, name = 'mixed';
        case 4, name = 'periodic';
        case 5, name = 'antiperiodic';
        case 6, name = 'age';
        case 7, name = 'age-antiperiodic';
        otherwise, name = 'dirichlet';
    end
end


function design = example_design ()
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
    design.thetam = (2 * pi / design.Poles) * 0.667;
    design.thetacg = 84.7661e-3;
    design.thetacy = 93.0181e-3;
    design.thetasg = 38.6622e-3;
    design.ls = 88.9e-3;
    design.ArmatureType = 'external';
    design.MagnetPolarisation = 'radial';
    design.MagFEASimMaterials.AirGap = 'Air';
    design.MagFEASimMaterials.Magnet = 'NdFeB 40 MGOe';
    design.MagFEASimMaterials.FieldBackIron = '1117 Steel';
    design.MagFEASimMaterials.ArmatureYoke = '1117 Steel';
    design.MagFEASimMaterials.ArmatureCoil = '36 AWG';
end
