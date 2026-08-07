function report = compare_radial_machine_static_results (reference, candidate, varargin)
%COMPARE_RADIAL_MACHINE_STATIC_RESULTS Compare static machine sweeps.
%
% report = compare_radial_machine_static_results (reference, candidate)
%
% A check passes when
%
%   max(abs(candidate-reference)) <= AbsoluteTolerance
%                                    + RelativeTolerance * signal scale
%
% The defaults allow the small discretisation changes expected when the
% rotor is redrawn and remeshed, while still detecting waveform regressions.

    options.RelativeTolerance = 0.03;
    options.FluxAbsoluteTolerance = 1e-7;
    options.FluxDensityAbsoluteTolerance = 1e-4;
    options = parse_options (options, varargin{:});

    assert (reference.schemaVersion == candidate.schemaVersion, ...
            'XFEMM:StaticMachineComparison', 'Result schema versions differ.');

    report.positions = exact_check (reference.positions, candidate.positions, ...
                                    100*eps, 'positions');
    report.fluxLinkage = tolerance_check (reference.fluxLinkage, ...
        candidate.fluxLinkage, options.RelativeTolerance, ...
        options.FluxAbsoluteTolerance, 'coil flux linkage');
    if isfield (reference, 'coilFluxDensity') ...
            && isfield (candidate, 'coilFluxDensity')
        report.coilFluxDensity = tolerance_check (...
            reference.coilFluxDensity, candidate.coilFluxDensity, ...
            options.RelativeTolerance, options.FluxDensityAbsoluteTolerance, ...
            'coil-region flux density');
    end
    if isfield (reference, 'circuitFluxLinkage') ...
            && isfield (candidate, 'circuitFluxLinkage')
        report.circuitFluxLinkage = tolerance_check (...
            reference.circuitFluxLinkage, candidate.circuitFluxLinkage, ...
            options.RelativeTolerance, options.FluxAbsoluteTolerance, ...
            'FEMM circuit flux linkage');
    end

    checks = fieldnames (report);
    report.passed = true;
    for ind = 1:numel (checks)
        if isstruct (report.(checks{ind}))
            report.passed = report.passed && report.(checks{ind}).passed;
        end
    end

    fprintf ('Static radial-machine comparison:\n');
    for ind = 1:numel (checks)
        check = report.(checks{ind});
        if isstruct (check)
            fprintf ('  %-34s max |delta| = %.6g, limit = %.6g, %s\n', ...
                     check.name, check.maxAbsoluteError, check.limit, ...
                     pass_text (check.passed));
        end
    end

    assert (report.passed, 'XFEMM:StaticMachineComparison', ...
            'Static radial-machine results exceeded one or more tolerances.');
end


function check = exact_check (reference, candidate, absoluteTolerance, name)
    assert_same_size (reference, candidate, name);
    check.name = name;
    check.maxAbsoluteError = max_or_zero (abs (candidate(:) - reference(:)));
    check.scale = max_or_zero (abs (reference(:)));
    check.limit = absoluteTolerance;
    check.passed = check.maxAbsoluteError <= check.limit;
end


function check = tolerance_check (reference, candidate, relativeTolerance, ...
                                   absoluteTolerance, name)
    assert_same_size (reference, candidate, name);
    check.name = name;
    check.maxAbsoluteError = max_or_zero (abs (candidate(:) - reference(:)));
    check.scale = max ([max_or_zero(abs(reference(:))), ...
                        max_or_zero(abs(candidate(:)))]);
    check.limit = absoluteTolerance + relativeTolerance * check.scale;
    check.relativeError = check.maxAbsoluteError / max (check.scale, eps);
    check.passed = check.maxAbsoluteError <= check.limit;
end


function assert_same_size (reference, candidate, name)
    assert (isequal (size (reference), size (candidate)), ...
            'XFEMM:StaticMachineComparison', ...
            '%s arrays have different sizes.', name);
end


function value = max_or_zero (values)
    if isempty (values)
        value = 0;
    else
        value = max (values);
    end
end


function value = pass_text (passed)
    if passed
        value = 'PASS';
    else
        value = 'FAIL';
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
