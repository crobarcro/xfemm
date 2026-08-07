function report = compare_radial_machine_result_files (referenceFile, candidateFile, varargin)
%COMPARE_RADIAL_MACHINE_RESULT_FILES Load and compare two saved sweeps.
    referenceData = load (referenceFile, 'results');
    candidateData = load (candidateFile, 'results');
    assert (isfield (referenceData, 'results'), ...
            'Reference MAT file does not contain results.');
    assert (isfield (candidateData, 'results'), ...
            'Candidate MAT file does not contain results.');
    report = compare_radial_machine_static_results (...
        referenceData.results, candidateData.results, varargin{:});
end
