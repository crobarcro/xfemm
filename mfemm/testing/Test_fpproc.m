%% Fpproc Unit tests

fprintf ('Running Test_fpproc\n');

% solve reference problem
femprob = fullfile (mfemmdeps.getmfilepath ('Test_fpproc.m'), 'reference.fem');
[ansfilename, femfilename] = analyse_mfemm(femprob);

solution = fpproc (ansfilename);

%% line integral

% make a contour
solution.clearcontour ()
solution.addcontour (-2, 2);
solution.addcontour (2, 2)
solution.addcontour (2, -2)
solution.addcontour (-2, -2)
solution.addcontour (-2, 2);

% perform some integrals
Bn = solution.lineintegral (0)
Ht = solution.lineintegral (1)
contour_length = solution.lineintegral (2)
force_from_stress_tensor = solution.lineintegral (3)
torque_from_stress_tensor = solution.lineintegral (4)
Bn2 = solution.lineintegral (5)

% Explicit deletion followed by variable cleanup must not try to destroy the
% native object a second time under Octave.
lastwarn('');
delete(solution);
clear solution;
[deleteWarning, ~] = lastwarn();
assert(isempty(deleteWarning), ...
       'MFEMM:fpproc:doubleDelete', ...
       'Explicit fpproc deletion caused a second native handle deletion.');

