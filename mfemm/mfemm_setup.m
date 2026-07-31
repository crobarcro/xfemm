function varargout = mfemm_setup(varargin)
% mfemm_setup: set up the matlab interface to xfemm
%
% Syntax
%
% mfemm_setup ()
% mfemm_setup ('Parameter', Value)
%
% Description
%
% mfemm_setup performs various operations to set up mfemm. At it's
% simplest, it will just add some directories to your Matlab/Octave path.
% At it's most complicated it will compile from scratch a set of libraries
% and matlab mex files. What it does, or attempts to do, will depend on
% what has already been done before, or supplied with the package, or what
% you tell it to do using various options.
%
% mfemm_setup has several optional arguments that may be supplied as
% Parameter-Value pairs. The possible options are shown below. Some may
% require further reading of the help below to fully understand.
%
% 'RunTests' - true or false flag. If true, a set of test runs of the mex
%   files will be performed. Defaults to false.
%
% 'ForceMexRecompile' - true or false flag. If true, forces the mex files
%   all be recompiled. Library files will not be recompiled unless they are
%   not present. Defaults to false.
%
% 'Debug' - true or false flag. Mex files will be built with debugging
%   symbols. Defaults to false if not supplied.
%
% 'Verbose' - true or false flag. Verbose build output if true. Defaults to
%   false if not supplied.
%
% 'Clean' - true/false flag indicating whether to run make clean before
%
% 'SessionOnly' - true/false flag. If true, compile the high-level analysis
%   session MEX gateway and its magnetic post-processing gateway dependency.
%   This is useful for focused development and CI testing without rebuilding
%   unrelated solver gateways. Defaults false.
%
%
% These may be supplied in any order. An example call the mfemm_setup might
% be:
%
%   mfemm_setup ('DoDebug', false, 'ForceMexRecompile', true)
%
% A large proportion of mfemm's functionality depends on mex interfaces to
% C++ code. This code must be compiled into mex functions to be accessible
% from the matlab or octave command line. In Matlab therefore, you must
% have previously run:
%
% mex -setup
%
% And selected a valid C++ compiler before attempting to run this
% mfemm_setup function.
%
%

    if nargin < 1

        if nargout > 0
            % return the subfunctions for testing/use
            varargout {1} = { @makelibs };
            return;
        end
    end
    
    % restore user's directory once we're done
    userdir = pwd ();
    CC = onCleanup (@() cd (userdir));

    % add the required directories to the path
    thisfilepath = which('mfemm_setup.m');

    if isoctave
        thisfilepath = canonicalize_file_name(thisfilepath);
    end

    thisfilepath = fileparts (thisfilepath);

    % add mfile paths
    addpath(thisfilepath);
    addpath (fullfile (thisfilepath, 'preproc'));
    addpath (fullfile (thisfilepath, 'postproc'));
    addpath (fullfile (thisfilepath, 'examples'));
    addpath (fullfile (thisfilepath, 'testing'));
    addpath (fullfile (thisfilepath, 'visualisation'));

    % parse inputs
    options.ForceMexRecompile = false;
    options.Debug = false;
    options.Verbose = false;
    options.RunTests = false;
    options.DoCrossBuildWin64 = false;
    options.W64CrossBuildMexLibsDir = '';
    options.Clean = true;
    options.Tidy = true;
    options.DebugSymbols = false;
    options.ExtraMEXFLAGS = '';
    options.SessionOnly = false;

    options = mfemmdeps.parseoptions (options, varargin);

    % make architecture specific directory for mex files if it doesn't
    % already exist
    warning off MATLAB:MKDIR:DirectoryExists
    if options.DoCrossBuildWin64
        mexdir = fullfile(thisfilepath, ['xfemm_mex_files_for_', 'win64']);
        thismexext = 'mexw64';
    else
        mexdir = fullfile(thisfilepath, ['xfemm_mex_files_for_' computer('arch')]);
        thismexext = mexext ();
    end
    mkdir (mexdir);
    warning on MATLAB:MKDIR:DirectoryExists

    % add it to the path
    addpath (mexdir);

    if options.SessionOnly
        needcompile = ~(exist('session_interface_mex', 'file') == 3) ...
                      || ~(exist('fpproc_interface_mex', 'file') == 3) ...
                      || options.ForceMexRecompile;
    else
        needcompile = ~(exist('mexfsolver', 'file') == 3) ...
                      || ~(exist('mexfmesher', 'file') == 3) ...
                      || ~(exist('fpproc_interface_mex', 'file') == 3) ...
                      || ~(exist('session_interface_mex', 'file') == 3) ...
                      || options.ForceMexRecompile;
    end

    if needcompile

        if options.Verbose
            fprintf('Compiling mex functions for mfemm.\n');
        end

        CC = onCleanup(@() cd(pwd));
        cd (thisfilepath);

        if ~exist (fullfile ('..', 'cfemm'), 'dir')
            error ( 'MFEMM:Build', ...
                    [ 'The cfemm directory was not found in the expected location,\n', ...
                      'you must preserve xfemm directory stucture for mfemm to build.\n', ...
                      'Compilation terminating' ] ...
                  );
        end

        % Unlike the CMake build, the MEX makefiles compile libfemm sources
        % directly.  Create CMake's generated version header when building
        % from a clean source checkout.
        ensurefemmversionheader (fullfile (thisfilepath, '..', 'cfemm'));

        makefilenames = {'MMakefile_fmesher.m', ...
                         'MMakefile_fsolver.m', ...
                         'MMakefile_session.m', ...
                         'MMakefile_fpproc.m', ...
                         'MMakefile_hpproc.m', ...
                         'MMakefile_hsolver.m' };

        if options.SessionOnly
            makefilenames = {'MMakefile_session.m', 'MMakefile_fpproc.m'};
        end

        if options.ForceMexRecompile
            % run make clean for all projects to force complete
            % recompilation
            if options.SessionOnly
                sessionmexfiles = {'session_interface_mex', 'fpproc_interface_mex'};
                for ind = 1:numel(sessionmexfiles)
                    sessionmex = fullfile (mexdir, [sessionmexfiles{ind}, '.', thismexext]);
                    if exist(sessionmex, 'file')
                        delete(sessionmex);
                    end
                end
            else
                delete (fullfile (mexdir, ['*.', thismexext]));
            end
            if options.Clean
                for ind = 1:numel(makefilenames)
                    mmake.make ('clean', makefilenames{ind});
                end
            end
        end

        % now invoke mmake for all files
        ws = warning( 'off', 'MATLAB:mex:GccVersion');
        for ind = 1:numel(makefilenames)
            if options.DoCrossBuildWin64
                mmake.make ( '', makefilenames{ind}, ...
                             'DoCrossBuildWin64', true, ...
                             'FcnMakeFileArgs', { 'DoCrossBuildWin64', true, ...
                                                  'W64CrossBuildMexLibsDir', options.W64CrossBuildMexLibsDir, ...
                                                  'Verbose', true, ...
                                                  'Debug', options.Debug, ...
                                                  'DebugSymbols', options.DebugSymbols, ...
                                                  'ExtraMEXFLAGS', options.ExtraMEXFLAGS } ...
                           );
            else
                mmake.make ( '', makefilenames{ind}, ...
                             'FcnMakeFileArgs',  { 'Debug', options.Debug, ...
                                                   'DebugSymbols', options.DebugSymbols, ...
                                                   'Verbose', options.Verbose, ...
                                                   'ExtraMEXFLAGS', options.ExtraMEXFLAGS })
            end
        end
        warning (ws);

        if options.Tidy
            for ind = 1:numel(makefilenames)
                mmake.make ('tidy', makefilenames{ind})
            end
        end

        % copy over the created mex files to the mex directory
        mexfilenames = {'mexfmesher', ...
                        'mexfsolver', ...
                        'session_interface_mex', ...
                        'fpproc_interface_mex', ...
                        'mexhsolver', ...
                        'hpproc_interface_mex'};

        if options.SessionOnly
            mexfilenames = {'session_interface_mex', 'fpproc_interface_mex'};
        end

        for ind = 1:numel(mexfilenames)
            movefile ([mexfilenames{ind}, '.',  thismexext], mexdir);
        end

         % force a path refresh
         rehash ()

    else

        if options.Verbose
            fprintf('Not compiling mex functions as they already exist.\n');
        end

    end

    if options.RunTests

        fmesher_test_file = fullfile (thisfilepath, '..', 'cfemm', 'fmesher', 'test', 'Temp.fem');

        fprintf ('Running mexfmesher on file:\n%s\n', fmesher_test_file);
        mexfmesher (fmesher_test_file);

        fsolver_test_file = fullfile (thisfilepath, '..', 'cfemm', 'fmesher', 'test', 'Temp');
        fprintf ('Running mexfsolver on file:\n%s\n', fsolver_test_file);
        mexfsolver (fsolver_test_file, double(true), double(false));

        fprintf ('Running fmesher on file:\n%s\n', fmesher_test_file);
        fmesher (fmesher_test_file);

        fsolver_test_file = fullfile (thisfilepath, '..', 'cfemm', 'fmesher', 'test', 'Temp.fem');
        fprintf ('Running fsolver on file:\n%s\n', fsolver_test_file);
        fsolver (fmesher_test_file);

        fpproc_test_file = fullfile (thisfilepath, '..', 'cfemm', 'fmesher', 'test', 'Temp.ans');
        fprintf ('Loading solution using fpproc from file:\n%s\n', fpproc_test_file);
        solution = fpproc (fpproc_test_file);
        vals = solution.getpointvalues(0.038, 0.207)

        % run test scripts
        run (fullfile (thisfilepath, 'testing', 'Test_fpproc.m'));

    end

    % print message about adding files to the path (octave currently
    % doesn't support displaying help in subfunctions)
    if ~isoctave

        fprintf(1, '\n\n');
        help mfemm_setup>addedtopathmsg

        % display some output
        fprintf (1, '\n');
        help mfemm_setup>printxfemmusermessage

    end

end

function ensurefemmversionheader (cfemmdir)

    headerfile = fullfile (cfemmdir, 'libfemm', 'femmversion.h');
    if exist (headerfile, 'file')
        return;
    end

    version = '0.0.0-dev';
    versionfile = fullfile (cfemmdir, 'VERSION');
    if exist (versionfile, 'file')
        version = strtrim (fileread (versionfile));
    end

    parts = regexp (version, '^(\d+)\.(\d+)(?:\.(\d+))?', 'tokens', 'once');
    if isempty (parts)
        error ('MFEMM:Build', 'Invalid xfemm version string: %s', version);
    end
    if numel (parts) < 3 || isempty (parts{3})
        parts{3} = '0';
    end

    template = fileread (fullfile (cfemmdir, 'libfemm', 'femmversion.h.in'));
    contents = strrep (template, '@XFEMM_VERSION_MAJOR@', parts{1});
    contents = strrep (contents, '@XFEMM_VERSION_MINOR@', parts{2});
    contents = strrep (contents, '@XFEMM_VERSION_PATCH@', parts{3});
    contents = strrep (contents, '@XFEMM_VERSION_STRING@', version);

    fid = fopen (headerfile, 'w');
    if fid < 0
        error ('MFEMM:Build', 'Could not create %s', headerfile);
    end
    cleanup = onCleanup (@() fclose (fid));
    fwrite (fid, contents);

end

function addedtopathmsg ()
% mfemm_setup has added the directories containing mfiles for mfemm to the
% search path. This change will not be saved when you close Matlab/Octave.
% You will need to have these directories on the path in future for mfemm
% to work when you restart Matlab/Octave. You can add them again by using
% the 'Set path' dialog in Matlab, or using the 'addpath' command, or by
% rerunning mfemm_setup, which calls addpath internally to do this.
% mfemm_setup will not recompile the mex functions if they already exist
% and will just modify your path.
%
% If you open the 'Set path' dialog now, you could also save the changes
% that have just been made.
%
% Run the command 'help path' for more information on the path. In general
% you can save changes to the Matlab path that persist across sessions
% using the 'Set path' dialog, or you can make a startup.m file which runs
% every time Matlab starts (.octaverc in Octave) containing the addpath
% commands.
%

end

function printxfemmusermessage ()
% --------    NOTE TO USERS  --------
%
% If you use xfemm, particularly for industrial work, but also academic, it
% will be greatly appreciated if you could write an email stating this and
% how it has supported your work. This is a low-cost way to ensure further
% development and maintenance will continue! Contact the authors on the
% discussion forum, or you will find an email address in the source files.
%
% If you wish to cite xfemm in your work, please use the following:
%
% Crozier, R, Mueller, M., "A New MATLAB and Octave Interface to a
% Popular Magnetics Finite Element Code", Proceedings of the 22nd
% International Conference on Electric Machines (ICEM 2016), September
% 2016.
%
% We would also suggest you cite the original FEMM program.
%
% ----------------------------------

end

function t = isoctave()
% ISOCTAVE.M
% ISOCTAVE  True if the operating environment is octave.
%    Usage: t=mfemmdeps.isoctave();
%
%    Returns 1 if the operating environment is octave, otherwise
%    0 (Matlab)
%
% ---------------------------------------------------------------
%
% COPYRIGHT : (c) NUHAG, Dept.Math., University of Vienna, AUSTRIA
%             http://nuhag.eu/
%             Permission is granted to modify and re-distribute this
%             code in any manner as long as this notice is preserved.
%             All standard disclaimers apply.

    if exist('OCTAVE_VERSION')
        % Only Octave has this variable.
        t=1;
    else
        t=0;
    end

end
