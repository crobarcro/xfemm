function [rules,vars] = MMakefile_fmesher (varargin)

    options.DoCrossBuildWin64 = false;
    options.W64CrossBuildMexLibsDir = '';
    options.Verbose = false;
    options.Debug = false;
    options.DebugSymbols = false;
    options.ExtraMEXFLAGS = '';

    options = mmake.parse_pv_pairs (options, varargin);

    if options.Debug
        options.DebugSymbols = true;
    end

    if isunix || mfemmdeps.isoctave
        ismscompiler = false;
    else
        cc = mex.getCompilerConfigurations ('C');
        if strncmpi (cc.Manufacturer, 'Microsoft', 9)
            ismscompiler = true;
        else
            ismscompiler = false;
        end
    end

    thisfilepath = mfemmdeps.getmfilepath (mfilename);

    thisfilepath = strrep (thisfilepath, '\', '/');

    if ispc || options.DoCrossBuildWin64
        trilibraryflag = '-DCPU86';
    else
        trilibraryflag = '-DLINUX';
    end

    %vars.LDFLAGS = sprintf('${LDFLAGS} -l%s ''-Wl,--no-undefined''', fullfile (matlabroot, 'sys', 'os', 'glnxa64', 'libstdc++.so.6'));
    %vars.LDFLAGS = '${LDFLAGS} -lstdc++ ''-Wl,--no-undefined''';
    vars.LDFLAGS = '${LDFLAGS} -static-libstdc++ ''-Wl,--no-undefined''';

    % flags that will be passed direct to mex
    vars.MEXFLAGS = ['${MEXFLAGS} -D_GLIBCXX_USE_CXX11_ABI=1 -I"../cfemm/fmesher" -I"../cfemm/fmesher/triangle" -I"../cfemm/libfemm" -I"../cfemm/libfemm/liblua" ', trilibraryflag, options.ExtraMEXFLAGS];
    %vars.MEXFLAGS = ['${MEXFLAGS} -I"../cfemm/fmesher" -I"../cfemm/fmesher/triangle" -I"../cfemm/libfemm" -I"../cfemm/libfemm/liblua" ', trilibraryflag];

    if options.Verbose
        vars.MEXFLAGS = [vars.MEXFLAGS, ' -v'];
    end
    
    if options.Debug
        vars.MEXFLAGS = [vars.MEXFLAGS, ' -DDEBUG '];
    else
        vars.MEXFLAGS = [vars.MEXFLAGS, ' -DNDEBUG '];
    end

    if isunix && ~mfemmdeps.isoctave ()
        if options.Debug
            vars.OPTIMFLAGS = '-OO';
            vars.MEXFLAGS = [vars.MEXFLAGS, ' CXXOPTIMFLAGS="-O0"'];
        else
            vars.OPTIMFLAGS = '-O2';
            vars.MEXFLAGS = [vars.MEXFLAGS, ' CXXOPTIMFLAGS="-O2 "'];
        end
    end

    vars.CXXFLAGS = '${CXXFLAGS} -std=c++17 ';

    if options.DebugSymbols
        vars.MEXFLAGS = [vars.MEXFLAGS, ' -g'];
    end

    if mfemmdeps.isoctave ()
        setenv('CXXFLAGS','-std=c++17');
    end

    if ~ismscompiler
        vars.CXXFLAGS = [vars.CXXFLAGS, ' -fpermissive'];
        vars.CFLAGS = '${CFLAGS} -fpermissive';
    end

    cfemmpath = fullfile (thisfilepath, '..', 'cfemm');
    fmesherpath = fullfile (cfemmpath, 'fmesher');
    trianglepath = fullfile (fmesherpath, 'triangle');
    libfemmpath = fullfile (cfemmpath, 'libfemm');
    libluacomplexpath = fullfile (libfemmpath, 'liblua');
    tanglesourcepath = gettanglesourcedir ();
    vars.MEXFLAGS = [vars.MEXFLAGS, ' -DTANGLE_AS_LIBRARY -I"', ...
                     strrep(tanglesourcepath, '\', '/'), '"'];

    [libluacomplex_sources, libluacomplex_headers] = getlibluasources ();

    libfemm_sources = getlibfemmsources ();

    fmesher_sources = { ...
        'fmesher.cpp', ...
        'nosebl.cpp', ...
        'TangleMeshConverter.cpp', ...
        'TangleMesherBackend.cpp', ...
        'TriangleMesher.cpp', ...
        'TriangleMesherBackend.cpp', ...
        'writepoly.cpp', ...
    };

    triangle_sources = {'triangle.c'};
    tangle_sources = {'tangle.cpp', 'float256.cpp'};

    [ libluacomplex_objs,  libluacomplex_rules ] = ...
         mmake.sources2rules ( libluacomplex_sources, ...
                               'SourceDir', libluacomplexpath );

    [ libfemm_objs, libfemm_rules ] = ...
         mmake.sources2rules ( libfemm_sources, ...
                               'SourceDir', libfemmpath );

    [ triangle_objs, triangle_rules ] = ...
         mmake.sources2rules ( triangle_sources, ...
                               'SourceDir', trianglepath );

    [ fmesher_objs, fmesher_rules ] = ...
         mmake.sources2rules ( fmesher_sources, ...
                               'SourceDir', fmesherpath );

    [ tangle_objs, tangle_rules ] = ...
         mmake.sources2rules ( tangle_sources, ...
                               'SourceDir', tanglesourcepath );

    vars.OBJS = [ libluacomplex_objs, ...
                  libfemm_objs, ...
                  triangle_objs, ...
                  tangle_objs, ...
                  fmesher_objs, ...
                  {'mexfmesher.cpp'}, ...
                ];

    if options.DoCrossBuildWin64

        vars = mfemmdeps.mmake_check_cross (options.W64CrossBuildMexLibsDir, vars);

    end

    [rules, vars] = mfemmdeps.mmake_rule_1 ( 'mexfmesher', vars, ...
                                            'DoCrossBuildWin64', options.DoCrossBuildWin64 );
                            
    rules(1).deps = vars.OBJS;

    rules = [ rules, libluacomplex_rules, libfemm_rules, triangle_rules, ...
              tangle_rules, fmesher_rules ];

    rules(end+1).target = 'tidy';
    rules(end).commands = { 'try; delete(''../cfemm/libfemm/liblua/*.${OBJ_EXT}''); catch; end;', ...
                            'try; delete(''../cfemm/libfemm/*.${OBJ_EXT}''); catch; end;', ...
                            'try; delete(''../cfemm/fmesher/*.${OBJ_EXT}''); catch; end;', ...
                            ['try; delete(''', tanglesourcepath, '/*.${OBJ_EXT}''); catch; end;'], ...
                            'try; delete(''*.${OBJ_EXT}''); catch; end;' };
    tidyruleind = numel (rules);

    rules(end+1).target = 'clean';
    rules(end).commands = [ rules(tidyruleind).commands, ...
                            {'try; delete(''*.${MEX_EXT}''); catch; end;'} ];

    % mexfmesher.${MEX_EXT}: ${OBJS}
    %     mex $^ -output $@
%     rules(1).target = {'crossw64'};
%     rules(1).deps = vars.OBJS;
%     rules(1).commands = 'mex ${MEXFLAGS} $^ -output mexfmesher.mexw64';

    if options.Verbose
        fprintf ('In %s, vars contents:\n', mfilename ());
        disp (vars);
        fprintf ('In %s, rules(1) command:\n', mfilename ());
        disp (rules(1).commands);
    end
    
end
