function [rules,vars] = MMakefile_fpproc (varargin)

    options.Verbose = false;
    options.Debug = false;
    options.DoCrossBuildWin64 = false;
    options.W64CrossBuildMexLibsDir = '';
    options.DebugSymbols = false;
    options.ExtraMEXFLAGS = '';

    options = mmake.parse_pv_pairs (options, varargin);

    if options.Debug
        options.DebugSymbols = true;
    end

%     mfemmdeps.getmfilepath (mfilename);

%     if ispc
%         trilibraryflag = '-DCPU86';
%     else
%         trilibraryflag = '-DLINUX';
%     end

    % flags that will be passed direct to mex
    %vars.MEXFLAGS = '${MEXFLAGS} -I"postproc" -I"../cfemm/fpproc" -I"../cfemm/libfemm" -I"../cfemm/libfemm/liblua" ';
    vars.MEXFLAGS = ['${MEXFLAGS} -D_GLIBCXX_USE_CXX11_ABI=1 -I"postproc" -I"../cfemm/fpproc" -I"../cfemm/libfemm" -I"../cfemm/libfemm/liblua" ', options.ExtraMEXFLAGS];
    
    if options.Debug
        vars.MEXFLAGS = [vars.MEXFLAGS, ' -DDEBUG '];
    else
        vars.MEXFLAGS = [vars.MEXFLAGS, ' -DNDEBUG '];
    end
    
    if isunix && ~mfemmdeps.isoctave ()
        if options.Debug
            vars.OPTIMFLAGS = '-OO';
            vars.MEXFLAGS = [vars.MEXFLAGS, ' CXXOPTIMFLAGS="-O0" '];
        else
            vars.OPTIMFLAGS = '-O2';
            vars.MEXFLAGS = [vars.MEXFLAGS, ' CXXOPTIMFLAGS="-O2" '];
        end
    end

    if options.Verbose
        vars.MEXFLAGS = [vars.MEXFLAGS, ' -v '];
    end

    if options.DebugSymbols
        vars.MEXFLAGS = [vars.MEXFLAGS, ' -g'];
    end

    vars.CXXFLAGS = '${CXXFLAGS} -std=c++17 ';

    if mfemmdeps.isoctave
        setenv('CXXFLAGS','-std=c++17');
    end


    %vars.LDFLAGS = '${LDFLAGS} -lstdc++ ''-Wl,--no-undefined''';
    vars.LDFLAGS = '${LDFLAGS} ''-Wl,--no-undefined''';

    [libluacomplex_sources, libluacomplex_headers] = getlibluasources ();

    libfemm_sources = getlibfemmsources ();

    fpproc_sources = { ...
        'fpproc.cpp', ...
        'makemask.cpp', ...
        'CMPointVals.cpp', ...
        'CPostProcMElement.cpp', ...
    };

    [ libluacomplex_objs,  libluacomplex_rules ] = ...
         mmake.sources2rules ( libluacomplex_sources, ...
                               'SourceDir', fullfile('..', 'cfemm', 'libfemm', 'liblua') );

    [ libfemm_objs, libfemm_rules ] = ...
         mmake.sources2rules ( libfemm_sources, ...
                               'SourceDir', fullfile('..', 'cfemm', 'libfemm') );

    [ fpproc_objs, fpproc_rules ] = ...
         mmake.sources2rules ( fpproc_sources, ...
                               'SourceDir', fullfile('..', 'cfemm', 'fpproc') );

    vars.OBJS = [ libluacomplex_objs, ...
                  libfemm_objs, ...
                  fpproc_objs, ...
                  {'postproc/fpproc_interface.${OBJ_EXT}', 'postproc/fpproc_interface_mex.${OBJ_EXT}' }
                ];

    if options.DoCrossBuildWin64

        vars = mfemmdeps.mmake_check_cross (options.W64CrossBuildMexLibsDir, vars);

    end

    [rules, vars] = mfemmdeps.mmake_rule_1 ( 'fpproc_interface_mex', vars, ...
                                             'DoCrossBuildWin64', options.DoCrossBuildWin64, ...
                                             'AddStdCppLib', true);
    rules(1).deps = vars.OBJS;

    rules = [ rules, libluacomplex_rules, libfemm_rules, fpproc_rules ];

    rules(end+1).target = 'postproc/fpproc_interface.${OBJ_EXT}';
    rules(end).deps = 'postproc/fpproc_interface.h';

    rules(end+1).target = 'postproc/fpproc_interface_mex.${OBJ_EXT}';
    rules(end).deps = 'postproc/class_handle.hpp';


    rules(end+1).target = 'tidy';
    rules(end).commands = { 'try; delete(''../cfemm/libfemm/liblua/*.${OBJ_EXT}''); catch; end;', ...
                            'try; delete(''../cfemm/libfemm/*.${OBJ_EXT}''); catch; end;', ...
                            'try; delete(''../cfemm/fpproc/*.${OBJ_EXT}''); catch; end;', ...
                            'try; delete(''postproc/*.${OBJ_EXT}''); catch; end;', ...
                            'try; delete(''*.${OBJ_EXT}''); catch; end;'};
    tidyruleind = numel (rules);

    rules(end+1).target = 'clean';
    rules(end).commands = [ rules(tidyruleind).commands, ...
                            {'try; delete(''*.${MEX_EXT}''); catch; end;'}];

                        
    if options.Verbose
        fprintf ('In %s, vars contents:\n', mfilename ());
        disp (vars);
        fprintf ('In %s, rules(1) command:\n', mfilename ());
        disp (rules(1).commands);
    end
end


