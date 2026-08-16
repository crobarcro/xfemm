function [objs, deps_rules] = sources2rules (sources, varargin)
% generates rules and object file names from a list of source files
%
% Syntax
%
% [objs, deps_rules] = sources2rules (sources)
% [objs, deps_rules] = sources2rules (..., 'Parameter', value)
%
% Description
%
% sources2rules takes a cell array of source file names and creates a cell
% array of object file targets and an array of structures containing
% dependancy rules for the source file headers, suitable for use by
% mmake.make to build the files. sources2rules looks for header files with
% the same base name as the source file but with the file extension '.h' or
% '.hpp', if it finds one, it adds an mmake rule to make the object file
% dependant on the header file. For example, if the supplied source file
% is 'mysource.cpp', and a header file named msource.h' is found, a rule
% like the following is added:
%
% rule.target = 'mysource.${OBJ_EXT}';
% rule.deps = 'mysource.h';
% rule.commands = []; % no commands associated with this rule
%
% Input
%
%  sources - cell array of strings containing the source files for which
%    object file name rules and dependancy rules will be created.
%
% Additional optional arguments may be supplied using Parameter-Value
% pairs, the available options are:
%
% 'SourceDirectory' - directory containing source files. If not supplied
%   the current working directory is assumed.
%
% 'HeaderDirectory' - directory containing header files. If not supplied
%   the same directory as SourceDirectory is assumed.
%
% Output
%
%  objs - cell array of strings containing the source file names, but with
%    the file extension replaced with '.${OBJ_EXT}'.
%
%  deps_rules - array of structures containing dependancy rules for the
%    any header files found for the corresponding source file. 
%
% See Also: 
%

    options.SourceDirectory = '';
    options.HeaderDirectory = '';
    
    options = mmake.parse_pv_pairs (options, varargin);
    
    if isempty (options.HeaderDirectory)
        options.HeaderDirectory = options.SourceDirectory;
    end

    objs = {};
    deps_rules = struct ('target', '', 'deps', '', 'commands', '');
    valid_header_suffix = {'.h', '.hpp'};
    valid_source_suffix = {'.cpp', '.cc', '.cxx', '.c', '.c++'};
    
    deps_rules_ind = 1;
    for objsind = 1:numel (sources)
        
        [sourcedir,name,ext] = fileparts ( sources{objsind} );
        
        if any (strcmpi (ext, valid_source_suffix))
            
            sourcefile = fullfile (options.SourceDirectory, sources{objsind});
            
            if ~exist (sourcefile, 'file')
                error ('Source file:\n%s\ndoes not exist.', sourcefile);
            end
            
            objs = [objs, { fullfile(options.SourceDirectory, sourcedir, [name, '.${OBJ_EXT}']) } ];
            
            for hsuffixind = 1:numel (valid_header_suffix)
                
                header_file = fullfile (options.HeaderDirectory, sourcedir, [name, valid_header_suffix{hsuffixind}]);
                
                if exist (header_file, 'file')

                    deps_rules(deps_rules_ind).target = objs{objsind};
                    deps_rules(deps_rules_ind).deps = header_file;
                    deps_rules(deps_rules_ind).commands = [];
                    deps_rules_ind = deps_rules_ind + 1;
                    break;
                    
                end
                
            end
        end
        
    end
    
end