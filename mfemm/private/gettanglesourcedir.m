function sourcedir = gettanglesourcedir ()
%GETTANGLESOURCEDIR Locate the pinned Tangle sources used by MEX builds.
% TANGLE_SOURCE_DIR can name an existing checkout for offline builds. Otherwise
% the exact revision is cached in the system temporary directory.

    revision = 'a808c624ec0584569e43593662f54890b602c6af';
    sourcedir = getenv ('TANGLE_SOURCE_DIR');
    configured = ~isempty (sourcedir);
    if ~configured
        sourcedir = fullfile (tempdir (), ['xfemm-tangle-', revision]);
    end
    required = {'tangle.cpp', 'float256.cpp', 'tangle_mesh.h'};
    if checkout_matches (sourcedir, revision, required)
        return;
    end
    if configured
        error ('xfemm:tangle:InvalidSource', ...
               'TANGLE_SOURCE_DIR is not a Tangle checkout at revision %s.', revision);
    end
    if exist (sourcedir, 'dir')
        rmdir (sourcedir, 's');
    end
    quoted = shellquote (sourcedir);
    command = sprintf (['git clone --quiet https://github.com/dcm3c/tangle.git %s', ...
                       ' && git -C %s checkout --quiet %s'], quoted, quoted, revision);
    [status, output] = system (command);
    if status ~= 0 || ~checkout_matches (sourcedir, revision, required)
        error ('xfemm:tangle:FetchFailed', ...
               'Could not retrieve pinned Tangle revision %s:\n%s', revision, output);
    end
end

function matches = checkout_matches (sourcedir, revision, required)
    matches = exist (sourcedir, 'dir') == 7;
    for ind = 1:numel (required)
        matches = matches && exist (fullfile (sourcedir, required{ind}), 'file') == 2;
    end
    if ~matches
        return;
    end
    [status, output] = system (sprintf ('git -C %s rev-parse HEAD', shellquote (sourcedir)));
    matches = status == 0 && strcmpi (strtrim (output), revision);
end

function quoted = shellquote (value)
    if ispc
        quoted = ['"', strrep(value, '"', '""'), '"'];
    else
        quoted = ['''', strrep(value, '''', '''"''"'''), ''''];
    end
end
