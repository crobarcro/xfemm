function sourcedir = gettanglesourcedir ()
%GETTANGLESOURCEDIR Locate the pinned Tangle sources used by MEX builds.
% TANGLE_SOURCE_DIR can name an existing checkout for offline builds. Otherwise
% the exact revision is cached in the system temporary directory. The repository
% and revision can be overridden with the XFEMM_TANGLE_REPOSITORY and
% XFEMM_TANGLE_REVISION environment variables so CI pins one version everywhere.

    % Development fork carrying the library options and boundary-match API.
    % Revert to https://github.com/dcm3c/tangle.git once those changes are
    % upstreamed, updating the revision to the merged upstream commit.
    repository = getenv ('XFEMM_TANGLE_REPOSITORY');
    if isempty (repository)
        repository = 'https://github.com/crobarcro/tangle.git';
    end
    revision = getenv ('XFEMM_TANGLE_REVISION');
    if isempty (revision)
        revision = 'b5d51ad8639526a83576b42aee5ce1ef82d96e9b';
    end
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
    command = sprintf (['git clone --quiet %s %s', ...
                       ' && git -C %s checkout --quiet %s'], ...
                       repository, quoted, quoted, revision);
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
