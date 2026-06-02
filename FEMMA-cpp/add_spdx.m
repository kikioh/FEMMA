function add_spdx()
% ADD_SPDX  Add SPDX license headers to the FEMMA C++ sources.
%
%   Files under a directory named  mma/  -> GPL-3.0-or-later
%                                           (translated from Svanberg's GPL code)
%   every other .cpp/.hpp/.h         -> MIT
%                                           (original FEMMA translation)
%
%   The script is idempotent: a file that already carries an
%   SPDX-License-Identifier near its top is left untouched, so re-running is
%   safe.  Run it from the C++ project root (the folder with src/ and
%   include/):
%
%       add_spdx

    GPL  = '// SPDX-License-Identifier: GPL-3.0-or-later';
    MIT  = '// SPDX-License-Identifier: MIT';
    exts = {'**/*.cpp','**/*.hpp','**/*.h','**/*.cc','**/*.cxx'};
    skipdirs = {'third_party','build','.git','out'};

    % collect all matching files
    files = [];
    for k = 1:numel(exts)
        files = [files; dir(exts{k})]; %#ok<AGROW>
    end

    n_mit = 0; n_gpl = 0; n_skip = 0;
    for i = 1:numel(files)
        f = files(i);
        fpath = fullfile(f.folder, f.name);

        % path parts relative to the current folder
        rel   = erase(fpath, [pwd filesep]);
        parts = regexp(rel, '[/\\]', 'split');

        % skip vendored / build / version-control directories
        if any(ismember(lower(parts), skipdirs))
            continue;
        end

        % read the file as raw text
        txt   = fileread(fpath);
        lines = splitlines(txt);
        head  = strjoin(lines(1:min(5,numel(lines))), newline);

        % idempotency: skip if already tagged near the top
        if contains(head, 'SPDX-License-Identifier')
            n_skip = n_skip + 1;
            continue;
        end

        % choose the tag: GPL if any parent folder is named "mma"
        dirparts = parts(1:max(end-1,0));
        if any(strcmpi(dirparts, 'mma'))
            tag = GPL; label = 'GPL'; n_gpl = n_gpl + 1;
        else
            tag = MIT; label = 'MIT'; n_mit = n_mit + 1;
        end

        % prepend the tag and rewrite the file
        fid = fopen(fpath, 'w');
        fwrite(fid, [tag newline txt]);
        fclose(fid);

        fprintf('  %s  %s\n', label, rel);
    end

    fprintf('\nDone: %d MIT, %d GPL, %d already tagged (skipped).\n', ...
            n_mit, n_gpl, n_skip);
end
