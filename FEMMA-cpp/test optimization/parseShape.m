function [amp,pha,hdr] = parseShape(filename)
% PARSESHAPE  Read a Bruker JCAMP-DX shape file written by FEMMA.
%
%   [amp,pha,hdr] = parseShape(filename)
%
%   amp  - amplitude column, per cent (constant 100 for a phase pulse)
%   pha  - phase column, degrees
%   hdr  - struct of header fields (e.g. hdr.NPOINTS, hdr.FEMMA_PULSE_FORM,
%          hdr.FEMMA_FIDELITY)
%
%   Header lines start with '##'; data lines are "amplitude, phase" pairs
%   between '##XYPOINTS=' and '##END='.

    fid = fopen(filename,'r');
    if fid<0, error('parseShape: cannot open %s',filename); end

    amp = []; pha = []; hdr = struct();
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line), continue; end
        if startsWith(line,'##END'), break; end
        if startsWith(line,'##')
            % header line of the form "##KEY= value" or "##$KEY= value"
            kv = regexp(line,'^##\$?([A-Z0-9_]+)=\s*(.*)$','tokens','once');
            if ~isempty(kv)
                hdr.(kv{1}) = strtrim(kv{2});
            end
            continue;
        end
        % data line "amplitude, phase"
        parts = str2double(strsplit(line,','));
        if numel(parts)==2 && all(~isnan(parts))
            amp(end+1,1) = parts(1); %#ok<AGROW>
            pha(end+1,1) = parts(2); %#ok<AGROW>
        end
    end
    fclose(fid);

    % sanity check against the declared point count
    if isfield(hdr,'NPOINTS')
        np = str2double(hdr.NPOINTS);
        if ~isnan(np) && numel(pha)~=np
            warning('parseShape: header NPOINTS=%d but read %d data points', ...
                    np, numel(pha));
        end
    end
end
