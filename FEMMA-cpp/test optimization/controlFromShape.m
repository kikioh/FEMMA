function control = controlFromShape(control,hdr)
% CONTROLFROMSHAPE  Update a FEMMA control struct from a shape-file header.
%
%   control = controlFromShape(control,hdr)
%
%   Reads the offset grid, power grid, duration/point count, and pulse form
%   from the header struct returned by parseShape and overwrites the matching
%   fields of control.  Spin-system fields (operators, states, type, method,
%   ...) are left untouched, so set those on the incoming control first.

    % number of points and total duration -> uniform pulse_dt
    np     = str2double(hdr.NPOINTS);
    dur_us = str2double(hdr.FEMMA_DURATION_US);
    dt     = dur_us*1e-6/np;
    control.pulse_dt = dt*ones(1,np);

    % offset grid: stored in kHz, control wants Hz
    [lo,hi,n] = parseRange(hdr.FEMMA_OFFSET_RANGE_KHZ);
    control.offsets = {linspace(lo*1e3,hi*1e3,n)};

    % power grid: nominal amplitude (kHz) scaled by a percentage range
    nominal_hz = str2double(hdr.FEMMA_PWR_NOMINAL_KHZ)*1e3;
    [plo,phi,pn] = parseRange(hdr.FEMMA_PWR_RANGE_PCT);
    control.pwr_levels = {2*pi*nominal_hz*linspace(plo/100,phi/100,pn)};

    % pulse form
    if isfield(hdr,'FEMMA_PULSE_FORM')
        control.pulse_form = hdr.FEMMA_PULSE_FORM;
    end
end

function [lo,hi,n] = parseRange(str)
% parse a "lo to hi , N points" range string
    tok = regexp(str,'([-\d.eE+]+)\s*to\s*([-\d.eE+]+)\s*,\s*(\d+)','tokens','once');
    if isempty(tok)
        error('controlFromShape: cannot parse range "%s"',str);
    end
    lo = str2double(tok{1});
    hi = str2double(tok{2});
    n  = str2double(tok{3});
end
