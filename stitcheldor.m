function [ x, y ] = stitcheldor(varargin)
% stitch together two ELDOR time traces
%
% USAGE:
% [x, y] = stitcheldor(xshort, yshort, xlong, ylong)
% [x, y] = stitcheldor(...,'Option',<value>)
%
% input y data should be complex to allow phase correction,
% otherwise phase-correction will be deactivated
%
% INPUTS:
% x/yshort/long: data to be stitched
%
% Options
% autophase: perform phase-correction, default true
% interp:    return interpolated equispaced x-axes (true) or
%            raw stitched (false)
% offset:    in ns, shift splitpoint manually, default 0
%            'offset' is defined negative, i.e. a positive value shifts the
%            splitpoint towards 0ns
%
p = inputParser;
p.addRequired('xshort', @(x)validateattributes(x,{'numeric'},{'vector'}));
p.addRequired('yshort', @(x)validateattributes(x,{'numeric'},{'vector'}));
p.addRequired('xlong',  @(x)validateattributes(x,{'numeric'},{'vector'}));
p.addRequired('ylong',  @(x)validateattributes(x,{'numeric'},{'vector'}));
p.addParameter('autophase', true, @(x)validateattributes(x,{'logical'},{'scalar'}));
p.addParameter('interp', true, @(x)validateattributes(x,{'logical'},{'scalar'}));
p.addParameter('offset', 0, @(x)validateattributes(x,{'numeric'},{'scalar'}));
p.FunctionName = 'stitcheldor';
p.parse(varargin{:});

% save x axes in output
xshort = p.Results.xshort;
xlong  = p.Results.xlong;

% find splitpoints for short.x and long.x
splitshort = iof(xshort, xshort(end)-p.Results.offset);
splitlong = iof(xlong, xshort(splitshort));
% shift by one if the found index is on the wrong side due to min(abs(...))
if xlong(splitlong) <= xshort(splitshort)
  splitlong = splitlong + 1;
end

% phase-correct y data if needed and data complex
if ~isreal(p.Results.yshort) && p.Results.autophase
  [ yshort, ~, phaseoffset] = autophase(p.Results.yshort);
  yshort = yshort - phaseoffset*1i;
else
  yshort = p.Results.yshort;
end
if ~isreal(p.Results.ylong) && p.Results.autophase
  [ ylong, ~, phaseoffset]  = autophase(p.Results.ylong);
  ylong = ylong - phaseoffset*1i;
else
  ylong  = p.Results.ylong;
  params.long.phase  = false;
end

% interpolate ylong at points xshort
yinterp = interp1(xlong, ylong, xshort);

% fit yinterp to yshort
f = @(x)sum((real(yshort) - (x(1)*real(yinterp) + x(2))).^2);
ab = fminsearch(f,[1 0]);

% scale ylong and stitch it
ystitched = ab(1) * ylong + ab(2);
ystitched = [ yshort(1:splitshort); ystitched(splitlong:end) ];
% generate stitched x axis
xstitched = [ xshort(1:splitshort); xlong(splitlong:end) ];

if p.Results.interp
    % return interpolated data
    xinterp = (xshort(1):xshort(2)-xshort(1):xlong(end))';
    yinterp = interp1(xstitched, ystitched, xinterp);
    y = yinterp;
    x = xinterp;
else
    % return raw data
    y = ystitched;
    x = xstitched;
end