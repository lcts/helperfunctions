function [ data ] = sym(varargin)
% sym symmetrises data passed to it. By default, the data is symmetrised
% around the center, but the symmetrisation point can also be specified.
%
% USAGE:
% symdata = sym(data)
% symdata = sym(data, dim)
% symdata = sym(data, dim, invpoint axis)
%
% data:     a vector or 2-dimensional array tobe symmetrised
% dim:      the dimension along which to symmetrise. If data is a vector, dim
%           is ignored
% invpoint: the point around which data is symmetrised. Requires axis to be
%           passed to the script as well.
% axis:     a vector of axis values for the dimension along which the array
%           should be symmetrised.
%
% invpoint and axis:
% If invpoint and axis are passed to the function, the script will first
% generate a new axis that has the same spacing as axis, but is symmetric
% around invpoint. It's maxima/minima will be +-abs(min(axis(1),axis(end)))
%
% VERSION 0.9: invpoint and axis not yet implemented!!
%

VERSION = '0.9';


p = inputParser;
p.addRequired('data', @(x)validateattributes(x,{'numeric'},{'2d'}));
p.addOptional('dim', 1, @(x)validateattributes(x,{'numeric'},{'scalar'}));
p.addOptional('invpoint', false, @(x)validateattributes(x,{'numeric'},{'scalar'}));
p.addOptional('axis', false, @(x)validateattributes(x,{'numeric'},{'vector'}));

p.FunctionName = 'sym';
p.parse(varargin{:});

if isrow(p.Results.data)
    dimension = 2;
else
    if iscolumn(p.Results.data)
        dimension = 1;
    else
        dimension = p.Results.dim;
    end
end

if ~p.Results.invpoint
    data = (p.Results.data + flipdim(p.Results.data, dimension))/2;
else
    if ~p.Results.axis
        error('sym:OptErr', 'Option flippoint requires option axis');
    else
        error('sym:OptErr', 'Not yet implemented');
    end
end
        