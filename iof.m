function [ index, element ] = iof(varargin)
% returns the index of element in vector closest to value
%
% USAGE:
% index = iof(vector, value)
% [index, element] = iof(vector, value, mode)
%
% mode can be one of
% 'closest': closest element to value (default)
% 'smaller': closest element smaller than value
% 'larger':  closest element larger than value
%
VERSION = '1.0';

p = inputParser;
p.addRequired('vector', @(x)validateattributes(x,{'numeric'},{'vector', 'real'}));
p.addRequired('value', @(x)validateattributes(x,{'numeric'},{'scalar', 'real'}));
p.addOptional('mode', 'closest', @(x)ischar(validatestring(x,{'closest', 'smaller', 'larger'})));
p.FunctionName = 'iof';
p.parse(varargin{:});

% find index of element in vector closest to value
[~, index] = min(abs(p.Results.vector - p.Results.value));

switch p.Results.mode
  case 'closest'
    element = p.Results.vector(index);
  case 'smaller'
    if p.Results.vector(index) >= p.Results.value; index = index - 1; end
    element = p.Results.vector(index);
  case 'larger'
    if p.Results.vector(index) <= p.Results.value; index = index + 1; end
    element = p.Results.vector(index);
end

