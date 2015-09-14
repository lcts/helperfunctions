function [ data, axis ] = sym(varargin)
% sym symmetrises data passed to it. By default, the data is symmetrised
% around the center, but the symmetrisation point can also be specified.
%
% USAGE:
% symdata = sym(data)
% symdata = sym(data, dim)
% [symdata, symaxis] = sym(data, dim, invpoint, axis)
%
% data:     a vector or 2-dimensional array to be symmetrised
% dim:      the dimension along which to symmetrise. If data is a vector, dim
%           is ignored
% invpoint: the point around which data is symmetrised. Requires axis to be
%           passed to the script as well.
% axis:     the axis vector of of the dimension along which the array
%           should be symmetrised. Has to be ordered.
%
% invpoint and axis:
% If invpoint and axis are passed to the function, the script will first
% generate a new axis that has the same spacing as axis, but is symmetric
% around invpoint. It's maxima/minima will be +-abs(min(axis(1),axis(end)))
% For convenience, this new axis is returned as well
%
p = inputParser;
p.addRequired('data', @(x)validateattributes(x,{'numeric'},{'2d'}));
p.addOptional('dim', 1, @(x)validateattributes(x,{'numeric'},{'scalar'}));
p.addOptional('invpoint', false, @(x)validateattributes(x,{'numeric'},{'scalar'}));
p.addOptional('axis', false, @(x)validateattributes(x,{'numeric'},{'2d'}));

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
        error('sym:OptErr', 'Option invpoint requires option axis');
    else
        if ~issorted(p.Results.axis)
            error('sym:OptErr', 'axis has to be ordered');
        else
            maxval  = min(max(p.Results.axis - p.Results.invpoint), max(-p.Results.axis + p.Results.invpoint));
            spacing = (max(p.Results.axis) - min(p.Results.axis)) / length(p.Results.axis);
            npoints = 2*maxval / spacing;
            axis = linspace(-maxval, maxval, npoints) + p.Results.invpoint;
            
        end
    end
end    

if ~p.Results.invpoint
    if p.Results.axis; warning('sym:OptWarn', 'invpoint not given, option axis ignored'); end
    if isrow(p.Results.data)
        dimension = 2;
    elseif iscolumn(p.Results.data)
        dimension = 1;
    else
        dimension = p.Results.dim;
    end
    data = (p.Results.data + flipdim(p.Results.data, dimension))/2;
    axis = NaN;
else
    if ~p.Results.axis
        error('sym:OptErr', 'Option invpoint requires option axis');
    elseif (isvector(p.Results.data) && ~isvector(p.Results.axis)) || ( ~isvector(p.Results.data) && isvector(p.Results.axis) )
        error('sym:DimMismatch', 'dimensions of data and axis must match');
    elseif isrow(p.Results.data)
        data = interp1(p.Results.axis,p.Results.data,axis);
        data = (data + flipdim(data, 2))/2;
    elseif iscolumn(p.Results.data)
        data = interp1(p.Results.axis,p.Results.data,axis);
        data = (data + flipdim(data, 1))/2;
    else
        if p.Results.dim == 1
            for ii = 1:size(p.Results.data,2)
                data(:,ii) = interp1(p.Results.axis,p.Results.data(:,ii),axis);
                data(:,ii) = (data(:,ii) + flipdim(data(:,ii), 1))/2;
            end
        else
            for ii = 1:size(p.Results.data,1)
                data(ii,:) = interp1(p.Results.axis,p.Results.data(ii,:),axis);
                data(ii,:) = (data(ii,:) + flipdim(data(ii,:), 2))/2;
            end
        end
    end
end