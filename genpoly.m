function [x, y, z, coeff] = genpoly(dim, order, limits, noiselvl, format)
% Generate 1d and 2d polynomial data with random coefficients and noise
%
% Useful for testing fit algorithms etc.
%
% [x, y, z] = genpoly(dim, order)
% [x, y, z, coeff] = genpoly(dim, order, limits, noiselvl, format)
%
% dim:      number of datapoints numel(x) (1d) or [numel(x) numel(y)] (2d).
%           Number of elements determines if 1d or 2d data is generated
% order:    order of the polynomial order(x) or [order(x) order(y)]
% limits:   axis limits [xmin xmax] or [xmin xmax ymin ymax],
%           default: [-1 1 -1 1]
% noiselvl: absolute noise level, default 0.1
% format:   'matrix' or 'vector', format of z data
%           'matrix': z(n,m) = f(x(n),y(m))
%           'vector': z(n) = f(x(n),y(n)
%           default 'vector'

if nargin < 5; format = 'vector'; end
if nargin < 4; noiselvl = 0.1; end
if nargin < 3 && length(dim) > 1; limits = [-1 1 -1 1]; end

% generate random data
x = linspace(limits(1), limits(2) ,dim(1))';
if length(dim) > 1
    y = linspace(limits(3), limits(4), dim(2))';
    xold = x; yold = y;
    x = x*ones(1,dim(2));
    y = ones(dim(1),1)*y';
    x = reshape(x,dim(1)*dim(2),1);
    y = reshape(y,dim(1)*dim(2),1);
else
    order(2) = 0;
    y = ones(length(x),1);
end

z = zeros(length(x),1);

% generate random params
a = 2*rand((max(order) + 1)*(max(order) + 2)/2 - ...
    (abs(order(1)-order(2))*(abs(order(1)-order(2))+1)/2),1) - 1;
poly = cell(length(a),1);

kk = 1;
if order(1) < order(2)
    for ii = 0:order(1)
        for jj = 0:order(2) - ii
            z = z + a(kk)*x.^ii.*y.^jj;
            poly{kk} = strcat('p',num2str(ii),num2str(jj));
            kk = kk+1;
        end
    end
else
    for jj = 0:order(2)
        for ii = 0:order(1) - jj
            z = z + a(kk)*x.^ii.*y.^jj;
            poly{kk} = strcat('p',num2str(ii),num2str(jj));
            kk = kk+1;
        end
    end
end
z = z + noiselvl*(2*rand(length(z),1) - 1);

if length(dim) > 1
    if strcmp(format,'matrix')
        x = xold;
        y = yold;
        z = reshape(z,dim(2),dim(1));
    end
else
    y = z;
end

for ii = 1:length(a)
    coeff.(poly{ii}) = a(ii);
    coeff.nameformat = 'pxorderyorder';
end