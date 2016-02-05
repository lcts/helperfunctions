function [x, y, z] = genpoly(dim, order, noiselvl, format)

if nargin < 4
    format = 'vector';
end
% generate random data
limits = sortrows(rand(2));
x = linspace(100*(2*limits(1,1)-1), 100*(2*limits(2,1)-1),dim(1))';
if length(dim) > 1
    y = linspace(100*(2*limits(1,2)-1), 100*(2*limits(2,2)-1),dim(2))';
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

kk = 1;
if order(1) < order(2)
    for ii = 0:order(1)
        for jj = 0:order(2) - ii
            z = z + a(kk)*x.^ii.*y.^jj;
            kk = kk+1;
        end
    end
else
    for jj = 0:order(2)
        for ii = 0:order(1) - jj
            z = z + a(kk)*x.^ii.*y.^jj;
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