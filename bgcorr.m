function [dataout, bgout] = bgcorr(x,y,varargin)


%% - TODO -------------------%
%                            %
% - Background indices       %
% - rescale data to order 1  %
% - 2D fit using matrix      %
%----------------------------%


%% ARGUMENT PARSING
% Check number of arguments and set defaults
p = inputParser;
p.addRequired('x', @(x)validateattributes(x,{'numeric'},{'vector','real'}));
p.addRequired('y', @(x)validateattributes(x,{'numeric'},{'vector','real'}));
p.addOptional('z', false, @(x)validateattributes(x,{'numeric'},{'2d','real'}));
p.addParameter('order', 1, @(x)validateattributes(x,{'numeric'},{'nonnegative','integer','vector'}));
p.addParameter('background',false, @(x)validateattributes(x,{'numeric'},{'2d','increasing'}));
p.addParameter('method','mldivide', @(x)ischar(validatestring(x,{'mldivide', 'polyfit', 'fit'})));
p.FunctionName = 'bgcorr';
parse(p,x,y,varargin{:});

% save data dimensions
lenx  = length(x);
sizey = size(y);
leny  = length(y);
z = p.Results.z;
sizez = size(z);
lenz  = max(sizez);
% save min/max values
xmin = min(x);
xrange = max(x) - min(x);
ymin = min(y);
yrange = max(y) - min(y);
zmin = min(min(z));
zrange = max(max(z)) - min(min(z));


% check if we have a z-axis or not
if islogical(z)
    % check arguments for 1D corrections
    % check order
    if length(p.Results.order) ~= 1
        error('MATLAB:bgcorr:parseError','''order'' must have length 1 for 1D correction.')
    end
    if p.Results.order > 9
        error('MATLAB:bgcorr:orderError','''order'' must be < 10 for 1D correction.')
    end
    
    % check x,y
    if lenx ~= leny
        error('MATLAB:bgcorr:dimagree','''x'' and ''y'' must be of same length.')
    end
    if isrow(x); x = x'; end
    if isrow(y); x = y'; end
    
    % check background
    if ~p.Results.background
        % use the default '25% from start/stop' as background.
        bg(1) = xmin;
        bg(2) = xmin + xrange/4;
        bg(3) = xmin + 3*xrange/4;
        bg(4) = xmin + xrange;
    elseif size(p.Results.background) == [4 1]
        bg = p.Results.background;
    else
        error('MATLAB:bgcorr:parseError',['Dimension of ''background'' must be [4 1] for 1D data, but is [', ...
            num2str(size(p.Results.background)), '].'])
    end
    
    % check method
    if strcmp(which(p.Results.method),'')
        error('MATLAB:bgcorr:invalidMethod','Unknown function %s. You might be lacking a toolbox.', ...
            p.Results.method)
    else
        method = p.Results.method;
    end
else
    % check arguments for 2D corrections
    % check order
    if length(p.Results.order) ~= 2
        error('MATLAB:bgcorr:parseError','''order'' must have length 2 for 2D correction.')
    end
    if max(p.Results.order) > 5
        error('MATLAB:bgcorr:orderError','''order'' must be < 6 for 2D correction.')
    end
    
    % check x,y,z
    if (isvector(z) && lenx ~= lenz) || ...
            (min(sizez) > 1 && (sizez(2) ~= lenx || sizez(1) ~= leny))
        error('MATLAB:bgcorr:dimagree','Dimensions of input matrices must agree.')
    end
    if isrow(x); x = x'; end
    if isrow(y); x = y'; end
    if isrow(z); x = z'; end
        
    % check background
    if ~p.Results.background
        % use the default '25% from start/stop' as background.
        bg(1,1) = xmin;
        bg(2,1) = xmin + xrange/4;
        bg(3,1) = xmin + 3*xrange/4;
        bg(4,1) = xmin + xrange;
        bg(1,2) = ymin;
        bg(2,2) = ymin + yrange/4;
        bg(3,2) = ymin + 3*yrange/4;
        bg(4,2) = ymin + yrange;
    elseif size(p.Results.background) == [4 2]
        bg = p.Results.background;
    else
        error('MATLAB:bgcorr:parseError',['Dimension of ''background'' must be [4 2] for 2D data, but is [', ...
            num2str(size(p.Results.background)), '].'])
    end
    
    % check method
    if strcmp(which(p.Results.method),'')
        error('MATLAB:bgcorr:invalidMethod','Unknown function %s. You might be lacking a toolbox.', ...
            p.Results.method)
    elseif strcmp(p.Results.method,'polyfit')
        error('MATLAB:bgcorr:invalidMethod','Invalid method ''%s'' for 2D datasets.', p.Results.method)
    else
        method = p.Results.method;
    end
end

%% PERFORM BACKGROUND CORRECTION
if islogical(z)
    % generate index mask for background region
    indexmask = ((x >= bg(1) & x <= bg(2)) | (x >= bg(3) & x <= bg(4)));
    % normalise data
    xnorm = (x - xmin)/xrange;
    ynorm = (y - ymin)/yrange;
    % 1D background correction
    switch method
        case 'fit'
            poly = strcat('poly',num2str(p.Results.order));
            result = fit(xnorm(indexmask),ynorm(indexmask),poly);
            bgout = result(xnorm) * yrange + ymin;
        case 'polyfit'
            [result, ~, mu] = polyfit(xnorm(indexmask),ynorm(indexmask),p.Results.order);
            bgout = polyval(result,xnorm,[],mu) * yrange + ymin;
        otherwise
            % generate matrix for fit
            Abg = zeros(length(x(indexmask)),p.Results.order+1);
            A   = zeros(length(x),p.Results.order+1);
            % go from highest order down to match polyval function
            for ii = p.Results.order:0
                Abg(:,ii) = xnorm(indexmask).^ii;
                A(:,ii) = xnorm.^ii;
            end
            % solve the system for y
            result = Abg\ynorm(indexmask);
            bgout = A*result * yrange + ymin;
    end
    % reshape result so that it has the same dimensions as the input
    bgout = reshape(bgout,sizey(1),sizey(2));
    y = reshape(y,sizey(1),sizey(2));
    dataout = y - bgout;
else
    % 2D background correction
    % reshape data to vectors if we've been given a matrix
    if min(sizez) > 1
        z = reshape(z,sizez(1)*sizez(2),1);
        x = reshape(x * ones(1,leny) ,sizez(1)*sizez(2),1);
        y = reshape(ones(lenx,1) * y',sizez(1)*sizez(2),1);
    end
    % generate index mask for background region
    indexmask = (((x >= bg(1,1) & x <= bg(2,1)) | (x >= bg(3,1) & x <= bg(4,1)))) & ...
        (((y >= bg(1,2) & y <= bg(2,2)) | (y >= bg(3,2) & y <= bg(4,2))));
    % normalise data
    xnorm = (x - xmin)/xrange;
    ynorm = (y - ymin)/yrange;
    znorm = (z - zmin)/zrange;
    switch method
        case 'fit'
            poly = strcat('poly',num2str(p.Results.order(1)),num2str(p.Results.order(2)));
            result = fit([xnorm(indexmask),ynorm(indexmask)],znorm(indexmask),poly);
            bgout = result(xnorm,ynorm) * zrange + zmin;
        otherwise
            % generate matrix for fit
            Abg = zeros(length(x(indexmask)), ...
                (max(p.Results.order) + 1)*(max(p.Results.order) + 2)/2 - ...
                (abs(p.Results.order(1)-p.Results.order(2))*(abs(p.Results.order(1)-p.Results.order(2))+1)/2));
            A = zeros(length(x), ...
                (max(p.Results.order) + 1)*(max(p.Results.order) + 2)/2 - ...
                (abs(p.Results.order(1)-p.Results.order(2))*(abs(p.Results.order(1)-p.Results.order(2))+1)/2));
            kk = 1;
            if p.Results.order(1) < p.Results.order(2)
                for ii = 0:p.Results.order(1)
                    for jj = 0:p.Results.order(2) - ii
                        Abg(:,kk) = x(indexmask).^ii.*y(indexmask).^jj;
                        A(:,kk) = x.^ii.*y.^jj;
                        kk = kk+1;
                    end
                end
            else
                for jj = 0:p.Results.order(2)
                    for ii = 0:p.Results.order(1) - jj
                        Abg(:,kk) = x(indexmask).^ii.*y(indexmask).^jj;
                        A(:,kk) = x.^ii.*y.^jj;
                        kk = kk+1;
                    end
                end
            end
            % solve the system for z
            result = Abg\znorm(indexmask);
            bgout = A*result * zrange + zmin;
    end
    % reshape result so that it has the same dimensions as the input
    bgout = reshape(bgout,sizez(1),sizez(2));
    z = reshape(z,sizez(1),sizez(2));
    dataout = z - bgout;
end
