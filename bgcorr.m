function [dataout, bgout] = bgcorr(x,y,varargin)
% Perform background correction
%
% bgcorr performs automatic background correction of 1d or 2d real data using
% polynomial functions with user-specifiable order, background area and
% method.
%
% Syntax
% dataout = bgcorr(x,y)
% dataout = bgcorr(x,y,z)
% [dataout, bgout] = bgcorr(x,y,'order',<value>,'background',<value>,'method',<value>)
% [dataout, bgout] = bgcorr(x,y,z,'order',<value>,'background',<value>,'method',<value>)
% [dataout, bgout] = bgcorr(...,'method','fit','opts',<fitopts>)
%
% Description
% dataout = bgcorr(x,y) returns vector y background-corrected using a
% poynomial of default order (1) fitted to the default area (left and right
% 25% of data)
%
% dataout = bgcorr(x,y,z) does the same, but for 2d data z, using order 
% [1 1]. if z is a matrix, the data is assumed to be formatted as
% z(n) = f(x(n),y(n)). If z is a matrix, the data is assumed to be
% formatted as z(n,m) = f(x(n), y(m)). Default order is [1 1], and default
% background area is the outer 25% in both dimensions.
%
% [dataout, bgout] = bgcorr(x,y,'order',<value>,'background',<value>,'method',<value>)
% [dataout, bgout] = bgcorr(x,y,z,'order',<value>,'background',<value>,'method',<value>)
% return background corrected data and the fitted background itself, using
% order, background and method specified by the user.
% - 'order' must be scalar and <10 for 1d and two-element and <6 for 2d data.
%    The elements denote polynomial order in x and y direction, respectively.
% - 'background' is formatted as [lowstart lowstop highstart highstop] for 1d
%    and [xlstart xlstop xhstart xhstop ylstart ylstop yhstart yhstop] for 2d
% - 'method' can be one of 'mldivide' (default), 'polyfit' (only for 1d) and 'fit'
%   - 'mldivide' uses the Matlab-internal mldivide function and requires
%      no toolboxes.
%   - 'polyfit' uses the Optimization Toolbox's 'polyfit' function. Works
%      only for 1d, but is slightly faster than 'mldivide'
%   - 'fit' uses the Curve-Fitting Toolbox's 'fit' function. It is
%      signifcantly slower, but allows advanced control over the fitting
%      process
% It is unlikely that setting this to anything other than the default will
% gain you anything
%
% [dataout, bgout] = bgcorr(...,'method','fit','opts',<fitopts>)
% Perform correction using the Curve-Fitting Toolbox's 'fit' function.
% 'opts' is passed to the function. See 'help fit' for details.
%

%% INFORM USER ABOUT SYNTAX CHANGE
% if user passed a 2xN or Nx2 array as x, they are likely using the syntax of
% the previous bgcorr version.
if isnumeric(x) && min(size(x)) == 2
    message = strcat('Invalid input format. ', ...
        '\nYou are likely using a deprecated input syntax (data passed as one Nx2 array). ', ...
        '\n\nThe current version of bgcorr has changed it''s syntax to accomodate correction of 2d datasets.', ...
        ' Please adjust your scripts accordingly. Type ''help bgcorr'' for the new syntax.');
    error('MATLAB:bgcorr:deprecated', message);
end


%% ARGUMENT PARSING & VALIDATION
% Perform basic validity checks and set defaults
p = inputParser;
p.addRequired('x', @(x)validateattributes(x,{'numeric'},{'vector','real'}));
p.addRequired('y', @(x)validateattributes(x,{'numeric'},{'vector','real'}));
p.addOptional('z', false, @(x)validateattributes(x,{'numeric'},{'2d','real'}));
p.addParameter('order', 1, @(x)validateattributes(x,{'numeric'},{'nonnegative','integer','vector'}));
p.addParameter('background',false, @(x)validateattributes(x,{'numeric'},{'vector'}));
p.addParameter('method','mldivide', @(x)ischar(validatestring(x,{'mldivide', 'polyfit', 'fit'})));
p.FunctionName = 'bgcorr';
parse(p,x,y,varargin{:});

% save data dimensions, b/c we might need to restore them later
z = p.Results.z;
sizey = size(y);
sizez = size(z);
% save min/max values so that we can restore them after normalisation
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
    validateattributes(p.Results.order,{'numeric'},{'scalar','<',10},'bgcorr','order for 1d data')
   
    % check x,y dimensions
    if numel(x) ~= numel(y)
        error('MATLAB:bgcorr:dimagree','''x'' and ''y'' must be of same length.')
    end
    % force columns
    if isrow(x); x = x'; end
    if isrow(y); x = y'; end
    
    % check background
    if ~p.Results.background
        % use the default '25% from start/stop' as background.
        bg(1) = xmin;
        bg(2) = xmin + xrange/4;
        bg(3) = xmin + 3*xrange/4;
        bg(4) = xmin + xrange;
    else
        % error if background does not have 4 increasing elements
        validateattributes(p.Results.background,{'numeric'},{'vector','increasing','numel',4},'bgcorr','background for 1d data')
        bg = p.Results.background;
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
    validateattributes(p.Results.order,{'numeric'},{'vector','numel',2,'<',6},'bgcorr','order for 2d data')
    
    % check x,y,z dimensions
    if (isvector(z) && numel(x) ~= numel(z)) || ...
            (min(sizez) > 1 && (sizez(2) ~= numel(x) || sizez(1) ~= numel(y)))
        error('MATLAB:bgcorr:dimagree','Dimensions of input matrices must agree.')
    end
    % force columns
    if isrow(x); x = x'; end
    if isrow(y); x = y'; end
    if isrow(z); x = z'; end
        
    % check background
    if ~p.Results.background
        % use the default '25% from start/stop' as background.
        bg(1) = xmin;
        bg(2) = xmin + xrange/4;
        bg(3) = xmin + 3*xrange/4;
        bg(4) = xmin + xrange;
        bg(5) = ymin;
        bg(6) = ymin + yrange/4;
        bg(7) = ymin + 3*yrange/4;
        bg(8) = ymin + yrange;
    else
        % error if background does not have 8 elements or if 1-4 (x) / 5-8
        % (y) are not increasing
        validateattributes(p.Results.background,{'numeric'},{'vector','numel',8},'bgcorr','background for 2d data')
        validateattributes(p.Results.background(1:4),{'numeric'},{'vector','increasing'},'bgcorr','elements 1:4 (x) of background')
        validateattributes(p.Results.background(5:8),{'numeric'},{'vector','increasing'},'bgcorr','elements 5:8 (y) of background')
        bg = p.Results.background;
    end
    
    % check method
    if strcmp(which(p.Results.method),'')
        error('MATLAB:bgcorr:invalidMethod','Unknown function %s. You might be lacking a toolbox.', ...
            p.Results.method)
    elseif strcmp(p.Results.method,'polyfit') % polyfit does not work for 2d data.
        error('MATLAB:bgcorr:invalidMethod','Invalid method ''%s'' for 2d data.', p.Results.method)
    else
        method = p.Results.method;
    end
end

%% BACKGROUND CORRECTION

%% 1d correction
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
        otherwise % mldivide is the default
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
    
%% 2d correction
else
    % reshape data to vectors if we've been given a matrix
    if min(sizez) > 1
        [x,y] = ndgrid(x,y);
        x = reshape(x,numel(z),1);
        y = reshape(y,numel(z),1);
        z = reshape(z,numel(z),1);
    end
    % generate index mask for background region
    indexmask = (((x >= bg(1) & x <= bg(2)) | (x >= bg(3) & x <= bg(4)))) & ...
                (((y >= bg(5) & y <= bg(6)) | (y >= bg(7) & y <= bg(8))));
    % normalise data
    xnorm = (x - xmin)/xrange;
    ynorm = (y - ymin)/yrange;
    znorm = (z - zmin)/zrange;
    switch method
        case 'fit'
            poly = strcat('poly',num2str(p.Results.order(1)),num2str(p.Results.order(2)));
            result = fit([xnorm(indexmask),ynorm(indexmask)],znorm(indexmask),poly);
            bgout = result(xnorm,ynorm) * zrange + zmin;
        otherwise % mldivide is the default
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
