function [dataout, bgout] = bgcorr(x, y, varargin)


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
p.addOptional('z', false, @(x)validateattributes(x,{'numeric'},{'2D','real'}));
p.addParameter('order', 1, @(x)validateattributes(x,{'numeric'},{'nonnegative','integer','vector'}));
p.addParameter('background',false, @(x)validateattributes(x,{'numeric'},{'positive','size',[1,4],'integer','increasing'}));
p.addParameter('method','mldivide', @(x)ischar(validatestring(x,{'mldivide', 'polyfit', 'fit'})));
p.FunctionName = 'bgcorr';
p.parse(x,y,varargin{:});

% save data dimensions
sizex = size(x);
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
zmin = min(z);
zrange = max(max(z)) - min(min(z));


% check if we have a z-axis or not
if islogical(z)
    % check arguments for 1D corrections
    
    % check background
    if ~p.Results.background
        % use the default '25% from start/stop' as background.
        background(1) = 1;
        background(2) = background(1)+ceil(length(x)*0.25);
        background(4) = length(x);
        background(3) = background(4)-ceil(length(x)*0.25);
    elseif p.Results.background <= length(x)
        background = p.Results.background;
    else
        error('MATLAB:bgcorr:parseError',['x has length ', num2str(length(x)), ...
            ', but indices', num2str(p.Results.background), ' were requested.'])
    end
    
    % check order
    if length(order) ~= 1
        error('MATLAB:bgcorr:parseError','''order'' should have one element but has %i.', length(order))
    end
    if order > 9
        error('MATLAB:bgcorr:orderError','''order'' must be < 10 for 1D correction.')
    end
    
    % check x,y
    if lenx ~= leny
        error('MATLAB:bgcorr:dimagree','Dimensions of input matrices must agree.')
    end
    if isrow(x); x = x'; end
    if isrow(y); x = y'; end
    
    % check model
    if strcmp(which(p.Results.model),'')
        error('MATLAB:bgcorr:invalidMethod','Unknown function %s. You might be lacking a toolbox.', ...
            p.Results.model)
    else
        model = p.Results.model;
    end
else
    % check arguments for 2D corrections
    
    % check background
    if ~p.Results.background
        % use the default '25% from start/stop' as background.
        background(1) = 1;
        background(2) = background(1)+ceil(length(x)*0.25);
        background(4) = length(x);
        background(3) = background(4)-ceil(length(x)*0.25);
    elseif p.Results.background <= length(x)
        background = p.Results.background;
    else
        error('MATLAB:bgcorr:parseError',['x has length ', num2str(length(x)), ...
            ', but indices', num2str(p.Results.background), ' were requested.'])
    end
    
    % check order
    if length(order) ~= 2
        error('MATLAB:bgcorr:parseError','''order'' should have one element but has %i.', length(order))
    end
    if max(order) > 5
        error('MATLAB:bgcorr:orderError','''order'' must be < 6 for 2D correction.')
    end
    
    % check x,y,z
    if (isvector(z) && lenx ~= lenz) || ...
            (min(sizez) > 1 && (sizez(1) ~= lenx || sizez(2) ~= leny))
        error('MATLAB:bgcorr:dimagree','Dimensions of input matrices must agree.')
    end
    if isrow(x); x = x'; end
    if isrow(y); x = y'; end
    if isrow(z); x = z'; end
    
    % check model
    if strcmp(which(p.Results.model),'')
        error('MATLAB:bgcorr:invalidMethod','Unknown function %s. You might be lacking a toolbox.', ...
            p.Results.model)
    elseif strcmp(p.Results.model,'polyfit')
        error('MATLAB:bgcorr:invalidMethod','Invalid method ''%s'' for 2D datasets.', p.Results.model)
    else
        model = p.Results.model;
    end
end

%% PERFORM BACKGROUND CORRECTION
if islogical(z)
    % normalise data
    xnorm = (x - xmin)/xrange;
    ynorm = (y - ymin)/yrange;
    % 1D background correction
    switch model
        case 'fit'
            poly = strcat('poly',num2str(order));
            result = fit(xnorm,ynorm,poly);
            bgout = result(xnorm) * yrange + ymin;
        case 'polyfit'
            [result, ~, mu] = polyfit(xnorm,ynorm,order);
            bgout = polyval(result,xnorm,[],mu) * yrange + ymin;
        otherwise
            % generate matrix for fit
            A = zeros(length(x),order+1);
            % go from highest order down to match polyval function
            for ii = order:0
                A(:,ii) = xnorm.^ii;
            end
            % solve the system for y
            result = A\ynorm;
            bgout = A*result * yrange + ymin;
    end
    % reshape result so that it has the same dimensions as the input
    bgout = reshape(bgout,sizey(1),sizey(2));
    y = reshape(y,sizey(1),sizey(2));
    dataout = y - bgout;
else
    % 2D background correction
    % normalise data
    xnorm = (x - xmin)/xrange;
    ynorm = (y - ymin)/yrange;
    znorm = (z - zmin)/zrange;
    % reshape data to vectors if we've been given a matrix
    if min(sizez) > 1
        z = reshape(z,sizez(1)*sizez(2),1);
        x = reshape(x * ones(1,length(y)) ,sizez(1)*sizez(2),1);
        y = reshape(ones(length(x),1) * y',sizez(1)*sizez(2),1);
    end
    switch model
        case 'fit'
            poly = strcat('poly',num2str(order(1)),num2str(order(2)));
            result = fit([x,y],z,poly);
            bgout = result(x,y);
        otherwise
            % generate matrix for fit
            %A = zeros(length(x),order*(order+1)/2+1);
            for ii = 0:order(1)
                for jj = 0:order(2) - ii
                    if order(2) < ii; break; end
                    A(:,ii) = x.^ii.y.^jj;
                end
            end
    end
    % reshape result so that it has the same dimensions as the input
    bgout = reshape(bgout,sizez(1),sizez(2));
    z = reshape(z,sizez(1),sizez(2));
    dataout = z - bgout;
end
