function [specs bgs params background] = bgcorr(data, varargin)
% Calculate the double integral of a spectrum.
%
% Syntax
% specs = bgcorr(data)
% specs = bgcorr(data, 'Option', Value, ...)
% [specs bgs params background] = bgcorr(data, ...)
%
% Description
% bgcorr performs polynomial background correction of data, as well as (optionally) (doubly) integrated data.
%
% Parameters & Options
% data       - 2-dimensional array of the data (data(x,y))
% 
% background - vector of indices [leftstart leftstop rightstart rightstop] delimiting the area used
%              for background fit. Take care to include as much background as possible in your spectrum but no signal.
%              If this parameter is not given, DoubleInt will use the left and right 25% of the data as background.
% order      - a vector the orders of the polynomials for background correction. The number of elements determines the number
%              of correction steps (1 or 2). Default [3 3].
% integrate  - number of integration steps, 0-2, default 0
%
% Output
% specs      - an array of the corrected data [xdata ydata firstint secondint]
%
% Additional Outputs
% bgs        - an array of the backgrounds used [xdata firstbg secondbg]
% params     - the coefficients used for background corrections for use in polyval
% background - the background indices

%% ARGUMENT PARSING
% Check number of arguments and set defaults
p = inputParser;
p.addRequired('data', @(x)validateattributes(x,{'numeric'},{'2d','real'}));
p.addParamValue('background',false, @(x)validateattributes(x,{'numeric'},{'positive','size',[1,4],'integer'}));
p.addParamValue('order',[3 3], @(x)validateattributes(x,{'numeric'},{'positive' 'row','integer'}));
p.addParamValue('integrate',0, @(x)validateattributes(x,{'numeric'},{'positive','scalar','integer','<=',2}));
p.FunctionName = 'bgcorr';
p.parse(data,varargin{:});

if ~p.Results.background  
  % use the default '25% from start/stop' as background.
  background(1) = 1;
  background(2) = background(1)+ceil(length(data(:,1))*0.25);
  background(4) = length(data(:,1));
  background(3) = background(4)-ceil(length(data(:,1))*0.25);
else
  background = p.Results.background;
  BGINVALID = false;
  if background(4) > length(data(:,1))
    background(4) = length(data(:,1));
    BGINVALID = true;
  end
  for i = 3:-1:1
    if background(i) > background(i+1)
      background(i) = background(i+1);
      BGINVALID = true;
    end
  end
  if BGINVALID
    warning('bgcorr:BGInvalid','Invalid background. Set to [%i %i %i %i].\n\n', background(1), background(2),background(3),background(4));
  end
end

if length(p.Results.order) > 2
   message = ['order has too many elements. A maximum of two background correction steps are supported.'];
   error('bgcorr:BackgroundSteps', message);   
end


%% INTEGRATE SPECTRUM
% save x-axis for integrated specs and backgrounds
specs(:,1) = data(:,1);
bgs(:,1)   = data(:,1);

% initial background correction
params(:,1) = polyfit(data([background(1):background(2) background(3):background(4)],1), ...
                      data([background(1):background(2) background(3):background(4)],2), p.Results.order(1));
bgs(:,2) = polyval(params(:,1), bgs(:,1));
specs(:,2) = data(:,2) - bgs(:,2);

% first integration step
if p.Results.integrate > 0
    specs(:,3) = cumtrapz(specs(:,1),specs(:,2));
    % if there is a second value in 'order'
    if length(p.Results.order) >= 2
        % perform second bg correction before second integration
        params(:,2) = polyfit(specs([background(1):background(2) background(3):background(4)],1), ...
                              specs([background(1):background(2) background(3):background(4)],2),p.Results.order(2));
        bgs(:,3) = polyval(params(:,2),bgs(:,1));
        specs(:,3) = specs(:,3) - bgs(:,3);
        % then integrate
    end
    % integrate a second time if called for
    if p.Results.integrate > 1
        specs(:,4) = cumtrapz(specs(:,1),specs(:,3));
    end
end
