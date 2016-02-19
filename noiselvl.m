function [noise, noiserms, localstd] = noiselvl(data, varargin)
% noiselvl calculates the local noise level by looking for minima in local
% standard deviations (i.e. flat areas of the graph).
%
% Necessary inputs:
% data       - a vector of the data to be analysed
% 
% Optional inputs:
% smoothing  - the area used for local variance, in +-points. Defaults to +-2.5% of the data.
%
% LocalNoise returns the following:
% noise      - the peak noise lvl
% noiserms   - the rms noise lvl
% localstd   - a vector containing the local standard deviation
% 
% the first and last <smoothing> points of localstd are set to the value of the first
% and last calculated points
%

% Check number of arguments and set defaults
p = inputParser;
p.addRequired('data', @(x)validateattributes(x,{'numeric'},{'vector','real'}));
p.addOptional('smoothing',ceil(length(data)*0.025), @(x)validateattributes(x,{'numeric'},{'nonnegative','scalar','integer'}));
p.FunctionName = 'noiselvl';
p.parse(data, varargin{:});

% calculate local standard deviation
localstd = zeros(1,length(data)-p.Results.smoothing);
for i = p.Results.smoothing+1:length(data)-p.Results.smoothing
  localstd(i) = std(data(i-p.Results.smoothing:i+p.Results.smoothing));
end
% pad localstd so that its length matches data
localstd(1:p.Results.smoothing) = localstd(p.Results.smoothing+1);
localstd(end:end+p.Results.smoothing) = localstd(end);
% minima in local noise correspond to flat areas of the data, i.e. contain
% actual noise and not signal
noiserms = min(localstd);
noise = noiserms*sqrt(2);
