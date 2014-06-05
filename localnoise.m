function [noiselvl, noiserms, localstd] = LocalNoise(data, area)
% LocalNoise calculates the local noise level by looking for minima in local
% standard deviations (flat areas of the graph).
%
% Necessary inputs:
% data      - a vector of the data to be analysed
% 
% Optional inputs:
% area      - the area used for local variance, in points. Defaults to 5% of the data.
%
% LocalNoise returns the following:
% noiselvl  - the peak noise lvl
% noiserms  - the rms noise lvl
% localstd  - a vector containing the local standard deviation
% 
% the first and last <area> points of localstd are set to the value of the first
% and last calculated points
%
VERSION = '1.0';

p = inputParser;
p.addRequired('data', @(x)validateattributes(x,{'numeric'},{'vector','real'}));
p.addOptional('area',ceil(length(data)*0.025), @(x)validateattributes(x,{'numeric'},{'nonnegative','scalar','integer'}));
p.FunctionName = 'localnoise';
p.parse(data, area);

% calculate local standard deviation
for i = area+1:length(data)-area
  localstd(i) = std(data(i-area:i+area));
end
% pad localstd so that its length matches data
localstd(1:area) = localstd(area+1);
localstd(end:end+area) = localstd(end);
% calculate local local standard deviation to find flat areas
for i = area+1:length(localstd)-area
  locallocalstd(i) = std(localstd(i-area:i+area));
end
% set standard deviation in flattest area as noiselvl
[~, imin] = min(locallocalstd);
noiserms = localstd(imin);
noiselvl = localstd(imin)*sqrt(2);
