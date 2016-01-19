function [ ydata, fwhm ] = sinc(varargin)
% calculate spectral shape of rectangular pulses, i.e. the sinc function. Real pulses aren't
% entirely rectangular and will have less high frequency content
%
% USAGE
% ydata = sinc(xdata, width)
% [ ydata, fwhm ] = sinc(xdata, width, offset)
%
% xdata:   the desired freqency axis, in Hz, required
% width:   the pulselength, in seconds, required
% offset:  the carrier frequency of the pulse, optional, default 0
%
% ydata:   the pulse shape
% fwhm:    the full width at half maximum of the pulse
%
p = inputParser;
p.addRequired('xdata', @(x)validateattributes(x,{'numeric'},{'vector', 'real'}));
p.addRequired('width', @(x)validateattributes(x,{'numeric'},{'scalar', 'real'}));
p.addOptional('offset', 0, @(x)validateattributes(x,{'numeric'},{'scalar', 'real'}));
p.FunctionName = 'sinc';
p.parse(varargin{:});

ydata =  sin(pi * (p.Results.xdata - p.Results.offset) * p.Results.width)./(pi * (p.Results.xdata - p.Results.offset));
fwhm  = 0; % yet to be implemented
end
