function factor = endoreff(varargin)
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
VERSION = '0.9';

p = inputParser;
p.addRequired('coupling', @(x)validateattributes(x,{'numeric'},{'vector', 'real'}));
p.addRequired('tau', @(x)validateattributes(x,{'numeric'},{'scalar', 'real'}));
p.FunctionName = 'endoreff';
p.parse(varargin{:});

% eta = p.Results.coupling .* (p.Results.tau) / (2 * pi));
% factor = sqrt(2)*eta./(eta.^2 + 1/2);
eta = p.Results.coupling* 2 * p.Results.tau;
factor = abs(1.4*eta./(eta.^2 + 0.49));
end
