function [ data, phase, offset, deviation ] = autophase(varargin)
% Perform automatic 0th-order phase correction of a complex vector by
% minimizing the signal in the imaginary channel
%
% USAGE:
% dataout = autophase(data)
% [dataout, phase] = autophase(data)
% [dataout, phase, offset, deviation] = autophase(data)
% [dataout, phase, offset, deviation] = autophase(data, 'units', '<value>', 'rot180', <boolean>)
%
% data:       complex data vector to be phase-corrected
% units:      'rad' or 'deg', return phase in degrees or rad
% rot180:     true or false, rotate an additional 180°. The script automatically adjusts
%             the phase so that real(data) has a positive integral. Set this to true if auto-
%             adjustment fails
%
% data:       phase-corrected data
% phase:      the phase used for correction
% offset:     the 0th order polynomial 
% deviation:  the deviation of the imaginary part from 0th order polynomial, normalized
%
p = inputParser;
p.addRequired('data', @(x)validateattributes(x,{'numeric'},{'vector'}));
p.addParamValue('units', 'rad', @(x)ischar(validatestring(x,{'rad', 'deg'})));
p.addParamValue('rot180', false, @(x)validateattributes(x,{'logical'},{'scalar'}));
p.FunctionName = 'autophase';
p.parse(varargin{:});

VERSION = '0.8';
fprintf('\nautophase v%s\n', VERSION);

% function for phase correction: Minimize signal intensity in imag. channel:
% Multiply data with phase angle:                     data*exp(i*x(1))
% take the imaginary part and 0-order bg correct it:  imag(...) - x(2)
% square it element-wise:                             (...).^2
% and sum over the resulting vector:                  sum(...)
% Define that as a function of phi:                   f = @(phi)...
f = @(x)sum((imag(p.Results.data * exp(i*x(1))) - x(2)).^2);

% find the minimum deviation from zero 
[ phase, deviation ] = fminsearch(f, [0 0]);
% save values and correct phases
offset = phase(2);
phase  = phase(1);
data   = p.Results.data * exp(i*phase);

% rotate 180° if necessary
% rotation not necessary when neither (rotate 0°) or both (rotate 360°) rot180 and
% trapz(real(data)) < 0 are true
if xor(trapz(real(data)) < 0, p.Results.rot180)
  data  = data * exp(i*pi);
  phase = phase + pi;
end

deviation = sqrt(deviation)/(length(data)*max(abs(data)));
if p.Results.units == 'deg'; phase = phase/pi*180; end
