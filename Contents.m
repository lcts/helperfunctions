% Data Evaluation Toolbox
% Version 1.2.0 4-Feb-2016
%
% A collection of functions that automate tasks often needed during data
% analysis/evaluation, like background correction etc.
% 
% Some functions are specific to the needs of EPR spectroscopy
% 
% Functions
%   autophase  - Perform automatic 0th-order phase correction of a complex vector by
%   bgcorr     - Calculate the double integral of a spectrum.
%   db2level   - converts dB into level
%   digitize   - 
%   endoreff   - 
%   fm         - 
%   iof        - returns the index of the element in 'vector' closest to a given
%   level2db   - converts level into dB
%   localnoise - LocalNoise calculates the local noise level by looking for minima in local
%   peakdet    - Detect peaks in a vector
%   sinc       - calculate spectral shape of rectangular pulses, i.e. the sinc function. Real pulses aren't
%   sym        - sym symmetrises data passed to it. By default, the data is symmetrised
