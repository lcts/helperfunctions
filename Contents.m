% Data Evaluation Toolbox
% Version 1.2.2 19-Feb-2016
%
% A collection of functions that automate tasks often needed during data
% analysis/evaluation, like background correction etc.
% 
% Some functions are specific to the needs of EPR spectroscopy
% 
% Functions
%
% Convenience Functions
%   db2level  - Convert dB into level
%   level2db  - Convert level into dB
%   iof       - Get the index of the element in 'vector' closest to a given value
%   genpoly   - Generate 1d and 2d polynomial data with random coefficients and noise
%
% Data treatment and analysis functions
%   autophase - Perform automatic 0th-order phase correction of complex data
%   bgcorr    - Perform background correction
%   digitize  - Get plotted data from an image file
%   fm        - Encode signal intensity as axis point density using frequency modulation
%   noiselvl  - Calculate noise level and noise-suppressed pseudo-derivative
%   peakdet   - Detect peaks in a vector
%   sym       - Symmetrise data
%
% EPR-specific functions
%   endoreff    - Calculate ENDOR efficiency
%   sinc        - Calculate spectral shape of rectangular pulses, i.e. the sinc function.
%   cpgm        - Extract decay data from Carr-Purcell Gill-Meiboom measurements
%   stitcheldor - Combine two ELDOR traces of different length into one
