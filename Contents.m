% Data Evaluation Toolbox
% Version 1.2.1 19-Feb-2016
%
% A collection of functions that automate tasks often needed during data
% analysis/evaluation, like background correction etc.
% 
% Some functions are specific to the needs of EPR spectroscopy
% 
% Functions
%   autophase - Perform automatic 0th-order phase correction of a complex vector by
%   bgcorr    - Perform background correction
%   db2level  - Convert dB into level
%   level2db  - Convert level into dB
%   digitize  - Get plotted data from an image file
%   endoreff  - Calculate ENDOR efficiency
%   fm        - Encode signal intensity as axis point density using frequency modulation
%   iof       - Get the index of the element in 'vector' closest to a given value
%   noiselvl  - Calculate the local noise level
%   peakdet   - Detect peaks in a vector
%   sinc      - Calculate spectral shape of rectangular pulses, i.e. the sinc function.
%   sym       - Symmetrises data
%   genpoly   - Generate 1d and 2d polynomial data with random coefficients and noise
%   cpgm      - Extract decay data from Carr-Purcell Gill-Meiboom measurements
