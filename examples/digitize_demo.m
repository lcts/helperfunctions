image_file = 'xepr-screenshot.png';

% read image as NxMx(3 uint8) RGB matrix
image_matrix = imread(image_file);
% get the data
% follow the instructions given in the command window
% selected positions can be fine-tuned using the arrow keys before
% confirming by hitting 'Return'/'Enter'
data = digitize(image_matrix);