clear all

% Files and directories
% enter directories with trailing slashes .../
indir  = '/home/cengelha/Daten/Aaron/2014-05-29/';
infile  = 'SOMD 5CS 324 359 Apo EK 6us.DSC';

outdir = '/home/cengelha/Daten/Aaron/2014-05-29/';
outfile = 'SOMD 5CS 324 359 Apo EK 6us phase corrected.dat';

% parameters
rotate = false; % rotate an additional 180°, true or false
units = 'deg';  % units of phase angle returned, 'deg' or 'rad'

% load data
[data.x, data.y] = eprload(strcat(indir, infile));

% phase correct it
[ data.ycorr, data.phase, data.offset, data.deviation ] = autophase(data.y, 'rot180', rotate, 'units', units);

% plot it
figure(1)
plot(data.x,real(data.y),data.x,imag(data.y))
figure(2)
plot(data.x,real(data.ycorr),data.x,imag(data.ycorr))

% print it
if units == 'deg'
  fprintf('\nphase: %.2f°\nstandard deviation: %e°\noffset: %e\n', data.phase, data.deviation, data.offset);
else
  fprintf('\nphase: %.3f\nstandard deviation: %e\noffset: %e\n', data.phase, data.deviation, data.offset);
end
% save it
out = [ data.x real(data.ycorr) imag(data.ycorr) ]';

outfile = strcat(outdir, outfile);

fid = fopen(outfile, 'w');
fprintf(fid,'%13.7e %13.7e %13.7e\n',out);
fclose(fid);
