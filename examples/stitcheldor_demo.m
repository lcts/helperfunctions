clear all

% Files and directories
% enter directories with trailing slashes .../
indir  = '';
outdir = '';
outfile = '';

file_short = '';
file_long  = '';

% parameters
autophase = true; % automatically phase-correct data, true or false
offset = 0;       % shift splitpoint towards 0ns
nosave = true;    % don't save result

% load data
[datashort.x, datashort.y] = eprload(strcat(indir, file_short));
[datalong.x,  datalong.y ] = eprload(strcat(indir, file_long ));

% stitch it
[data.x, data.y] = stitchELDOR(datashort.x, datashort.y, datalong.x, datalong.y, 'offset', offset, 'autophase', autophase);

% plot it
figure(1)
plot(data.short.x,real(data.short.y),data.long.x,real(data.long.y))
figure(2)
plot(data.short.x,imag(data.short.y),data.long.x,imag(data.long.y))

figure(3)
hold on
plot(data.stitched.x, real(data.stitched.y), data.interp.x, real(data.interp.y + 0.2 * max(data.interp.y)));
plot(data.stitched.x, imag(data.stitched.y), data.interp.x, imag(data.interp.y + 0.2 * max(data.interp.y)));
hold off

% save it
if ~nosave
	% save interpolated data
	out = [ data.interp.x real(data.interp.y) imag(data.interp.y) ]';
	% save raw stitched data
	% out = [ data.stitched.x real(data.stitched.y) imag(data.stitched.y) ]';

	outfile = strcat(outdir, outfile);

	fid = fopen(outfile, 'w');
	fprintf(fid,'%13.7e %13.7e %13.7e\n',out);
end
fclose(fid);
