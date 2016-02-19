function [xout,yout] = cpgm(infile,varargin)
% pick peaks from CPGM dataset
%
% [xout,yout] = cpgm(infile)
% [xout,yout] = cpgm(infile,window)
% [xout,yout] = cpgm(...,'extrap',<true|false>)
%
% infile: string, filename of CPGM dataset
% window: window in ns over which to integrate the echos, optional
% extrap: extrapolate into noise until end of data or use only detected maxima
%
% if window is not passed, cpgm returns peak values instead

p = inputParser;
p.addRequired('infile', @(x)validateattributes(x,{'char'},{'vector'}));
p.addOptional('window', false, @(x)validateattributes(x,{'numeric'},{'scalar','real'}));
p.addParameter('extrap', false, @(x)validateattributes(x,{'logical'},{'scalar'}));
p.FunctionName = 'cpgm';
parse(p,infile,varargin{:});

% load data
[x,y] = eprload(infile);

% autophase it and make it real
y = autophase(y);
y = real(y);
% find peaks that are above the noise
noise = noiselvl(y);
[maxtab, ~] = peakdet(y,3*noise);

% filter found peaks
% 'real' peaks will have the same distance, so we find the distance that is
% most common
diffvec = diff(maxtab(:,1));
dist = mode(diffvec);
% and record the first two times this distance occurs
firstindex = find(diffvec == dist,2);

% then we calculate the total intensity of the series, with each peak from
% the first to the second one with proper distance
% the series with largest intensity contains the correct peaks
total = -Inf;
if p.Results.extrap
    lastindex = length(x);
else
    lastindex = maxtab(end,1);
end
for ii = 1:firstindex(2)
    currentindex = maxtab(ii,1):dist:lastindex;
    if total < sum(y(currentindex))
        total = sum(y(currentindex));
        index = currentindex;
    end
end

% prepare figures
hMain = findobj('Tag',mfilename);
if ~isempty(hMain); close(hMain); end

hMain = figure('Tag',mfilename);
if ~p.Results.window
    hAxes = axes('Parent',hMain);
    xlabel(hAxes, 'time / ns');
    ylabel(hAxes, 'intensity / a.u.');
else
    
    hAxes(1) = axes('Parent', hMain, 'Tag', 'ax1',...
        'Layer', 'top', ...
        'XAxisLocation', 'Bottom', 'YAxisLocation', 'Left', ...
        'XColor', 'k', 'YColor', 'k');
    xlabel(hAxes(1), 'time / ns');
    ylabel(hAxes(1), 'intensity / a.u.');
    
    hAxes(2) = axes('Parent', hMain, 'Tag', 'ax2', ...
        'Position', get(hAxes(1),'Position'), ...
        'Layer', 'top', ...
        'XAxisLocation', 'Top', 'YAxisLocation', 'Right', ...
        'XTickLabel', '', ...
        'Color', 'none', 'XColor', 'k', 'YColor', [.8 0 0]);
    ylabel(hAxes(2), 'integrated echo intensity / a.u.');
    linkaxes(hAxes);
end

% generate new data vectors and plot result
if ~p.Results.window
    % from maxima
    xout = x(index);
    yout = y(index);
    plot(x,y,xout,yout,'o');
else
    % from integral over y(max +- window/2)
    line('XData', x,'YData', y, ...
        'LineWidth', 0.5, 'LineStyle', '-', 'Color', [0 0 .8], ...
        'Parent', hAxes(1));
    
    xout = x(index);
    yout = zeros(length(index),1);
    for ii = 1:length(index)
        lower = iof(x,x(index(ii)) - p.Results.window/2,'mode','smaller');
        upper = iof(x,x(index(ii)) + p.Results.window/2,'mode','larger');
        yout(ii) = sum(y(lower:upper));
        if yout(ii) > 0
            rectangle('Position',[x(lower),0,x(upper)-x(lower),yout(ii)],...
                'EdgeColor','none', ...
                'FaceColor', [.8 .0 .0], ...
                'Parent', hAxes(2));
        elseif yout(ii) < 0
            rectangle('Position',[x(lower),yout(ii),x(upper)-x(lower),abs(yout(ii))],...
                'EdgeColor','none', ...
                'FaceColor', [.8 .0 .0], ...
                'Parent', hAxes(2));
        end
    end
    xlim(hAxes(2),'auto');
    ylim(hAxes(2),'auto');
end