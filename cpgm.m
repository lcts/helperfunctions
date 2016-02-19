function [xout,yout] = cpgm(infile,window)
% pick peaks from CPGM dataset
%
% [xout,yout] = cpgm(infile)
% [xout,yout] = cpgm(infile,window)
%
% infile: Filename of CPGm dataset
% window: window in ns over which to integrate the echos, optional
%
% if window is not passed, cpgm returns peak values instead

if nargin < 2
    window = false;
end

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
firstfit = find(diffvec == dist,2);

% then we calculate the total intensity of the series, with each peak from
% the first to the second one with proper distance
% the series with largest intensity contains the correct peaks
total = -Inf;
for ii = 1:firstfit(2)
    index = maxtab(ii,1):dist:maxtab(end,1);
    if total < sum(y(index))
        total = sum(y(index));
        bestindex = index;
    end
end

% prepare figures
hMain = findobj('Tag',mfilename);
if isempty(hMain)
    hMain = figure('Tag',mfilename);
    if ~window
        hAxes = axes('Parent',hMain);
        xlabel(hAxes, 'time / ns');
        ylabel(hAxes, 'intensity / a.u.');
    else
        
        hAxes(1) = axes('Parent', hMain, ...
            'Layer', 'top', ...
            'XAxisLocation', 'Bottom', 'YAxisLocation', 'Left', ...
            'XColor', 'k', 'YColor', 'k');
        xlabel(hAxes(1), 'time / ns');
        ylabel(hAxes(1), 'intensity / a.u.');
        
        hAxes(2) = axes('Parent', hMain, ...
            'Position', get(hAxes(1),'Position'), ...
            'Layer', 'top', ...
            'XAxisLocation', 'Top', 'YAxisLocation', 'Right', ...
            'XTickLabel', '', ...
            'Color', 'none', 'XColor', 'k', 'YColor', [.8 0 0]);
        ylabel(hAxes(2), 'integrated echo intensity / a.u.');
        linkaxes(hAxes);
    end
else
    figure(hMain);
end

% generate new data vectors and plot result
if ~window
    % from maxima
    xout = x(bestindex);
    yout = y(bestindex);
    plot(x,y,xout,yout,'o');
else
    % from integral over y(max +- window/2)
    line('XData', x,'YData', y, ...
        'LineWidth', 0.5, 'LineStyle', '-', 'Color', [0 0 .8], ...
        'Parent', hAxes(1));
    
    xout = x(bestindex);
    yout = zeros(length(bestindex),1);
    for ii = 1:length(bestindex)
        lower = iof(x,x(bestindex(ii))-window/2);
        upper = iof(x,x(bestindex(ii))+window/2);
        yout(ii) = sum(y(lower:upper));
        rectangle('Position',[x(lower),0,x(upper)-x(lower),yout(ii)],...
            'EdgeColor','none', ...
            'FaceColor', [.8 .0 .0], ...
            'Parent', hAxes(2));
    end
    xlim(hAxes(2),'auto');
    ylim(hAxes(2),'auto');
end