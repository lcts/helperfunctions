function [ data ] = digitize(varargin)

p = inputParser;
p.addRequired('imatrix', @(x)validateattributes(x,{'numeric'},{'3d'}));
p.addParamValue('mode', 'precise', @(x)ischar(validatestring(x,{'quick', 'precise'})));
p.addParamValue('yrange', 0.05, @(x)validateattributes(x,{'numeric'},{'scalar','<=',1,'>=',0}));
p.addParamValue('crange', 0.05, @(x)validateattributes(x,{'numeric'},{'scalar','<=',1,'>=',0}));
p.addParamValue('snaprange', 0.01, @(x)validateattributes(x,{'numeric'},{'scalar','<=',1,'>=',0}));

p.FunctionName = 'digitize';
p.parse(varargin{:});

A = double(p.Results.imatrix)/255;
yrange = round(p.Results.yrange*size(A,1));
snaprange = round(p.Results.snaprange*size(A,1));
crange = p.Results.crange;

% plot the image
close(findobj('type','figure','name','ImageFigure'))
hImageFigure = figure('name','ImageFigure');
hImageAxes = axes('Parent',hImageFigure);

image(A,'Parent',hImageAxes);
xlabel(hImageAxes, 'x coordinate / pixel');
ylabel(hImageAxes, 'y coordinate / pixel');

% set datacursormode and get cursor object
dcm_obj = datacursormode(hImageFigure);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on')

% select the data range
disp('Select the lower left corner of the plot area, then hit ''Enter''.')
pause
cinfo = getCursorInfo(dcm_obj);
llx = cinfo.Position(1);
lly = cinfo.Position(2);

disp('Select the upper right corner of the plot area, then hit ''Enter''.')
pause
cinfo = getCursorInfo(dcm_obj);
urx = cinfo.Position(1);
ury = cinfo.Position(2);

limits(1) = min(urx,llx);
limits(2) = min(ury,lly);
limits(3) = abs(urx-llx);
limits(4) = abs(ury-lly);

% and display it
rectangle('Position',limits, 'EdgeColor','w','FaceColor', 'none','Parent', hImageAxes)

% select background color
disp('Select a point not on the graph, then hit ''Enter''.')
pause
cinfo = getCursorInfo(dcm_obj);
bgx = cinfo.Position(1);
bgy = cinfo.Position(2);

% select the graph
disp('Select a point on the graph, then hit ''Enter''.')
pause
cinfo = getCursorInfo(dcm_obj);
refx = cinfo.Position(1);
refy = cinfo.Position(2);
if refx <= limits(1) || ...
        refx  >= limits(1) + limits(3) || ...
        refy  <= limits(2) || ...
        refy  >= limits(2) + limits(4)
    error('The selected point is outside the chosen graph boundary.')
end
tic
% convert picture to CIELab
A = rgb2lab(A);
% create color map, contrast gets more weight than hue
cmap = squeeze(sqrt(4 * A(:,:,1).^2 + A(:,:,2).^2 + A(:,:,3).^2));
% average the color in snaprange
bgval = mean(mean(cmap(bgy-snaprange:bgy+snaprange,bgx-snaprange:bgx+snaprange,:),2),1);

% create map of colors that are sufficiently far away from background
snapmap = abs(cmap(refy-snaprange:refy+snaprange,refx-snaprange:refx+snaprange,:) - bgval);

% look for point in snaprange that deviates most from bgcolor
[~,index] = max(snapmap(:));
[y,x] = ind2sub(size(snapmap),index);
refx = x + refx - snaprange;
refy = y + refy - snaprange;
refval = A(refy,refx,:);

% initialise graph search
x = refx;
y = refy;
newref = refval;
kk = 1;

% go from selected point to the right
for ii = x:urx
    % determine next colour range
    colmin = newref - crange;
    colmax = newref + crange;
    % scan yrange
    currentrun = 0;
    for jj = y-yrange:y+yrange
        % if outside the limits
        if jj > limits(2) + limits(4) || jj < limits(2)
            % move on
            continue
        end
        % if colour is within range
        if A(jj,ii,:) <= colmax & A(jj,ii,:) >= colmin
            % add a datapoint
            dright(1,kk) = ii;
            dright(2,kk) = jj;
            kk = kk + 1;
            % increase number of values found
            currentrun = currentrun + 1;
            % and add the value to the array of current matches
            matches(currentrun) = jj;
        end
    end
    % set the starting coordinate and colour for next x value
    if currentrun ~= 0
        y = min(matches) + ceil((max(matches) - min(matches))/2);
        newref = A(y,ii,:);
        clear('matches')
    end
end

% same thing again, but now we go from starting point to the left
x = refx;
y = refy;
newref = refval;
kk = 1;

for ii = x:-1:llx
    colmin = newref - crange;
    colmax = newref + crange;
    currentrun = 0;    
    for jj = y-yrange:y+yrange
        if jj > limits(2) + limits(4) || jj < limits(2)
            continue
        end
        if A(jj,ii,:) <= colmax & A(jj,ii,:) >= colmin
            % add a datapoint
            dleft(1,kk) = ii;
            dleft(2,kk) = jj;
            kk = kk + 1;
            % increase number of values found
            currentrun = currentrun + 1;
            % and add the value to the array of current matches
            matches(currentrun) = jj;
        end
    end
    % set the starting coordinate and colour for next x value
    if currentrun ~= 0
        y = min(matches) + ceil((max(matches) - min(matches))/2);
        newref = A(y,ii,:);
        clear('matches')
    end
end

% put the two parts together
dleft = flipdim(dleft,2);
data = [dleft dright];

% shift identical x-values around to make a strictly monotonic x axis
% calculate the mean spacing between x-values
meanspacing = abs(data(1,1) - data(1,end))/size(data,2);
% initialise counter
currentrun = 1;
% loop over the x-axis
for ii = 2:size(data,2)
    % if the current x-value is equal to the previous one
    if data(1,ii) == data(1,ii-1)
        % increment the counter
        currentrun = currentrun + 1;
    % if not
    else
        % if we just left a sequence of identical values
        if currentrun ~= 1
            % separate all values in the preceding sequence by 1/100 of the
            % mean spacing
            for jj = 1:currentrun
                data(1,ii-jj) = data(1,ii-jj) - (jj-1) * meanspacing/100;
            end
            % and reset the sequence counter
            currentrun = 1;
        end
    end
end
% we might have left the loop inside a sequence of identical numbers
if currentrun ~= 1
    % if so, shift those values
    for jj = 1:currentrun - 1;
        data(1,ii-jj) = data(1,ii-jj) - jj * meanspacing/100;
    end
end

% display the detected line
line('XData',data(1,:),'YData',data(2,:),'Color','w','LineWidth',1.5, 'Parent', hImageAxes)