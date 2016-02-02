function fm_demo

% Demo program for point density 

% x-axis
x = linspace(-0.05,0.05,201);

% yaxis
% Use triangl function to generate message signal
Ta = 0.01;
Tl = 0.01;
Tr = 0.01;
y=triangl((x + Tl)/Ta)-triangl((x - Tr)/Ta);

% or use a gaussian
offset = 0;
width = 0.0001;
y = exp(-(x - offset).^2/width);

% fm it
% DEFAULTS
% carrier frequency relative to the x-axis "frequency", default 0.5
% a value of 0.5 keeps the number of points close to the original data
% higher values result in more points
carrier = 0.3;
% modulation level as a fraction of the carrier, default 0.95
% 0 <= level <= 1
% 0 yields no effect, 1 yields zero points for the minimal y values
% 0.9-0.95 is probably good
level = 0.95;
% sample rate, relative to carrier + modulation, default 2.2
% has to obey Nyquist, so must be >2
% higher values might be needed for a smooth distribution
srt = 5;
% mode
% 'normal' (default) or 'absolute': point density proportional to signal or
% absolute signal
[x2,y2] = fm(x,y,carrier,level,srt,'mode','normal');

% plot it
figure(2)
plot(x,y,'--',x2,y2,'.','MarkerSize',10)



%%% Defining triangl function  used in above code
% triangl(t)=1-|t| , if |t|<1
% triangl(t)=0 , if |t|>1
function y = triangl(t)
y=(1-abs(t)).*(t>=-1).*(t<1); % i.e. setting y to 1 -|t|  if  |t|<1 and to 0 if not
end

end