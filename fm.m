function [xout, yout] = fm(x,y,varargin)

p = inputParser;
p.addRequired('x', @(x)validateattributes(x,{'numeric'},{'vector'}));
p.addRequired('y', @(x)validateattributes(x,{'numeric'},{'vector'}));
% DEFAULTS
% carrier frequency relative to the x-axis "frequency"
% a value of 0.5 keeps the number of points close to the original data
% higher values result in more points
p.addOptional('carrier', 0.5, @(x)validateattributes(x,{'numeric'},{'scalar'}));
% modulation level as a fraction of the carrier
% 0 <= level <= 1
% 0 yields no effect, 1 yields zero points for the minimal y values
% 0.9-0.95 is probably good
p.addOptional('level', 0.95, @(x)validateattributes(x,{'numeric'},{'scalar','<=',1,'>=',0}));
% sample rate for interpolation, relative to carrier + modulation
% has to obey Nyquist, so must be >2
% higher values might be needed for a smooth distribution
p.addOptional('srt', 2.2, @(x)validateattributes(x,{'numeric'},{'scalar','>',2}));
% point density proportional to value ('normal') or absolute value ('absolute')?
p.addParameter('mode', 'normal', @(x)ischar(validatestring(x,{'normal', 'absolute'})));

p.FunctionName = 'fm';
p.parse(x,y,varargin{:});

% calculate actual parameter values
carrier = p.Results.carrier * length(x) / 2;
level = p.Results.level * carrier;
srt = p.Results.srt * (carrier + level);
kf = level * 2 * pi;

ts = 1/srt;   % sampling interval
t  = -1:ts:1; % new x-axis

% interpolate data with signal rate
xout = linspace(x(1),x(end),length(t));
yout = interp1(x,y,xout);

% generate signal function (renormalise to [-1 1]
if strcmp(p.Results.mode, 'absolute')
    m_sig = abs(yout);
    m_sig = 2 * (m_sig - min(m_sig)) / (max(m_sig) - min(m_sig)) - 1;
else
    m_sig = 2 * (yout - min(yout)) / (max(yout) - min(yout)) - 1;
end

% frequency modulation
m_intg = kf * ts * cumsum(m_sig);
s_fm = cos(2*pi*carrier*t + m_intg);




% find zero crossings
n = length(t);
t1 = s_fm(1:n-1);
t2 = s_fm(2:n);
tt = t1.*t2;
indx = find(tt < 0);

% filter zero crossings from resampled data
xout = xout(indx); yout = yout(indx);

% plot it
limits = [-1 1 -1.2 1.2];

figure(1)
subplot(311); plot(t,m_sig);
axis(limits) % set x-axis and y-axis limits 
xlabel('t (s)'); ylabel('m (t)')
title('signal')

subplot(312); plot(t,s_fm);
axis(limits) % set x-axis and y-axis limits 
title('frequency modulated signal')
xlabel('t (s)'); ylabel('s_{FM}(t)')

subplot(313); plot(t(indx),m_sig(indx),'.');
axis(limits) % set x-axis and y-axis limits 
title('filtered signal')
xlabel('t (s)'); ylabel('m(t)')

end