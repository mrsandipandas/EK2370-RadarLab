tic
[Y, fs] = audioread('Data/sdr_wo_pa_i.wav');
c = 3e8; % Speed of light
fc = 2.424e9; % Carrier frequency
T = 100e-3; % Pulse length
N = fs*T; % Number of samples in each pulse

% Reformating data structure for taking into account the first channel only
Y = Y(:,1);

% Clipping part of the data to consider full pulse lengths
row = size(Y,1);
row = row - mod(row,N);
Y = Y(1:row,:);
Y = reshape(Y,N,[])';

% Take mean and reduce ground clutter
Y = Y - mean(Y,1);

% FFT across the rows for multiple pulses
% 4*N padding for smoother visualization
f = abs(fft(Y,4*N,2));
f = 20*log10(f);

% Only considering lower half of the frequency domain
f = f(:,1:size(f,2)/2);
f = f - max(max(f));
delta_f = linspace(0,fs/2,size(f,2));
vel = (delta_f*c)/(2*fc);
time = linspace(1, T*size(f,1), size(f,1));

% Visualization of Velocity information
imagesc(vel,time,f);
caxis([-35 0]);
colorbar;
set(gca,'XLim',[0 40]);
xlabel('Velocity [m/sec]'); ylabel('Time [sec]');
toc












