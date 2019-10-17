clc
clear all
[Y, fs] = audioread('Data/Umer_Running.m4a');
c = 3e8; % Speed of light
fc = 2.424e9; % Carrier frequency
T = 100e-3; % Pulse length
N = fs*T; % Number of samples in each pulse

% Reformating data structure for taking into account the first channel only
Y = Y(:,1);
len = size(Y,1);
row = floor(len/N);

caxis_range = -35;
caxis([caxis_range 0]);
colorbar;
set(gca,'XLim',[0 40]);
set(gca,'YLim',[0 row*T]);
xlabel('Velocity [m/sec]'); ylabel('Time [sec]');

f = zeros(1, 2*N);
delta_f = linspace(0,fs/2,2*N);
vel = (delta_f*c)/(2*fc);
time = linspace(1, T*row, row);
hold on

% Visualization of Velocity information
hImage = imagesc(vel,time,f);

cdata = caxis_range*ones(row, 2*N);
itr = 1;
for i = 1:N:len
    if i+N-1 < len
        data = (Y(i:i+N-1,:))';
    end
    % Take mean and reduce ground clutter
    data = data - mean(data);
    % FFT across the rows for multiple pulses
    % 4*N padding for smoother visualization
    f = abs(fft(data,4*N,2));
    f = 20*log10(f);
    
    % Only considering lower half of the frequency domain
    f = f(:,1:size(f,2)/2);
    f = f - max(f);
    delta_f = linspace(0,fs/2,size(f,2));
    cdata(itr,:) = f;
    itr = itr + 1;
    % Visualization of Velocity information
    set(hImage, 'CData', cdata);
    pause(0.1);
end

