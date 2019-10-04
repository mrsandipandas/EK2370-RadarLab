%%%%% This script obtains the range of a target from an audio file

close all
clear all
clc
tic
[Y,fs] = audioread('Data/FMCW_test_05.wma');           % Load audio file

c = 2.997e8;                                    % Speed of light in air (m/s)
f_start = 2.4e9;                                % Start Frequency (Hz)
f_stop = 2.5e9;                                 % Stop Frequency (Hz)
BW = f_stop-f_start;                            % Bandwidth
Tsample = 20e-3;                                % Pulse width(s)
N = fs*Tsample;                                 % Number of samples in each pulse
rr = c/(2*BW);                                  % Range resolution (m)
range_max = rr*N/2;                             % Maximum range due to sampling frequency

data = Y(:,1);                                  % Backscattered data from radar
trig = Y(:,2);                                  % Trigger signal from modulator
len = size(Y,1);
row = floor(len/N);

thresh = 0;                                     % Treshold level for the trigger signal
start = (trig > thresh);                        % Set pulse start

caxis([-30 0]);
set(gca,'XLim',[0 100]);
set(gca,'YLim',[0 row*Tsample]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
range = linspace(0,range_max,4*N);              % Create frequency array
t = linspace(0, row*Tsample, 893);
f = zeros(1, 2*N);

hold on

hImage = imagesc(range,t,f);
cdata = -30*ones(row, 2*N);
count = 0;
prev = 0;
i = 100;
sample_old = zeros(N,1);
while i < (size(start,1)-N)
    if start(i) == 1 && mean(start(i-11:i-1)) == 0
        count = count + 1;
        sample = (data(i:i+N-1))';
        sample = sample - mean(sample);
        del = sample - sample_old;
        sample_old = sample;
        if count > 1
            f = abs(fft(del,4*N,2));
            f = 20*log10(f);
            f = f(:,1:size(f,2)/2);
            f = f - max(f);
            t(count) = i/fs;
            cdata(count,:) = f;
            
            % Visualization of Velocity information
            set(hImage, 'CData', cdata);
            %set(hImage, 'YData', t(1:count));
            pause(0.001);
        end        
        i = i + N;
    end
    i = i + 1;
end


