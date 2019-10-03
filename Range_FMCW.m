%%%%% This script obtains the range of a target from an audio file

close all;
clear all;
clc
tic
[y,fs] = audioread('Data/FMCW_test_05.wma');           % Load audio file

c = 2.997e8;                                    % Speed of light in air (m/s)
f_start = 2.4e9;                                % Start Frequency (Hz)
f_stop = 2.5e9;                                 % Stop Frequency (Hz)
BW = f_stop-f_start;                            % Bandwidth
Tsample = 20e-3;                                % Pulse width(s)
N = fs*Tsample;                                 % Number of samples in each pulse
rr = c/(2*BW);                                  % Range resolution (m)
range_max = rr*N/2;                             % Maximum range due to sampling frequency

data = y(:,1);                                  % Backscattered data from radar
trig = y(:,2);                                  % Trigger signal from modulator

thresh = 0;                                     % Treshold level for the trigger signal
start = (trig > thresh);                        % Set pulse start

count = 0;
for i = 100:(size(start)-N)
    if start(i) == 1 && mean(start(i-11:i-1)) == 0
        count = count +1;
        s(count,:) = data(i:i+N-1);
        time(count) = i/fs;
    end
end
ave = mean(s,1);
for j = 1:size(s,1)
        s(j,:) = s(j,:)-ave;
end
sfft = fft(s,4*N,2);                            % Calculate DFT of y using zero padding
sfft = 20*log10(abs(sfft));
sfft = sfft(:,1:end/2);                         % Remove upper part of DFT
sfft = sfft-max(max(sfft));                     % Normalize DFT results
range = linspace(0,range_max,4*N);              % Create frequency array

imagesc(range,time,sfft);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-50 0]);
toc