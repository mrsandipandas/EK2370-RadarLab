%%%%%% This script obtains the SAR image from a FMCW radar measurement
%%%%%% moving alogn a linear rail
close all
clear all
clc

[y,fs] = audioread('Data/SAR.m4a');             % Load audio file

c = 2.997e8;                                    % Speed of light in air (m/s)
f_start = 2.4e9;                                % Start Frequency (Hz)
f_stop = 2.5e9;                                 % Stop Frequency (Hz)
BW = f_stop-f_start;                            % Bandwidth
f_c = f_start + BW/2;                           % Center frequency
Tsample = 20e-3;                                % Pulse width(s)
Trp = 250e-3;                                   % Range profile time duration
N = fs*Tsample;                                 % Number of samples in each pulse
Nrp = fs*Trp;

data = y(:,1);                                  % Backscattered data from radar
trig = y(:,2);                                  % Trigger signal from modulator

rpstart = abs(trig) > mean(abs(trig));                        % Set pulse start
count = 0;

for i = Nrp+1:(size(rpstart,1)-Nrp)
    if rpstart(i) == 1 && sum(rpstart(i-Nrp:i-1)) == 0
        count = count +1;
        RP(count,:) = data(i:i+Nrp-1);
        RPtrig(count,:) = trig(i:i+Nrp-1);
    end
end

thresh = 0.08;

for j= 1:size(RP,1)
    clear SIF;
    SIF = zeros(N,1);
    start = (RPtrig(j,:) > thresh);
    count=0;
    for i = 12:(size(start,2)-2*N)
        [Y,I] = max(start(1,i:i+2*N));
        if mean(start(i-10:i-2)) == 0 && I == 1
            count = count + 1;
            SIF = RP(j,i:i+N-1)'+SIF;
        end
    end
    SI = SIF/count;
    FF = ifft(SI);
    clear SI;
    sif(j,:) = fft(FF(size(FF,1)/2+1:size(FF,1)));
end

cr = BW/Tsample;                        % Chirp rate
Rs = 0;                                 % Down scene distance to range center
delta_x = c/f_c/2;                       % Radar spacing bewtween each range profile aquisition
L = delta_x*(size(sif,1));
Xa = linspace(-L/2,L/2,(L/delta_x));
time = linspace(0, Tsample, size(sif,2));
Kr = linspace(((4*pi/c)*(f_start)),((4*pi/c)*(f_stop)),(size(time,2)));

for j = 1:size(sif,1)
sif(j,:) = sif(j,:) - mean(sif,1);
end

clear N;
N = size(sif,2);
for i = 1:N
    H(i) = 0.5 + 0.5*cos(2*pi*(i-N/2)/N);
end

for i = 1:size(sif,1)
    sif_h(i,:) = sif(i,:).*H;
end
sif = sif_h;

fig_count=1;
figure(fig_count);
S_image = angle(sif);
imagesc(Kr, Xa, S_image);
colormap('default');
xlabel('K_r (rad/m)');
ylabel('SAR Position, Xa (m)');
colorbar;
fig_count = fig_count + 1;

zpad = 2048;
s_zeros = zeros(zpad, size(sif,2));
for i = 1:size(sif,2)
index =round((zpad-size(sif,1))/2);
s_zeros(index+1:(index + size(sif,1)),i) = sif(:,i);
end
sif = s_zeros;

S = fftshift(fft(sif,[],1),1);
Kx = linspace((-pi/delta_x),(pi/delta_x),(size(S,1)));

figure(fig_count);
S_image = 20*log10(abs(S));
imagesc(Kr, Kx, S_image, [max(max(S_image))-40,max(max(S_image))]);
colormap('default');
xlabel('K_r (rad/m)');
ylabel('K_x (rad/m)');
colorbar;
fig_count =fig_count+1;

figure(fig_count);
S_image = angle(S);
imagesc(Kr, Kx, S_image);
colormap('default');
xlabel('K_r (rad/m)');
ylabel('K_x (rad/m)');
colorbar;
fig_count = fig_count + 1;

if Rs==0
S_matched=S;
end

kstart=floor(min(Kr)-20);
kstop=ceil(max(Kr));
Ky_e=linspace(kstart,kstop,1024);

count = 0;
for ii = 1:zpad
count = count + 1;
Ky(count,:) = sqrt(Kr.^2 - Kx(ii)^2);
S_st(count,:) = (interp1(Ky(count,:), S_matched(ii,:), Ky_e));
end

S_st(find(isnan(S_st))) = 1E-30;

figure(fig_count);
S_image = angle(S_st);
imagesc(Ky_e, Kx, S_image);
colormap('default');
xlabel('K_y (rad/m)');
ylabel('K_x (rad/m)');
colorbar;
fig_count = fig_count + 1;

v = ifft2(S_st,(size(S_st,1)*4),(size(S_st,2)*4));

bw = c*(kstop-kstart)/(4*pi);
max_range = (c*size(S_st,2)/(2*bw));
figure(fig_count);
S_image = v;
S_image = fliplr(rot90(S_image));
cr1 = -20; % depends on the Kx of the Stolt Interpolation
cr2 = 20; % depends on the Kx of the Stolt Interpolation
dr1 = 1;
dr2 = 100;

dr_index1 = round((dr1/max_range)*size(S_image,1));
dr_index2 = round((dr2/max_range)*size(S_image,1));
cr_index1 = round(( (cr1+zpad*delta_x/(2*1)) /(zpad*delta_x/1))*size(S_image,2));
cr_index2 = round(((cr2+zpad*delta_x/(2*1)) /(zpad*delta_x/1))*size(S_image,2));
trunc_image = S_image(dr_index1:dr_index2,cr_index1:cr_index2);
downrange = linspace(-1*dr1,-1*dr2, size(trunc_image,1));
crossrange = linspace(cr1, cr2, size(trunc_image, 2));
for ii = 1:size(trunc_image,2)
trunc_image(:,ii) = (trunc_image(:,ii)').*(abs(downrange*1)).^(3/2);
end
trunc_image = 20 * log10(abs(trunc_image));
imagesc(crossrange, downrange, trunc_image, [max(max(trunc_image))-30,max(max(trunc_image))-0]);
colormap('default');
axis equal;
colorbar;