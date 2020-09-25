clear all;
clc;

% Defining constants
[y, fs] = audioread('SAR_Test_File.m4a');
c = 299792458;
fc = 2.43e9; % center frequency
f_start = 2.408e9;
f_stop = 2.495e9;
BW = f_stop - f_start;
Ts = 0.02;
Trp = 0.250;
Ns = fs*Ts;
Nrp = fs*Trp;

% Labeling data
data = -y(:,1);
trig = -y(:,2);

% Identify range profile acquisition and parse data and sync signal
rp_start = abs(trig) > mean(abs(trig));
count = 0;
RP = [];
RP_trig = [];
for i = Nrp+1:(size(rp_start,1) - Nrp)
    if rp_start(i) == 1 && sum(rp_start(i-Nrp:i-1)) == 0
        count = count + 1;
        RP(count,:) = data(i:i+Nrp-1);
        RP_trig(count,:) = trig(i:i+Nrp-1);
    end
end

% Parse data according to rising edge trigger
thresh = 0.08;
for j = 1:size(RP,1)
    clear SIF;
    SIF = zeros(Ns,1);
    start = (RP_trig(j,:) > thresh);
    count = 0;
    for i = 12:(size(start,2) - 2*Ns)
        [Y, I] = max(start(1,i:i+2*Ns));
        if mean(start(i-10:i-2)) == 0 && I == 1
            count = count + 1;
            SIF = RP(j,i:i+Ns-1)' + SIF;
        end
    end
    SI = SIF/count;
    FF = ifft(SI);
    clear SI;
    sif(j,:) = fft(FF(size(FF,1)/2+1:size(FF,1)));
end

% Replace NaN with 1e-30 in sif
[row, column] = size(sif);
for m = 1:row
    for n = 1:column
        if isnan(sif(m,n))
            sif(m,n) = 1e-30;
        end
    end
end

% Define more constants
cr = BW/Ts;
Rs = 0;
lambda = c/fc;
delta_x = lambda/2;
L = delta_x*(size(sif,1));
Xa = linspace(-L/2, L/2, (L/delta_x));
time = linspace(0, Tsample, size(sif,2));
Kr = linspace(((4*pi/c)*(fc - BW/2)), ((4*pi/c)*(fc + BW/2)), (size(time,2)));

for j = 1:size(sif,1)
    sif(j,:) = sif(j,:) - mean(sif,1);
end

% Mean subtraction for ground clutter reduction
%Clear N;
N = size(sif,2);
H = [];
for i= 1:N
    H(i) = 0.5 + 0.5*cos(2*pi*(i-N/2)/N);
end

% Hanning window
sif_h = [];
for i = 1:size(sif,1)
    sif_h(i,:) = sif(i,:).*H;
end
sif = sif_h;

% Plot the phase of the SAR data matrix before along track FFT
fig_count = 1;
figure(fig_count);
S_image = angle(sif);
imagesc(Kr, Xa, S_image);
colormap('default');
xlabel('K_r(rad/m)');
ylabel('SAR Position, Xa(m)');
colorbar;
figcount = figcount + 1;

% Wow ok!