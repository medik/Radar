clear
%load velocity_test_file.mat
[y, fs] = audioread('sdr_multi_target.wav');

% pause(2);
% fs=44100;
% recorder1 = audiorecorder(fs,16,2,6);
% recordblocking(recorder1, 60);
% y = getaudiodata(recorder1);

ch1 = y(:,1);

 figure(1)
 plot(ch1);
% title('Captured data plot, CW radar');
% xlabel('Samples')
% ylabel('Amplitude')
% grid MINOR

%%
ch1_mean = mean(ch1);
ch1 = ch1 - abs(ch1_mean);

% sync
%ch2 = -y(:,2);

alpha=4;

deltaT = 1/fs;
N=4410;
%N=2000;
Tdwell = N*deltaT;

f_max = fs/2;
deltaf = f_max/N;
c=299792458;
f=2.43e9; % center frequency
f=5.8e9
freq_vec = [0:N-1]*deltaf/alpha;
vel_vec = freq_vec * (c/(2*f));

time_vec = [];
matrix = [];
matrix_2 = [];
tnow=0;
max_db=0;
%%
while ~isempty(ch1)
    if length(ch1) >= N
        temp = ch1(1:N)';
        ch1 = ch1(N+1:end);
    else
        temp = ch1(1:end)';
        temp = [temp zeros(1,N-length(temp))];
        ch1 = [];
    end
    
    temp_fft = fft(temp, alpha*N);
    temp_max = mag2db(max(abs(temp_fft)));
    
    if max_db < temp_max
        max_db=temp_max;
    end
    
    temp_fft_db = mag2db(abs(temp_fft));
    

    matrix = [matrix; temp_fft_db];
    matrix_2 = [matrix_2; temp;];
    time_vec = [time_vec tnow];
    tnow = tnow + Tdwell;
end
%%

temp_mx = [];
[b,d] = size(matrix_2);
for i = [1:b]
    row = matrix_2(i,:);
    temp_mx = [temp_mx; mag2db(abs(fft(row,alpha*N)))];
end


%%

c=299792458;
f=2.43e9; % center frequency
f=5.8e9;
freq_vec = [0:N-1]*deltaf/alpha;
vel_vec = freq_vec * (c/(2*f));

max_db = max(max(temp_mx));
matrix=temp_mx-max_db;
M=100;
figure(3)
imagesc(vel_vec(1:M), time_vec, matrix(:,1:M), [-45 0]);
colorbar
xlabel('Velocity [m/s]')
ylabel('Time [s]')
title('Velocity vs time Spectogram (lower if)')


% %% 
% x = 1;
% max_index = [];
% target_velocity = [];
% while x <= size(matrix,1)
%     [~,index] = max(matrix(x,1:M));
%     target_velocity = [target_velocity; vel_vec(index)];
%     x = x + 1;
% end
% figure(3)
% plot(time_vec, target_velocity);
% title('Velocity vs time, multi target')
% xlabel('Time [s]')
% ylabel('Velocity [m/s]')
% grid MINOR
