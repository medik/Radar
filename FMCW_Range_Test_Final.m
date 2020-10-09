clear

[y, fs] = audioread('FMCW_single_target.wav');
% pause(2);
% fs=44100;
% recorder1 = audiorecorder(fs,16,2,6);
% recordblocking(recorder1, 60);
% y = getaudiodata(recorder1);

sig = -y(:,1)';
sync = -y(:,2)' >= 0.1;

figure(1)
plot(sig)
title('Data plot')
xlabel('Samples')
ylabel('Signal magnitude')

figure(2)
plot(-y(:,2)')
title('Sync plot (before digitization)')
xlabel('Samples')
ylabel('Signal magnitude')

deltaT = 1/(fs);
N=(20e-3)/(deltaT);

time_vec=[];
matrix=[];

sync_start = [];
matrix_down = [];

n=1;
while n <= length(sync)
    if sync(n)
        sync_start = [sync_start n];
        time_vec=[time_vec n/fs];
        
        if n+N-1 <= length(sig)
            temp = sig(n:n+N-1);
        else
            temp = [sig(n:end)];
            temp = [temp zeros(1,N-length(temp))];
        end
        
        matrix = [matrix; temp;];
        n = n+N-1;
        if n <= length(sig) && 1
             while sync(n) > 0
                 n=n+1;
             end
        end
        
        n_temp = n;
         
        if n+N-1 <= length(sig)
            temp = sig(n:n+N-1);
        else
            temp = [sig(n:end)];
            temp = [temp zeros(1,N-length(temp))];
        end
        
        matrix_down = [matrix_down; temp;];
        n = n_temp;
    end
    n=n+1;
end

%%
% MS Clutter rejection
clutter_rej=1;
if clutter_rej
    [mx_N, mx_M] = size(matrix);

    for i=[1:mx_M]
        col=matrix(:,i);
        col_mean = mean(col);
        matrix(:,i)=col-col_mean;
    end
end

temp_mx = [];
[b,d] = size(matrix);
for i = [1:b]
    row = matrix(i,:);
    temp_mx = [temp_mx; ifft(row,4*N)];
end
matrix = temp_mx;

%%
% MS Clutter rejection
clutter_rej=1;
if clutter_rej
    [mx_N, mx_M] = size(matrix_down);

    for i=[1:mx_M]
        col=matrix_down(:,i);
        col_mean = mean(col);
        matrix_down(:,i)=col-col_mean;
    end
end

temp_mx = [];
[b,d] = size(matrix_down);
for i = [1:b]
    row = matrix_down(i,:);
    temp_mx = [temp_mx; ifft(row,4*N)];
end
matrix_down = temp_mx;

% 
% %2-pulse MTI
% k = 2;
% [P,~] = size(matrix);
% temp_MTI = [];
% while k <= P
%     temp = matrix(k,:) - matrix(k-1,:);
%     temp_MTI = [temp_MTI; temp];
%     k = k + 1;
% end
% matrix = temp_MTI;

%3-pulse MTI
k = 3;
[P,~] = size(matrix);
temp_MTI = [];
temp_MTI(1,:) = matrix(1,:);
temp_MTI(2,:) = matrix(2,:);
while k <= P
    temp = matrix(k,:) - 2*matrix(k-1,:) + matrix(k-2,:);
    temp_MTI = [temp_MTI; temp];
    k = k + 1;
end
matrix = temp_MTI;

% 
% %2-pulse MTI
% k = 2;
% [P,~] = size(matrix);
% temp_MTI = [];
% while k <= P
%     temp = matrix(k,:) - matrix(k-1,:);
%     temp_MTI = [temp_MTI; temp];
%     k = k + 1;
% end
% matrix = temp_MTI;

%3-pulse MTI
k = 3;
[P,~] = size(matrix_down);
temp_MTI = [];
temp_MTI(1,:) = matrix_down(1,:);
temp_MTI(2,:) = matrix_down(2,:);
while k <= P
    temp = matrix_down(k,:) - 2*matrix_down(k-1,:) + matrix_down(k-2,:);
    temp_MTI = [temp_MTI; temp];
    k = k + 1;
end
matrix_down = temp_MTI;

c=299792458;

%% Finding beat frequency for each up and down chirp

[P,~] = size(matrix);
vel_vec=[];
for i=[1:P]
    [~, i_up] = max(matrix(i,1:100));
    [~, i_down] = max(matrix_down(i,1:100));
    
    % highest frequency is N*fs/2
    freq_up = i_up*(fs/(2*4*N));
    freq_down = i_down*(fs/(2*4*N));
    delta_f = (freq_up-freq_down)/2;

    vel_vec = [vel_vec delta_f*c/( (2.45e9) )];
    %vel_vec = [vel_vec delta_f]
end


% FFT and normalize
matrix_fft = abs(matrix);
matrix_max_db = mag2db(max(max(matrix_fft)));
matrix_fft_db = mag2db(matrix_fft)-matrix_max_db;



%%
figure(3)
f_max = fs/2;
deltaf = 2.495e9-2.408e9;

deltaR = c/(2*deltaf);
Rmax = 4*N*deltaR/2;
M=100;
R_vec = [1:N/2]*deltaR/4;
imagesc(R_vec(1:M), time_vec, matrix_fft_db(:,1:M), [-50 0]);
xlabel('Range [m]')
ylabel('Time [s]')
colorbar
title('Range vs time spectogram')


%%

figure(4)
clf
hold on
for n=sync_start
    if n+N-1 <= length(sig)
        plot([n:n+N-1],sig(n:n+N-1), 'b')
        plot([n:n+N-1],sync(n:n+N-1), 'r')
    end
end
legend('Parsed up chirp data','Sync')
hold off
title('Captured parsed up-chirp data')
xlabel('Samples')
ylabel('Signal magnitude')

%%
x = 1;
max_index = [];
target_range = [];
matrix=matrix_fft_db;
while x <= size(matrix,1)
    [~,index] = max(matrix(x,1:M));
    target_range = [target_range; R_vec(index)];
    x = x + 1;
end
figure(5)
plot(time_vec, target_range);
title('Range vs time')
xlabel('Time [s]')
ylabel('Range [m]')
grid MINOR


%%
figure(6)
del_x = diff(target_range');
del_t = diff(time_vec);

vel = -1*movmean(del_x./del_t, 70);

plot(time_vec(1:end), [0 vel])
title('Velocity vs time using two point finite difference method')
xlabel('Time [s]')
ylabel('Velocity [m/s]')

