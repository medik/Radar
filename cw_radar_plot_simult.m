fs=44100;
%Td=0.5; % dwell time

recorder1 = audiorecorder(fs,16,2,6);

deltaT = 1/fs;
N=4410;
%N=2000;
Td = N*deltaT;

%N=Td*fs;
alpha=4;
f_max = fs/2;
deltaf = f_max/N;
c=299792458;
f=2.43e9; % center frequency
freq_vec = [0:N-1]*deltaf/alpha;
vel_vec = freq_vec * (c/(2*f));

Time_Max = 60;
time_vec = [0:Td:Time_Max];

time_now=0;

matrix = zeros(length(time_vec), alpha*N);


max_db=0
while 1
    recordblocking(recorder1, Td);
    y = getaudiodata(recorder1);
    
    ch1 = -y(:,1);
    ch1_mean = mean(ch1);
    ch1 = ch1 - abs(ch1_mean);
    
    %% FFT
    
    temp_fft = fft(ch1, alpha*N);
    temp_max = mag2db(max(abs(temp_fft)));
    
    if max_db < temp_max
        max_db=temp_max;
    end
    
    matrix(1+mod(time_now,length(time_vec)),:) = temp_fft;
    
    figure(8)
    M=100;
    imagesc(vel_vec(1:M), time_vec, matrix(:,1:M)-max_db, [-45 0]);
    colorbar
    xlabel('Velocity [m/s]')
    ylabel('Time [s]')
    title('3N padding')
    
    time_now=time_now+1;
end


