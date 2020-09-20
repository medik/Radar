fs=44100;
Td=0.1; % dwell time

recorder1 = audiorecorder(fs,16,2,0);

N=Td*fs;
alpha=4;
f_max = fs/2;
deltaf = f_max/N;
c=299792458;
f=2.43e9; % center frequency
freq_vec = [0:alpha*N-1]*deltaf;
vel_vec = freq_vec * (c/(2*f));

Time_Max = 60;
time_vec = [0:Td:Time_Max];



time_now=0;

matrix = zeros(length(time_vec), alpha*N);

while 1
    recordblocking(recorder1, Td);
    y = getaudiodata(recorder1);
    
    ch1 = -y(:,1);
    ch1_mean = mean(ch1);
    ch1 = ch1 - abs(ch1_mean);
    
    % FFT
    
    temp_fft = fft(ch1, alpha*N);
    temp_max = mag2db(max(abs(temp_fft)));
    temp_fft = mag2db(abs(temp_fft))-temp_max;
    
    matrix(1+mod(time_now,length(time_vec)),:) = temp_fft;
    
    imagesc(vel_vec, time_vec, matrix, [-45 0]);
    time_now=time_now+1;
end


