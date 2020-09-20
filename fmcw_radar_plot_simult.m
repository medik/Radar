clear
disp('Radar on')
fs=44100;
Td=50e-3; % dwell time

recorder1 = audiorecorder(fs,16,2,6);

Tchirp=20e-3;
Nchirp=Tchirp*fs;

% Number of column of the matrix is equal to the chirp time for FMCW
N=Tchirp*fs;
alpha=4;
f_max = fs/2;
deltaf = f_max/N;
c=299792458;
%f=2.43e9; % center frequency
%freq_vec = [0:alpha*N-1]*deltaf;

deltaf = 2.495e9-2.408e9;
deltaR = c/(2*deltaf);
%Rmax = alpha*N*deltaR/2;
R_vec = [1:N/2]*deltaR;

Time_Max = 60;
time_vec = [0:Tchirp:Time_Max];

time_now=0;

matrix = zeros(length(time_vec), alpha*N);
disp('Recording')
record(recorder1);
pause(1);
while 1
    tic
    disp('Retrieving data')
    pause(Td);
    y = getaudiodata(recorder1);
    toc
    tic
    disp('Processing data')
    sig = -y(:,1);
    sync = -y(:,2);

    % Finding sync 

    n=1;
    while n <= length(sync)
        if sync(n)
            if n+N-1 <= length(sig)
                temp = sig(n:n+N-1);
            else
                temp = [sig(n:end)];
                temp = [temp zeros(1,N-length(temp))];
            end

            % FFT
            temp_fft = fft(temp, alpha*N);
            if time_now <= length(time_vec)
                matrix(time_now, :) = temp_fft;
                time_now = time_now+1;
            else
                matrix = [matrix; temp_fft];
                time_vec = [time_vec; time_vec(end)+Tchirp];
            end    
            n = n+N-1;
        end
        n=n+1;
    end
    toc
    
    tic
    disp('MS Clutter rej')
    % MS clutter rejection
    clutter_rej=1;
    if clutter_rej
        [mx_N, mx_M] = size(matrix);

        for i=[1:mx_M]
            col=matrix(:,i);
            col_mean = mean(col);
            matrix(:,i)=col-col_mean;
        end
    end
    toc

    tic
    disp('MTI')  
    % MTI
    mti=2;
    [P,Q] = size(matrix);
    temp_MTI = zeros(P,Q);
    if mti==2
        %2-pulse MTI
        k = 2;
        temp_MTI(1,:) = matrix(1,:);
        while k <= P
            temp_MTI(k,:) = matrix(k,:) - matrix(k-1,:);
            k = k + 1;
        end
    else
        %3-pulse MTI
        k = 3;
        temp_MTI(1,:) = matrix(1,:);
        temp_MTI(2,:) = matrix(2,:);
        while k <= P
            temp_MTI(k,:) = matrix(k,:) - 2*matrix(k-1,:) + matrix(k-2,:);
            k = k + 1;
        end
    end
    toc


    tic
    % Normalization
    disp('Normalization')
    matrix_fft = abs(matrix);
    matrix_max_db = mag2db(max(max(matrix_fft)));
    matrix_fft_db = mag2db(matrix_fft)-matrix_max_db;
    toc
    
    M=50;
    T_window=50/Tchirp;

    tic
    disp('Displaying')

    % Time window

    if T_window > length(time_vec)
      imagesc(R_vec(1:M), time_vec(end-T_window:end), matrix_fft_db(end-T_window:end,1:M), [-50 0]);
    else
      imagesc(R_vec(1:M), time_vec(1:end), matrix_fft_db(1:end,1:M), [-50 0]);
    end
    toc
    
    time_now=time_now+1;
    toc
end


