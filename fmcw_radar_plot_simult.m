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

time_step=0;

matrix = zeros(length(time_vec), alpha*N);
mean_vec = [];

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

    % Retrieving signal
    sig = -y(:,1);

    % Digitize the sync
    sync = -y(:,2) >= 0.1; 

    % Finding the upchirp (=1)
    n=1;
    while n <= length(sync)
      if sync(n)
	time_step = time_step+1;
	
	% Check if the length is exactly Tchirp, add padding otherwise
        if n+N-1 <= length(sig)
          temp = sig(n:n+N-1);
        else
          temp = [sig(n:end)];
          temp = [temp; zeros(N-length(temp),1)];
        end

	% Take an FFT of the time domain signal (with length Tchirp)
        temp_fft = fft(temp, alpha*N);

	% Making sure that it adds more rows when the matrix is filled
        if time_step <= length(time_vec)
          matrix(time_step, :) = temp_fft;
        else
          matrix = [matrix; temp_fft']; % temp_fft is a column vector, transpose to make it a row
          time_vec = [time_vec; time_vec(end)+Tchirp];
        end

	%% Update mean vector
	if time_step == 1
	  mean_vec=temp_fft;
	else
	  mean_vec=(mean_vec+temp_fft)/time_step;
	end

	% We have already processed Tchirp, move ahead of time. 
        n = n+N-1;
      end
        n=n+1;
    end
    toc

    [P,Q] = size(matrix);
    matrix_clut_rej=zeros(P,Q);
    
    tic
    disp('MS Clutter rej')
    % MS clutter rejection
    clutter_rej=1;
    if clutter_rej
        [mx_N, mx_M] = size(matrix);

        for i=[1:mx_M]
            col=matrix(:,i);
            col_mean = mean_vec(i);
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
    tocS    
end


