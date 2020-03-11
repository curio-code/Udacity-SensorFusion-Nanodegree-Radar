clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
R = 110;
v = -20;

%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.


%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

%speed of light
c = 3e8;       

%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

%Range Resolution
d_res = 1;

%BandWidth
B = c/(2*d_res);

%Chirp Time
MaxRange = 200;
TChirp = 5.5*(2*MaxRange/c);

%Slope
slope = B/TChirp;

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*TChirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));



%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
        
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    
    r_t(i) = R + v * t(i);
    td(i) = 2*r_t(i)/c;
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*(fc*t(i)+0.5*slope*t(i)^2));
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i))+0.5*slope*(t(i)-td(i))^2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
    
end

%% RANGE MEASUREMENT
len = B*TChirp; % lenght of the signal

 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix=reshape(Mix,[Nr,Nd]);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
signal_fft = fft(Mix,Nr);

 % *%TODO* :
% Take the absolute value of FFT output
signal_fft = abs(signal_fft/len);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
signal_fft = signal_fft(1:len/2+1);


f = B*(0:(len/2))/len;
 % plot FFT output
figure ('Name','Range from First FFT')
subplot(2,1,1)
plot(f,signal_fft)

R = (c*TChirp*f)/(2*B);
 % *%TODO* :
 
 %plotting the range
figure('Name','Range from First FFT')
plot(R,signal_fft)
 
axis ([0 200 0 0.5]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
T_h = 12;
T_v = 6;


% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
G_h = 6;
G_v = 3;


% *%TODO* :
% offset the threshold by SNR value in dB
offset=5;



signal_cfar = zeros(Nr/2,Nd);

% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.
train_cells = (2*T_v+2*G_v+1)*(2*T_h+2*G_h+1) - (2*G_v+1)*(2*G_h+1);


for i   = T_h+G_h+1:Nr/2-T_h-G_h
    for j = T_v+G_v+1:Nd-T_v-G_v
        S1 = sum(db2pow(RDM(i-T_h-G_h:i+T_h+G_h, j-T_v-G_v:j+T_v+G_v)),[1 2]);
        
        S2 = sum(db2pow(RDM(i-G_h:i+G_h, j-G_v:j+G_v)),[1 2]);
        
        train_total = S1-S2;
        
        
        threshold = train_total/train_cells;
        threshold = pow2db(threshold) + offset;
        threshold = db2pow(threshold);
        
        signal = db2pow(RDM(i, j));
        
        if (signal <= threshold)
            signal = 0;
        else 
            signal = 1;
        end
        
        signal_cfar(i,j) = signal;
        
    end
end
figure,surf(doppler_axis,range_axis,signal_cfar);
colorbar;


 
 