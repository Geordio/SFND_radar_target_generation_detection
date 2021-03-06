clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%max range
R_max = 200;
%speed of light = 3e8
c = 3e8;
% fc = 77e9;
d_res = 1
max_v = 100;

%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
 
% set the range of the target vehicle
range = 100;

% set the velocity of target
v = 50;


%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

% d_res = c / (2*Bsweep);
B = c / (2*d_res);

% sweep time (chirp time)
Tchirp = 5.5 * R_max * 2 / c; % use 5.5 as per lesson

% Slope=Bandwidth/Tchirp
slope = B / Tchirp

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq
fprintf("B: %d, Tchirp: %d, slope: %d\n", B, Tchirp, slope);
                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


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
    % delay time can be calculated from range and speed of light
    r_t(i) = range + v*t(i);
    td(i) = 2*r_t(i)/c;
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*(fc*t(i) + slope*t(i)^2/2));
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i)) + slope*(t(i)-td(i))^2/2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);    
end

%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix, [Nr, Nd]);
disp(size(Mix));

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
signal_fft = fft(Mix,Nr);

 % *%TODO* :
% Take the absolute value of FFT output, after normalising
signal_fft = abs(signal_fft/Nr);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
signal_fft = signal_fft(1:Nr/2)

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

 % *%TODO* :
 % plot FFT output 
plot(signal_fft);

 
axis ([0 200 0 1]);

% return;


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
Tr = 12
Td = 8;

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 6;
Gd = 4;

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 12;

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1); 

% https://uk.mathworks.com/help/matlab/ref/zeros.html
rdm_size = size(RDM);
disp(rdm_size);
signal_cfar = zeros(rdm_size);

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


   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR
   
num_t_cells = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1) - (2*Gr+1)*(2*Gd+1);
% location of the cell under test within the sample window / roi
cut_col = Gd+Td+1;
cut_row = Gr+Tr+1;

fprintf("CUT location: %u,%u\n", cut_row, cut_col);

% slide a region of interest / kernel around the radar doppler map
for r_idx=1:(Nr/2 -2*(Gr+Tr))
    for d_idx=1:(Nd-2*(Gd+Td))
               
        % extract the block of elements containing the  training guard and
        % CUT
        % roi is the region of interest, a window that I slide across the
        % RDM in order to test each cell by sampling the noise
        % CUT is at teh centre of the ROI. ROI is mad up of the training
        % and guard cells in each dimension and CUT
        roi = RDM( r_idx:r_idx+2*(Gr+Tr),d_idx:d_idx+2*(Gd+Td));
     
        % get the value of the CUT before I start to overwrite the roi
        cut_val = roi(cut_row, cut_col);
        
        % set the guard and CUT to 0 so that they are not counted in sum to
        % get the noise level
        roi(Tr+1:Tr+1+2*Gr,Td+1:Td+1+(2*Gd)) = 0;
        
        % note that sum only return a vetcor of the column sums, so need
        % to do sum(sum()) in order to get all
        noise_level = sum(sum(db2pow(roi)));
        mean_of_t_cells = noise_level/num_t_cells;
        threshold = pow2db(mean_of_t_cells*offset);
        %         fprintf("r_idx: %u, d_idx: %u, sum: %u, mean: %d, threshold: %f\n", r_idx, d_idx,noise_level, mean_of_t_cells, threshold);
        

        if (cut_val<threshold)
            cut_val =0;
        else
            cut_val =1;
        end
        % insert the thresholded value into the results matrix
        signal_cfar(r_idx+cut_row, d_idx+cut_col) = cut_val;
    end
end



% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 
% by creating a results matrix signal_cfar the same size as the RDM, initialise to zero in the following manner:
% rdm_size = size(RDM);
% disp(rdm_size);
% signal_cfar = zeros(rdm_size);
% any elements that I do not overwrite in the loop remain at 0.


% disp(size(signal_cfar));

% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,signal_cfar);
colorbar;


 
 