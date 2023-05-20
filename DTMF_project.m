clearvars
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%'''Initialization'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = 'Enter phone number: ';
xstr=num2str(input(prompt,'s'));  
xarray=zeros(1,numel(xstr)); 
Fs = 8000;
N = 800;
L = 0;
M = 799;
A = 0;  %zero padding to the left
P = 0;  %zero padding to the left
Q = 0;  %zero padding to the right
W = 0;  %zero padding to the right
R = numel(xstr); %no. of symbols
m = 0;
phoneNum = [];
counter = 0;
counter1 = R+1;
n = R-1;  %no. of guards
len = 160; %length of guard
t2 = (0:((M*R)+(R-1)+(n*len)))/Fs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%'''PART2'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numel(xstr)
    xarray(i)=str2num(xstr(i));
    t = (L:M)/Fs;
    z = string(xarray(i));
    x_t = sym2TT(z);
    d = eval(x_t);
    counter = counter + 1;
    counter1 = counter1 - 1;
    A = 1;
    W = (R-counter)*960;
    P = (R-counter1)*960;
    d = [zeros(A,P) d zeros(A,W)];
    L = M+1+160;
    M = M+800+160;
    m = m + d;
    d = 0;
end
x_t = m;    %signal x(t)
figure;
plot(t2*1e3,x_t);
set(gca, 'Xlim', [0 ((100*R)+(20*n))]);
ylabel('Amplitude');
xlabel('Time (ms)');
title('x(t)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%'''PART3'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variance = 0.1;
W = sqrt(variance).*randn(1,size(x_t,2));   %WGN generation
y_t = x_t + W;                              %signal y(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%'''PART4'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load handel.mat

if max(x_t)>abs(min(x_t))                  %normalizing x(t) so no clipping occurs when making audio file
    normalized_x_t = (x_t)/max(x_t);
else
    normalized_x_t = (x_t)/abs(min(x_t));
end

filename = 'C:\Users\Abdelrhman Mohamed\Downloads\part4_sound_sample_2.wav'; %write your Own path

audiowrite(filename,normalized_x_t,Fs);
% sound(normalized_x_t,Fs); %sound of x(t) (uncomment thisline if you want to listen to x(t))

if max(y_t)>abs(min(y_t))                  %normalizing y(t) so no clipping occurs when making audio file
    normalized_y_t = (y_t)/max(y_t);
else
    normalized_y_t = (y_t)/abs(min(y_t));
end

filename = 'C:\Users\Abdelrhman Mohamed\Downloads\part4_sound_sample_2.wav'; %write your Own path
audiowrite(filename,normalized_y_t,Fs);
sound(normalized_y_t,Fs);                
clear normalized_y_t Fs
%%%%%%%%%%%%%%%%%%%%%%%%%%%'''PART5'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(t2*1e3,y_t);                        
title('Signal x(t) with WGN [y(t)]');
set(gca, 'Xlim', [0 ((100*R)+(20*n))]);
ylabel('Amplitude');
xlabel('Time (ms)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%'''PART6'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = 8000;  
nfft = length(y_t);
nfft2 = 2^nextpow2(nfft);
f = Fs*(1:(nfft2))/nfft2;
Y_f = abs(fft(y_t,nfft2));
figure;
plot(f,Y_f);                                 
set(gca, 'Xlim', [600 1700]);
ylabel('Amplitude');
xlabel('Frequency (Hz)');
title('Y(f)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%'''PART7'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spectrogram (y_t,16,8,2^14,Fs)
figure
spectrogram(y_t,rectwin(16),8,2^14,Fs,'yaxis');
title('Signal Y(t) with rectangular window size 16');
figure
spectrogram(y_t,rectwin(64),32,2^14,Fs,'yaxis');
title('Signal Y(t) with rectangular window size 64');
figure
spectrogram(y_t,rectwin(256),128,2^14,Fs,'yaxis');
title('Signal Y(t) with rectangular window size 256');
figure
spectrogram(x_t,rectwin(1024),512,2^14,Fs,'yaxis');
title('Signal Y(t) with rectangular window size 1024');
figure
spectrogram(x_t,rectwin(4096),2048,2^14,Fs,'yaxis');
title('Signal Y(t) with rectangular window size 4096');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
spectrogram(y_t,blackman(16),8,2^14,Fs,'yaxis');
title('Signal Y(t) with blackman window size 16');
figure
spectrogram(y_t,blackman(64),32,2^14,Fs,'yaxis');
title('Signal Y(t) with blackman window size 64');
figure
spectrogram(y_t,blackman(256),128,2^14,Fs,'yaxis');
title('Signal Y(t) with blackman window size 256');
figure
spectrogram(x_t,blackman(1024),512,2^14,Fs,'yaxis');
title('Signal Y(t) with blackman window size 1024');
figure
spectrogram(x_t,blackman(4096),2048,2^14,Fs,'yaxis');
title('Signal Y(t) with blackman window size 4096');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:R
A=1;
W = (R-i)*960;
C = [zeros(A,((i-1)*960)),ones(A,(800)),zeros(A,W)];
O = y_t.*C;
V=[];
buttonArray =['1' '2' '3' 
              '4' '5' '6'
              '7' '8' '9' 
              '*' '0' '#'];
original_f = [697 770 852 941 1209 1336 1477];
L = round(original_f/8000*205);  % Indices of the DFT
start=(((i-1)*960)+1);
tone = O(start:start+204);
TONE = abs(goertzel(tone,L+1)); % takes DFT of digit using goertzel algorithm
vert = TONE(1:4)>50;
horz = TONE(5:7)>20;
V = buttonArray(vert,horz);
phoneNum = cat(2,phoneNum,V);
end
 phoneNum 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''PART1'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = sym2TT(S)
Symbol = {'1','2','3','A','4','5','6','B','7','8','9','C','*','0','#','D'};
lfg = [697 770 852 941]; % Low frequency group
hfg = [1209 1336 1477 1633];  % High frequency group
f  = [];
for c=1:4,
    for r=1:4,
        f = [ f [lfg(c);hfg(r)] ];
    end
end
Fs  = 8000;       % Sampling frequency 8 kHz
N = 800;          % Tones of 100 ms
t1   = (0:N-1)/Fs; % 800 samples at Fs
syms t
pit = 2*pi*t; 
pit1 = 2*pi*t1;
IndexC = strfind(Symbol,S);
Index = find(not(cellfun('isempty',IndexC)));
x = sum(sin(f(:,Index)*pit));
x1 = sum(sin(f(:,Index)*pit1));

%%%% if we want to plot a certain symbol in time domain %%%%

%     plot(t1*1e3,x1);
%     title(['Symbol "', Symbol{Index},'": [',num2str(f(1,Index)),',',num2str(f(2,Index)),']';string(x)])
%     set(gca, 'Xlim', [0 25]);
%     ylabel('Amplitude');
%     xlabel('Time (ms)');
end 