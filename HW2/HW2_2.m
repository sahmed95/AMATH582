% Part 2
%----------------------------------------------------------------------

% Piano

tr_piano=16; % record time in seconds
y= audioread('music1.wav'); 
Fs=length(y)/tr_piano; % sampling rate
figure(1)
subplot(2,1,1)
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano) signal'); 



y = y'/2; 
n = length(y);  % length of signal
L = tr_piano;    % final time 
t = linspace(0,L, n); % discretized time vector
k = (2*pi/L)*[0:n/2-1 -n/2:-1];  %rescale 
ks = fftshift(k); %shift

yt = fft(y);  %Fourier transform

subplot(2,1,2)
plot(ks, abs(fftshift(yt))/max(abs(yt)));
xlabel('Frequency'), ylabel('FFT(y)')
title('FFT of signal')

%Gabor filtering 

a = 500; % width parameter
dt = 0.25;
ygt_spec = [];
tslide = 0:dt:t(end-1);
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); % Gabor
    yg=g.*y; ygt=fft(yg);
    ygt_spec=[ygt_spec; abs(fftshift(ygt))];
    
end

figure(2)
pcolor(tslide,ks/(2*pi),(ygt_spec).'), shading interp
xlabel('Time(sec)'); ylabel('Frequency (Hz)')
title('Spectogram (Piano), a =500, dt = 0.1')
set(gca, 'Fontsize', [15]),
ylim([0, 950])
colormap(hot)

%-----------------------------------------------------------------------

% Recorder

tr_rec=14; % record time in seconds
y=audioread('music2.wav'); 
Fs=length(y)/tr_rec;
figure(3)
subplot(2,1,1)
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder) signal');
y = y'/2; 
n = length(y);  % length of signal
L = tr_rec;    % final time 
t = linspace(0,L, n); % discretized time vector
k = (2*pi/L)*[0:n/2-1 -n/2:-1];  %rescale 
ks = fftshift(k); %shift

yt = fft(y);  %Fourier transform

subplot(2,1,2)
plot(ks, abs(fftshift(yt))/max(abs(yt)));
xlabel('Frequency'), ylabel('FFT(y)')
title('FFT of signal')

%Gabor filtering 

a = 500; % width parameter
dt = 0.1;
ygt_spec = [];
tslide = 0:dt:t(end-1);
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); % Gabor
    yg=g.*y; ygt=fft(yg);
    ygt_spec=[ygt_spec; abs(fftshift(ygt))];
    
end

figure(4)
pcolor(tslide,ks/(2*pi),(ygt_spec).'), shading interp
xlabel('Time(sec)'); ylabel('Frequency (Hz)')
title('Spectogram (Recorder), a =500, dt = 0.1')
set(gca, 'Fontsize', [15]),
ylim([0, 2000])
colormap(hot)

