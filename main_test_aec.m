clear;
clc; close all;
frameSize = 2048;
%% Load near speech signal
% [v1, fs] = audioread('timitcorpus/timit/timit/dr8-mbcg0/sa1.wav');
% [x1, fs] = audioread('timitcorpus/timit/timit/dr8-fbcg1/sa2.wav');
% [v2, fs] = audioread('timitcorpus/timit/timit/dr8-mbcg0/si486.wav');
% [x2, fs] = audioread('timitcorpus/timit/timit/dr8-fbcg1/si982.wav');
% [v3, fs] = audioread('timitcorpus/timit/timit/dr8-mbcg0/si957.wav');
% [x3, fs] = audioread('timitcorpus/timit/timit/dr8-fbcg1/si1612.wav');
% [v4, fs] = audioread('timitcorpus/timit/timit/dr8-mbcg0/si2217.wav');
% [x4, fs] = audioread('timitcorpus/timit/timit/dr8-fbcg1/si2242.wav');
% [v5, fs] = audioread('timitcorpus/timit/timit/dr8-mbcg0/sx57.wav');
% [x5, fs] = audioread('timitcorpus/timit/timit/dr8-fbcg1/sx82.wav');
% [v6, fs] = audioread('timitcorpus/timit/timit/dr8-mbcg0/sx147.wav');
% [x6, fs] = audioread('timitcorpus/timit/timit/dr8-fbcg1/sx172.wav');
% farspeech = v1;
% farspeech(end+1:end+length(v2)) = v2;
% farspeech(end+1:end+length(x1)) = 0;
% farspeech(end+1:end+length(x2)) = 0;
% farspeech(end+1:end+length(v3)) = v3;
% farspeech(end+1:end+length(v4)) = v4;
% farspeech(end+1:end+length(x3)) = 0;
% farspeech(end+1:end+length(x4)) = 0;
% farspeech(end+1:end+length(v5)) = v5;
% farspeech(end+1:end+length(v6)) = v6;
% farspeech(end+1:end+length(x5)) = 0;
% farspeech(end+1:end+length(x6)) = 0;
% farspeech = farspeech';
% sound(nearspeech, fs);
% pause(10);

% nearspeech(1:length(v1)) = 0;
% nearspeech(end+1:end+length(v2)) = 0;
% nearspeech(end+1:end+length(x1)) = x1;
% nearspeech(end+1:end+length(x2)) = x2;
% nearspeech(end+1:end+length(v3)) = 0;
% nearspeech(end+1:end+length(v4)) = 0;
% nearspeech(end+1:end+length(x3)) = x3;
% nearspeech(end+1:end+length(x4)) = x4;
% nearspeech(end+1:end+length(v5)) = 0;
% nearspeech(end+1:end+length(v6)) = 0;
% nearspeech(end+1:end+length(x5)) = x5;
% nearspeech(end+1:end+length(x6)) = x6;
% sound(farspeech, fs);
% pause(10);
load nearspeech;
nearspeech = v;
player = audioDeviceWriter('SupportVariableSizeInput', true, ...
                                    'BufferSize', 512, 'SampleRate', fs);
%% Load far speech signal
load farspeech;
farspeech = x;
%% Generate echo of the far speech signal
delay = 0.2; % seconds
alpha = 0.1; % echo strength
farSpeechEcho = echo_gen(farspeech, alpha, delay, fs);

% Impulse response = delta function
% roomImpulseResponse = zeros(1, M); roomImpulseResponse(1) = 1;
% room = dsp.FIRFilter('Numerator', roomImpulseResponse);
% farSpeechEcho = room(farspeech);

% Impulse response = RIR
% [B, A] = cheby2(4,20,[0.1 0.7]);
% impulseResponseGenerator = dsp.IIRFilter(...
%     'Numerator', [zeros(1,6) B], ...
%     'Denominator', A);
% roomImpulseResponse = impulseResponseGenerator( ...
% (log(0.99*rand(1,M)+0.01).*sign(randn(1,M)).*exp(-0.002*(1:M)))');
% roomImpulseResponse = roomImpulseResponse/norm(roomImpulseResponse)*4;
% room = dsp.FIRFilter('Numerator', roomImpulseResponse'); 
% farSpeechEcho = room(farspeech);
%% Generate mic signal
micSignal = farSpeechEcho + nearspeech + 0.001*randn(length(nearspeech), 1);
sound(micSignal, fs);
pause(10);
%% Adaptive filter
alpha = 0.002;
c = 0.001;
w = zeros(1, length(farspeech));
for i = 1:length(farspeech)
    tic
    mu = alpha/(c+(x(i)'*x(i)));
    e(i) = micSignal(i) - w(i)'* x(i);
    w(i+1) = w(i) + mu * e(i) * x(i);
    toc
    timeForEachIteration(i) = toc;
end
sound(e);
pause(10);

% filter_order = 32;
% mu = 0.005;
% echoCanceller = dsp.LMSFilter('Method','Normalized LMS' , ...
%     'Length',filter_order, ...
%     'StepSize',mu);
% echoCanceller    = dsp.FrequencyDomainAdaptiveFilter('Length', 2048, ...
%                     'StepSize', 0.025, ...
%                     'InitialPower', 0.01, ...
%                     'AveragingFactor', 0.98, ...
%                     'Method', 'Unconstrained FDAF');

%% Echo Return Loss Enhancement (ERLE) 
% Measure of how much the echo has been attenuated
diffAverager = dsp.FIRFilter('Numerator', ones(1,1024));
farEchoAverager = clone(diffAverager);
AECScope1   = dsp.TimeScope(4, fs, ...
                'LayoutDimensions', [4,1], ...
                'TimeSpan', length(nearspeech)/fs, 'TimeSpanOverrunAction', 'Scroll', ...
                'BufferLength', length(nearspeech));      
AECScope1.ActiveDisplay = 1;
AECScope1.ShowGrid      = true;
AECScope1.YLimits       = [-1 1];
AECScope1.Title         = 'Near-End Speech Signal';

AECScope1.ActiveDisplay = 2;
AECScope1.ShowGrid      = true;
AECScope1.YLimits       = [-1 1];
AECScope1.Title         = 'Microphone Signal';

AECScope1.ActiveDisplay = 3;
AECScope1.ShowGrid      = true;
AECScope1.YLimits       = [-1 1];
AECScope1.Title         = 'Output of Acoustic Echo Canceller';

AECScope1.ActiveDisplay = 4;
AECScope1.ShowGrid      = true;
AECScope1.YLimits       = [0 50];
AECScope1.YLabel        = 'ERLE (dB)';
AECScope1.Title         = 'Echo Return Loss Enhancement';

num_frames = floor(length(nearspeech)/frameSize);

% for i = 1:num_frames
%     start_idx = (i-1)*frameSize+1;
%     end_idx = i*frameSize;
%     near = nearspeech(start_idx:end_idx)';
%     far = farspeech(start_idx:end_idx)';
%     mic = micSignal(start_idx:end_idx)';
%     far_echo = farSpeechEcho(start_idx:end_idx)';
%     % Apply NLMS filter
%     [y, e] = echoCanceller(far, mic);
%     player(e);
%     % Compute ERLE
%     erle = diffAverager((e - near).^2)./ farEchoAverager(far_echo.^2);
%     erledB = -10*log10(erle);
%     AECScope1(near, mic, e, erledB);
% end
% release(AECScope1);
e = e';
erle = diffAverager((e - nearspeech).^2)./ farEchoAverager(farSpeechEcho.^2);
erledB = -10*log10(erle);
AECScope1(nearspeech, micSignal, e, erledB);
release(AECScope1);

%% Functions
function y = echo_gen(x, alpha, delay, fs)
    D = delay * fs;  
    y = zeros(size(x));  
    y(1:D) = x(1:D);     
     for i=D+1:length(x)  
       y(i) = x(i) + alpha*x(i-D);  
     end  
end