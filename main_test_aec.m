clear;
clc; close all;
frameSize = 2048;
%% Load near and far speech signal
% Audio 1
% fs = 8000;
% load nearspeech;
% nearspeech = v;
% load farspeech;
% farspeech = x;

% Audio 2
fs = 16000;
v1 = audioread('timitcorpus/timit/timit/dr8-mbcg0/sa1.wav');
x1 = audioread('timitcorpus/timit/timit/dr8-fbcg1/sa2.wav');
v2 = audioread('timitcorpus/timit/timit/dr8-mbcg0/si486.wav');
x2 = audioread('timitcorpus/timit/timit/dr8-fbcg1/si982.wav');
v3 = audioread('timitcorpus/timit/timit/dr8-mbcg0/si957.wav');
x3 = audioread('timitcorpus/timit/timit/dr8-fbcg1/si1612.wav');
v4 = audioread('timitcorpus/timit/timit/dr8-mbcg0/si2217.wav');
x4 = audioread('timitcorpus/timit/timit/dr8-fbcg1/si2242.wav');
v5 = audioread('timitcorpus/timit/timit/dr8-mbcg0/sx57.wav');
x5 = audioread('timitcorpus/timit/timit/dr8-fbcg1/sx82.wav');
v6 = audioread('timitcorpus/timit/timit/dr8-mbcg0/sx147.wav');
x6 = audioread('timitcorpus/timit/timit/dr8-fbcg1/sx172.wav');
farspeech = v1;
farspeech(end+1:end+length(v2)) = v2;
farspeech(end+1:end+length(x1)) = 0;
farspeech(end+1:end+length(x2)) = 0;
farspeech(end+1:end+length(v3)) = v3;
farspeech(end+1:end+length(v4)) = v4;
farspeech(end+1:end+length(x3)) = 0;
farspeech(end+1:end+length(x4)) = 0;
farspeech(end+1:end+length(v5)) = v5;
farspeech(end+1:end+length(v6)) = v6;
farspeech(end+1:end+length(x5)) = 0;
farspeech(end+1:end+length(x6)) = 0;
farspeech = farspeech*2;
x = farspeech;

nearspeech(1:length(v1)) = 0;
nearspeech(end+1:end+length(v2)) = 0;
nearspeech(end+1:end+length(x1)) = x1;
nearspeech(end+1:end+length(x2)) = x2;
nearspeech(end+1:end+length(v3)) = 0;
nearspeech(end+1:end+length(v4)) = 0;
nearspeech(end+1:end+length(x3)) = x3;
nearspeech(end+1:end+length(x4)) = x4;
nearspeech(end+1:end+length(v5)) = 0;
nearspeech(end+1:end+length(v6)) = 0;
nearspeech(end+1:end+length(x5)) = x5;
nearspeech(end+1:end+length(x6)) = x6;
nearspeech = nearspeech*4;
nearspeech = nearspeech';
v = nearspeech;

% sound(farspeech, fs);
% pause(length(farspeech)/fs + 1);
% sound(nearspeech, fs);
% pause(length(nearspeech)/fs + 1);

%% Generate echo of the far speech signal
% delay = 0.05; % seconds
% alpha = 0.1; % echo strength
% farSpeechEcho = echo_gen(farspeech, alpha, delay, fs);
% [b,a] = fir1(20, 0.5);
% farSpeechEcho = filter(b, a, farspeech);

% Impulse response = RIR
M = fs/2 + 1;
[B, A] = cheby2(4,20,[0.1 0.7]);
impulseResponseGenerator = dsp.IIRFilter(...
    'Numerator', [zeros(1,6) B], ...
    'Denominator', A);
roomImpulseResponse = impulseResponseGenerator( ...
(log(0.99*rand(1,M)+0.01).*sign(randn(1,M)).*exp(-0.05*(1:M)))');
roomImpulseResponse = roomImpulseResponse/norm(roomImpulseResponse)*4;
figure;plot(roomImpulseResponse); title('Room Impulse Response'); xlabel('samples'); ylabel('Amplitude');grid on;
room = dsp.FIRFilter('Numerator', roomImpulseResponse'); 
farSpeechEcho = room(farspeech);
%% Generate mic signal
micSignal = farSpeechEcho + nearspeech + 0.001*randn(length(nearspeech), 1);
% sound(micSignal, fs);
% pause(length(micSignal)/fs +1);

%% Echo Return Loss Enhancement (ERLE) setup
% Measure of how much the echo has been attenuated
diffAverager = dsp.FIRFilter('Numerator', ones(1,1024));
farEchoAverager = clone(diffAverager);
diffAverager_total = clone(diffAverager);
farEchoAverager_total = clone(diffAverager);

player = audioDeviceWriter('SupportVariableSizeInput', true, ...
                                    'BufferSize', 512, 'SampleRate', fs);

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
AECScope1.Title         = 'Clean signal from AEC';

AECScope1.ActiveDisplay = 4;
AECScope1.ShowGrid      = true;
AECScope1.YLimits       = [0 50];
AECScope1.YLabel        = 'ERLE (dB)';
AECScope1.Title         = 'Echo Return Loss Enhancement';

%% Adaptive filter
num_frames = floor(length(nearspeech)/frameSize);
% orders = [2, 4, 8, 16, 32, 64, 128, 256, 512];
% step_sizes = 0.05:0.05:1;
orders = [128];
step_sizes = [0.05];
erle_mat = zeros(length(orders), length(step_sizes));
for i = 1:length(orders)
    for j = 1:length(step_sizes)
        filter_order = orders(i);
        mu = step_sizes(j);
        echoCanceller = dsp.LMSFilter('Method','Normalized LMS' , ...
            'Length',filter_order, ...
            'StepSize',mu);
        % echoCanceller    = dsp.FrequencyDomainAdaptiveFilter('Length', 64, ...
        %                     'StepSize', 0.025, ...
        %                     'InitialPower', 0.01, ...
        %                     'AveragingFactor', 0.98, ...
        %                     'Method', 'Unconstrained FDAF');
        clean_sig = zeros(length(nearspeech), 1);
        for k = 1:num_frames
            start_idx = (k-1)*frameSize+1;
            end_idx = k*frameSize;
            near = nearspeech(start_idx:end_idx);
            far = farspeech(start_idx:end_idx);
            mic = micSignal(start_idx:end_idx);
            far_echo = farSpeechEcho(start_idx:end_idx);
            % Apply NLMS filter
            [y, e] = echoCanceller(far, mic);
            player(e);
            clean_sig(start_idx:end_idx) = e;
            % Compute ERLE
            erle = diffAverager((e - near).^2)./ farEchoAverager(far_echo.^2);
            erledB = -10*log10(erle);
            AECScope1(near, mic, e, erledB);
        end
        release(AECScope1);
        erle = diffAverager_total((clean_sig - nearspeech).^2)./ farEchoAverager_total(farSpeechEcho.^2);
        erledB = -10*log10(erle);
        erle_mat(i,j) = mean(erledB(find(erledB>0)));
    end
end
% T = array2table(erle_mat, 'VariableNames', string(step_sizes), 'RowNames', string(orders));
% writetable(T, 'nlms_params_audio2.csv', 'WriteRowNames',true);
audiowrite('clean_signal_audio2.wav', clean_sig, fs);
audiowrite('micSignal_audio2.wav', micSignal, fs);

max_val = max(max(erle_mat));
[iopt, jopt] = find(erle_mat == max_val);
disp("Optimum filter order");disp(orders(iopt));
disp("Optimum step size");disp(step_sizes(jopt));

% e = e';
% erle = diffAverager((e - nearspeech).^2)./ farEchoAverager(farSpeechEcho.^2);
% erledB = -10*log10(erle);
% AECScope1(nearspeech, micSignal, e, erledB);
% release(AECScope1);

%% Functions
% function y = echo_gen(x, alpha, delay, fs)
%     D = delay * fs;  
%     y = zeros(size(x));  
%     y(1:D) = x(1:D);     
%      for i=D+1:length(x)  
%        y(i) = x(i) + alpha*x(i-D);  
%      end  
% end