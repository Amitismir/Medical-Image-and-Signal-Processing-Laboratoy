clc;
clear;
load("E:\6th Semester\MISP Lab\MyLab\LAB4\FiveClass_EEG.mat");
%% Part 3.1.1.
% Bandpass filter between 1 - 4
fs = 256;      % Sampling Frequency 
N      = 4;    % Order 
Fpass1 = 1;  % First Passband Frequency 
Fpass2 = 43;   % Second Passband Frequency 
Apass  = 1;    % Passband Ripple (dB) 
  
% Construct an FDESIGN object and call its CHEBY1 method. 
h  = fdesign.bandpass('N,Fp1,Fp2,Ap', N, Fpass1, Fpass2, Apass, fs); 
Hd = design(h, 'cheby1'); 
for c=1:30 
    Delta_X(:,c) = filter(Hd,X(:,c)); 
end 
% Bandpass filter between 4 - 8 
fs = 256;      % Sampling Frequency 
N      = 4;    % Order 
Fpass1 = 4;  % First Passband Frequency 
Fpass2 = 8;   % Second Passband Frequency 
Apass  = 1;    % Passband Ripple (dB) 
  
% Construct an FDESIGN object and call its CHEBY1 method. 
h  = fdesign.bandpass('N,Fp1,Fp2,Ap', N, Fpass1, Fpass2, Apass, fs); 
Hd = design(h, 'cheby1'); 
for c=1:30 
    Theta_X(:,c) = filter(Hd,X(:,c)); 
end 

% Bandpass filter between 8 - 13 
fs = 256;      % Sampling Frequency 
N      = 4;    % Order 
Fpass1 = 8;  % First Passband Frequency 
Fpass2 = 13;   % Second Passband Frequency 
Apass  = 1;    % Passband Ripple (dB) 
  
% Construct an FDESIGN object and call its CHEBY1 method. 
h  = fdesign.bandpass('N,Fp1,Fp2,Ap', N, Fpass1, Fpass2, Apass, fs); 
Hd = design(h, 'cheby1'); 
for c=1:30 
    Alpha_X(:,c) = filter(Hd,X(:,c)); 
end 

% Bandpass filter between 13 - 30 
fs = 256;      % Sampling Frequency 
N      = 4;    % Order 
Fpass1 = 13;  % First Passband Frequency 
Fpass2 = 30;   % Second Passband Frequency 
Apass  = 1;    % Passband Ripple (dB) 
  
% Construct an FDESIGN object and call its CHEBY1 method. 
h  = fdesign.bandpass('N,Fp1,Fp2,Ap', N, Fpass1, Fpass2, Apass, fs); 
Hd = design(h, 'cheby1'); 
for c=1:30 
    Beta_X(:,c) = filter(Hd,X(:,c)); 
end 
%% Part 3.1.2

t = 5;
n = fs * t;
dur = linspace(0 , t, n);

figure("Position",[0,0,1000,800]);
subplot(5,1,1);
plot(dur, X(1:n,1));
xlabel('Time (s)', 'Interpreter','latex');
ylabel('Amplitude (a.u.)', 'Interpreter','latex');
title('Original EEG Signal of Channel 1', 'Interpreter','latex');
grid on;
subplot(5,1,2);
plot(dur, Delta_X(1 : n, 1));
xlabel('Time (s)', 'Interpreter','latex');
ylabel('Amplitude (a.u.)', 'Interpreter','latex');
title('Delta-Filtered EEG Signal of Channel 1', 'Interpreter','latex');
grid on;
subplot(5,1,3);
plot(dur, Theta_X(1 : n, 1));
xlabel('Time (s)', 'Interpreter','latex');
ylabel('Amplitude (a.u.)', 'Interpreter','latex');
title('Theta-Filtered EEG Signal of Channel 1', 'Interpreter','latex');
grid on;
subplot(5,1,4);
plot(dur, Alpha_X(1:n, 1));
xlabel('Time (Samples)', 'Interpreter','latex');
ylabel('Amplitude (a.u.)', 'Interpreter','latex');
title('Alpha-Filtered EEG Signal of Channel 1', 'Interpreter','latex');
grid on;
subplot(5,1,5);
plot(dur,Beta_X(1:n, 1));
xlabel('Time (Samples)', 'Interpreter','latex');
ylabel('Amplitude (a.u.)', 'Interpreter','latex');
title('Beta-Filtered EEG Signal of Channel 1', 'Interpreter','latex');
grid on;
sgtitle('EEG Signal of First 5 Seconds From Channel 1 in Different Frequency Bands');

%% Part 3.2 
for i=1:200 
   Delta_Trials(:,:,i) = Delta_X(trial(i): trial(i)+256*10,:); 
   Theta_Trials(:,:,i) = Theta_X(trial(i): trial(i)+256*10,:); 
   Alpha_Trials(:,:,i) = Alpha_X(trial(i): trial(i)+256*10,:); 
   Beta_Trials(:,:,i) = Beta_X(trial(i): trial(i)+256*10,:); 
end 
%% Part 3.3
for i = 1:200
    for j = 1:30
        Power_Delta(:,j,i) = (Delta_Trials(:,j,i).^2);
        Power_Theta(:,j,i) = (Theta_Trials(:,j,i).^2);
        Power_Alpha(:,j,i) = (Alpha_Trials(:,j,i).^2);
        Power_Beta(:,j,i) = (Beta_Trials(:,j,i).^2);
    end
end
%% Part 3.4
clc;

for class = 1:5
    Delta_X_Avg(:,:,class) = mean(Power_Delta(:,:,(y==class)),3);
    Theta_X_Avg(:,:,class) = mean(Power_Theta(:,:,(y==class)),3);
    Alpha_X_Avg(:,:,class) = mean(Power_Alpha(:,:,(y==class)),3);
    Beta_X_Avg(:,:,class) = mean(Power_Beta(:,:,(y==class)),3);
end
%% Part 3.5
newWin = ones(1,200) / sqrt(200);
for class = 1:5
    for ch = 1:30
        Delta_Avg_Smooth(:,ch,class) = conv(Delta_X_Avg(:,ch, class), newWin, "same");
        Theta_Avg_Smooth(:,ch,class) = conv(Theta_X_Avg(:,ch, class), newWin, "same");
        Alpha_Avg_Smooth(:,ch,class) = conv(Alpha_X_Avg(:,ch, class), newWin, "same");
        Beta_Avg_Smooth(:,ch,class) = conv(Beta_X_Avg(:,ch, class), newWin, "same");
    end
end

%% Part 3.6
% CPz channel = 16
band_data = {
    Delta_Avg_Smooth, 'Delta';
    Theta_Avg_Smooth, 'Theta';
    Alpha_Avg_Smooth, 'Alpha';
    Beta_Avg_Smooth,  'Beta';
};
t = 10;
n_samples = size(Delta_Avg_Smooth, 1);   % = 2561
duration = linspace(0, t, n_samples);    % time vector in seconds
% Loop over each band
for b = 1:size(band_data,1)
    data = band_data{b,1};
    band_name = band_data{b,2};
    
    figure('Position',[0,0,1000,750]);
    for cls = 1:5
        subplot(5,1,cls);
        plot(duration,data(:,16,cls));  % 16 = channel index
        %xline(768, Color='r', LineWidth=1.5);
        xlabel('Time (s)', 'Interpreter','latex');
        ylabel('Amplitude (a.u.)', 'Interpreter','latex');
        title(['Class ' num2str(cls)], 'Interpreter','latex');
    end
    
    sgtitle([band_name ' Band - Channel 16'], 'Interpreter','latex');
end
