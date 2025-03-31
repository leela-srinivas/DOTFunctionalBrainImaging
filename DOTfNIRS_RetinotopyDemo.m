%% ESE 488 DOT/fNIRS Retinotopy Data set - Intro

close all; clear
load('NeuroDOT_Data_Sample_CCW1.mat'); % data, info, flags

%% plot grid
figure;
hold on;
title("3D Source-Detector Layout Map")
plot3(info.optodes.dpos3(:,1),info.optodes.dpos3(:,2),info.optodes.dpos3(:,3),'bx', "MarkerSize", 3)
plot3(info.optodes.spos3(:,1),info.optodes.spos3(:,2),info.optodes.spos3(:,3),'ro', "MarkerSize", 3)
text(info.optodes.spos3(1,1),info.optodes.spos3(1,2),info.optodes.spos3(1,3),'Src 1 (lower left)')
legend(["Detector Positions", "Source Positions"])
axis image, rotate3d 

%% show data traces, light falloff
figure('Position',[100 100 750 400])
subplot(1,2,1)
semilogy(data(info.pairs.r3d < 30,:)')
title("Data Traces, source-detector distance < 30 mm ")
ylabel('Intensity, \muW'), xlabel('t, samples')

subplot(1,2,2)
semilogy(info.pairs.r3d,mean(data,2),'.')
title("Light Falloff")
xlabel("source-detector distance, mm")
ylabel("Mean Signal Intensity")

%% log-ratio differential data
y = -log(bsxfun(@times,data,1./mean(data,2)));

maxDistance = 40; % this plot only
figure('Position',[150 150 750 400])
subplot(1,2,1), plot(y(info.pairs.r3d<maxDistance,:)') 
title(['Differential Log-Ratio Signals, r_{sd} < ' num2str(maxDistance) ' mm'])
xlabel("Time, samples"), ylabel('y = - log(I/\mu_I)')

pause

% Same signals displayed as an image 
subplot(1,2,2), imagesc(y(info.pairs.r3d<maxDistance,:)), caxis([-.2 .2]), colorbar
title(['Differential Log-Ratio Signals, r_{sd} < ' num2str(maxDistance) ' mm'])
xlabel("Time, samples"), ylabel('Measurements')

%% Spectrum - all measurements
fs = info.system.framerate;
f = [0:fs/size(y,2):fs/2]; % frequency vector 
Y = fft(y,[],2); % fft along 2nd dimension of y, all measurements 

figure('Position',[75 75 500 500])
subplot(2,2,1), plot(f,abs(Y(:,1:length(f))))
title('Spectrum - all measurements')
xlabel('f (Hz)'), ylabel('|Y_{sd}(f)|')
pause

% Average measurement spectrum
subplot(2,2,2), plot(f,mean(abs(Y(:,1:length(f))),1))
xlabel('f (Hz)'), ylabel('|Y_{sd}(f)|'), title('Average magnitude spectrum')
pause
set(gca,'XScale','log')
pause

% just WL 2, NN 2
subplot(2,2,3), plot(f,mean(abs(Y(info.pairs.NN==2 & info.pairs.WL == 2,1:length(f))),1))
xlabel('f (Hz)'), ylabel('|Y_{sd}(f)|'), title('Average spectrum: 2nd distance meas, 850 nm only')
set(gca,'XScale','log')
pause

% image plot, these meas
subplot(2,2,4), 
imagesc(y(info.pairs.NN==2 & info.pairs.WL == 2,:)), caxis([-.2 .2])
title(['Differential Log-Ratio Signals, 2nd distance meas, 850 nm only'])
xlabel("Time, samples"), ylabel('Measurements')
hold on
%show stimulus repetition start times 
for i = 1:length(info.paradigm.Pulse_2), xline(info.paradigm.synchpts(info.paradigm.Pulse_2(i))), end





