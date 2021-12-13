clear all
l = 10; % second
clean_file = 'input/clean.m4a';
noise_file = 'input/music.m4a';
samples = [l*48000*1.5,2.5*l*48000-1];
[clean,fs] = audioread(clean_file,samples);
clean= (clean(:,1)-0.00065);
[noise,fs] = audioread(noise_file,samples);
noise = noise(:,1)*0.05;

input = clean + noise;
snr_before = 10*log10(sum(clean.^2)/sum(noise.^2))
% sound(input,fs)

t = 0:1/fs:l-1/fs;
iter = 10;
snr_iter = zeros(iter,1);

choose_channel = [1,2;1,2;2,3;2,3;3,5;3,5;3,5;3,5;3,5;3,5];
for i = 1:10
    mra_data = MRA(input, i, l,fs);
    xd = zeros(480000,1);
    for j = 1:2
        xd = xd + mra_data(choose_channel(i,j),:)';
    end
    %sound(xd, fs)
    % enlarge the breathing sessions
    %xd(abs(xd)>0.00045) = xd(abs(xd)>0.00045)*2;
    xd = xd*3;
    if i==10
        sound(xd,fs)
        figure
        subplot(2,1,1)
        plot(t, input,t,clean);
        legend('input','clean')
        xlabel('time') 
        ylabel('Magnitude') 
        ylim([-5*10^-3 10*10^-3])
        subplot(2,1,2)
        plot(t, input,t,xd);
        legend('input','output')
        xlabel('time') 
        ylabel('Magnitude') 
        ylim([-5*10^-3 10*10^-3])
    end
    noise_after = abs(xd - clean);
    snr_after = 10*log10(sum(clean.^2)/sum(noise_after.^2))
    snr_iter(i) = snr_after;
end
snr_iter = [snr_before; snr_iter];

figure
plot(0:iter, snr_iter);
tex = text(0:iter,snr_iter,num2str(snr_iter));
xlabel('iter') 
ylabel('SNR') 
title('SNR vs iterations')

% figure
% plot(t, input, t, xd);
% legend('input','output')
% xlabel('time') 
% ylabel('Magnitude') 

%sound(xd,fs)


function mra = MRA(data, iter, l, fs)
% Separating Signal Components in Time
t = 0:1/fs:l-1/fs;
% Multiresolution analysis
mra = modwtmra(modwt(data,'coif5',iter),'coif5');
helperMRAPlot(data,mra,t,'wavelet','Wavelet MRA',[2 3 4 9])
end
function xd = wavelet_denoising(input,lev)
wname = 'coif5'; 
tree = wpdec(input,lev,wname); 
det1 = wpcoef(tree,2014);
sigma = median(abs(det1))/0.6745; 
alpha = 1.8;
thr = wpbmpen(tree,sigma,alpha);
keepapp = 1;
xd = wpdencmp(tree,'s','nobest',thr,keepapp);
end


