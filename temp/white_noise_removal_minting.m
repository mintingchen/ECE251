l = 10; % second
% how the white noise signal is generated?
% how to calculate the SNR?
% Perform 1-10 iteration and plot the output
% may need to record the clean sound again
file = 'input/clean.m4a';
samples = [l*48000*1.5,2.5*l*48000-1];
[x_data,fs] = audioread(file,samples);
t = 0:1/fs:l-1/fs;
x_data= x_data(:,1);
plot(t,x_data)
sound(x_data,fs)
%%
lev = 10; 
snr_list = zeros(10,1);

for i = 1:lev
    y = awgn(x_data,5,'measured'); % Add white Gaussian noise. figure(1)
    %plot(t,y,t,x_data)
    wname = 'coif5'; 
    tree = wpdec(y,lev,wname); 
    det1 = wpcoef(tree,2014);
    sigma = median(abs(det1))/0.6745; 
    alpha = 1.8;
    thr = wpbmpen(tree,sigma,alpha);
    keepapp = 1;
    xd = wpdencmp(tree,'s','nobest',thr,keepapp);
    noise = xd-x_data;
    snr_list(i)= 10*log10(sum(x_data.^2)/sum(noise.^2));
end
%%
noise_xy = y-x_data;
ori_snr = 10*log10(sum(x_data.^2)/sum(noise_xy.^2))
%%
index = linspace(1,10,10);
plot(index, snr_list)

% figure
% plot(t,x_data,t,xd,'g') 
% xlabel('time') 
% ylabel('Amplitude') 
% title('signal denoise');
%% plot bar diagram 
[breath_clean,fs] = audioread(file,samples);
[breath_movie,fs] = audioread('movie.m4a');
[breath_music,fs] = audioread('music.m4a');
[breath_white_low1,fs] = audioread('white.m4a');
%[breath_white_low2,fs] = audioread('breath_white_lowersnr2.m4a');
%[breath_white_low3,fs] = audioread('breath_white_lowersnr3.m4a');
standard_len = length(breath_clean);
audio = {x_data,breath_movie,breath_music,breath_white_low1,2/3*breath_white_low1,1/2*breath_white_low1,y};
for i = 2:length(audio)-1
    audio{i} = audio{i}(:,1);
    audio{i} = audio{i}(1:standard_len);
    audio{i} = audio{i}*3 + x_data;
end
audio_snr = [];
for i = 2:length(audio)
    noise = abs(audio{i} - x_data);
    snr = 10*log10(sum(x_data.^2)/sum(noise.^2));
    audio_snr = [audio_snr,snr];
end
%%
audio_snr_denoise = [];
audio_denoise = {};
for i = 2:length(audio)
    wname = 'coif5'; 
    lev = 10
    tree = wpdec(audio{i},lev,wname); 
    det1 = wpcoef(tree,2014);
    sigma = median(abs(det1))/0.6745; 
    alpha = 1.8;
    thr = wpbmpen(tree,sigma,alpha);
    keepapp = 1;
    audio_denoise{i} = wpdencmp(tree,'s','nobest',thr,keepapp);
    noise = abs(audio_denoise{i}-x_data);
    audio_snr_denoise(i)= 10*log10(sum(x_data.^2)/sum(noise.^2));
end
%%
X = categorical({'movie','music','low white noise','medium white noise','high white noise','synthesized white noise'});
X = reordercats(X,{'movie','music','low white noise','medium white noise','high white noise','synthesized white noise'});
snr_bar = [audio_snr(1) audio_snr_denoise(2); audio_snr(2) audio_snr_denoise(3);audio_snr(3) audio_snr_denoise(4);...
    audio_snr(4) audio_snr_denoise(5);audio_snr(5) audio_snr_denoise(6);audio_snr(6) audio_snr_denoise(7)];
bar(X,snr_bar)

