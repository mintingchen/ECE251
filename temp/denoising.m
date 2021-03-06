length = 10; % second
iter = 10; % number of channals

file = "input/breath_in_movie.m4a";
samples = [1,length*48000];
[x_data,fs] = audioread(file,samples);
t = 0:1/fs:length-1/fs;
x_data= x_data(:,1);
%sound(x_data,fs)

% Multiresolution analysis
mra_data = MRA(x_data, iter, length, fs); 
y_data = sum(mra_data(1:3,:));
figure
plot(t, x_data, t, y_data);
legend('input','output')
xlabel('time') 
ylabel('Magnitude') 

sound(y_data*3,fs)

function mra = MRA(data, iter, length, fs)
% Separating Signal Components in Time
t = 0:1/fs:length-1/fs;
% Multiresolution analysis
mra = modwtmra(modwt(data,iter),1);
helperMRAPlot(data,mra,t,'wavelet','Wavelet MRA',[2 3 4 9])
end

% Empirical Wavelet Transform
% t = 0:1/fs:length-1/fs;
% [mra_ewt,cfs,wfb,info] = ewt(x_data,'MaxNumPeaks',iter);
% helperMRAPlot(x_data,mra_ewt,t,'EWT','Empirical Wavelet Transform',[1 2 3 5])
% y_data = sum(mra_ewt(:,8),2);
% figure
% plot(t, x_data, t, y_data);
% legend('input','output')
% xlabel('time') 
% ylabel('Magnitude') 
