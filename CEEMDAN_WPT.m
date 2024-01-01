clc;
clear;
close all;
%%
Total = readtable('……\148047_HBA_Probe1_Deoxy.csv');
ch1=Total{1:1250,7};
Fs=10;     % 采集频率
T=1/Fs;        % 采集时间间隔

%信号长度最好为偶数
N=length(ch1); % 采集信号的长度

y = ch1; %将数据放入y向量

t=(0:1:N-1)*T;   % 定义整个采集时间点
t=t';            % 转置成列向量

%% 绘制时域信号
figure(1)
plot(t,y,'Color','r');
xlabel('Time(s)')
ylabel('Amplitude')
title('x(t)')

%% 绘制频域信号


delta_f = 1*Fs/N;
 x=ch1;
X = fftshift(abs(fft(x)))/N;
X_angle = fftshift(angle(fft(x)));
f = (-N/2:N/2-1)*delta_f;
 
figure(2);
subplot(3,1,1);
plot(t,x);
title('原信号');
subplot(3,1,2);
plot(f,X);
grid on;
title('原信号频谱幅度特性');
subplot(3,1,3);
plot(f,X_angle);
title('原信号频谱相位特性');

%% CEEMDAN-WPT
% CEEMDAN分解
Nstd = 0.2;                % 高斯白噪声标准差
NR = 100;                  % 加入噪声的次数
MaxIter = 200;            % 最大迭代次数
[imf, its]=ceemdan(y,Nstd,NR,MaxIter);
[m, n]=size(imf);  % 分解出来的imf数量不稳定，在这里会出现14或者15个
CC=zeros(1,m);  % 相关系数


%画高频分量
figure;title('high-frequency IMF component');
for i=1:5
    subplot(5,1,i);plot(t,imf(i,:));ylabel(['IMF',num2str(i)]);xlabel('Time(s)'); %title('high-frequency IMF component')
    CC(i)=corr(imf(i,:)',y,'type','Pearson');   % 相关系数
end

%画低频分量
figure;title('low-frequency IMF component');
for i=6:10
    subplot(5,1,i-5);plot(t,imf(i,:),'Color','#D95319');ylabel(['IMF',num2str(i)]);xlabel('Time(s)');%title('low-frequency IMF component')
    CC(i)=corr(imf(i,:)',y,'type','Pearson');   % 相关系数
end

figure;
plot(CC,'k-o');xlabel('分解阶次');ylabel('相关系数');

%对前3个IMF分量进行小波阈值去噪
imf=imf';
idx=5;
Y_imf=zeros(N,idx);
for i=1:idx
    y_imf=imf(:,i);
    % 用db4小波对含噪信号进行3层分解并提取系数
    [c,l]=wavedec(y_imf,3,'db4'); 
    
    %取第3层低频近似系数
    ca3=appcoef(c,l,'db4',3);
    
    %取各层高频细节系数
    cd3=detcoef(c,l,3);
    cd2=detcoef(c,l,2);
    cd1=detcoef(c,l,1);
    
    % 阈值获取
    thr=thselect(y,'rigrsure');
    
    % 进行软阈值处理
    ysoft3=wthresh(cd3,'s',thr);
    ysoft2=wthresh(cd2,'s',thr);
    ysoft1=wthresh(cd1,'s',thr);
    c1=[ca3;ysoft3;ysoft2;ysoft1];
    
    % 小波重构
    Y_imf(:,i)=waverec(c1,l,'db4');
end

%小波重构后的信号
figure;title('high-frequency IMF component after denoising)');
for i=1:5
    subplot(5,1,i);plot(t,Y_imf(:,i),'Color','b');ylabel(['IMF',num2str(i)]);xlabel('Time(s)'); %title(['high-frequency IMF component after denoising)',1])
%     CC(i)=corr(imf(i,:)',y,'type','Pearson');   % 相关系数
end

% 舍弃imf9及其后面的分量，将去噪完成的模态分量(imf1—imf4)与保留的未经处理的模态分量(imf5—imf8)重构得到纯净信号Y(t)
ceemdY=sum(Y_imf,2)+sum(imf(:,6:m),2);
% Y=sum(imf(:,6:8),2);

%%  画对比图
figure;
plot(t,ceemdY,'Color','m')
xlabel('Time(s)')
ylabel('Amplitude')
title('x’(t)')

figure;
subplot(2,1,1);plot(t,y);xlabel('Time(s)');ylabel('Amplitude');title('x(t)');axis([0 130 min(y)-0.01 max(y)+0.01]);
subplot(2,1,2);plot(t,ceemdY);xlabel('Time(s)');ylabel('Amplitude');title('x’(t)');axis([0 130 min(ceemdY)-0.01 max(ceemdY)+0.01]);

%% 对比
figure;
plot(t,y,'LineWidth',1.5);xlabel('Time(s)');ylabel('Amplitude');
hold on
plot(t,ceemdY,'red','LineWidth',1.5);
legend('x(t)','x’(t)')
%% SNR
sigPower = sum(abs(y).^2) / length(y);             % 求出信号功率
noisePower = sum(abs(ceemdY - y).^2) / length(ceemdY - y);   % 求出噪声功率
SNR= 10*log10(sigPower/noisePower) 

%% RMSE
rmse =sqrt(mean(abs(ceemdY - y).^2))