% For the LSE, the MSE performance versus the iteration for 
% B = 1,2,3 and infinity case
% Code is written by Jiang Zhu, Qi Zhang and Xiangming Meng, 2019.08.05
% This code is built based on the VALSE (written by Mihai-Alin Badiu)
clear;
close all;
clc;  
rng(1)
%% Generate a noisy superposition of K complex sinusoids with angular frequencies in [-pi,pi)
N    = 100;         % size of the full data
M    = 100;          % number of measurements with indices randomly chosen from 0,...,N-1
K    = 3;           % number of sinusoidal components
d    = 2*pi/N;      % minimum separation of the angular frequencies
SNR  = 0;          % signal to noise ratio (dB)
B_all = [1;2;3;inf];
% B_all = 1;

% Algorithm parameters
tau = [];
yy_min = [];
alpha = [];
Iter_max = 1000;
T_VALSE = Iter_max;
MC = 2;
MSE_Iter = nan(Iter_max,length(B_all),MC);
MSE_VALSE = nan(Iter_max,length(B_all),MC);
% Generate random measurement indices
tmp  = randperm(N-1)';
Mcal = [sort(tmp(1:M-1),'ascend'); N]-1;    % indices of the measurements
% Generate K frequencies with minimum separation d
% method_EP = 'uniform_EP' or 'diagonal_EP', which uses the scalar EP or diagonal EP;
% For low SNR scenario with small number of measurments, 'scalar EP ' is
% more prefereable
method_EP = 'diagonal_EP';
for mc = 1:MC
    waitbar(mc/MC)
    omega = zeros(K,1);
    omega(1) = pi * (2*rand - 1);
    for k = 2:K
        th = pi * (2*rand - 1);
        while min(abs((wrapToPi(th-omega(1:k-1))))) < d
            th = pi * (2*rand - 1);
        end
        omega(k) = th;
    end
%     omega = [-1;2];
    A   = exp(1i*(0:1:Mcal(end)).'*omega.');         % matrix with columns a(omega(1)),..., a(omega(K))
    r   = 1 + .2.*randn(K,1);                        % magnitudes of the complex amplitudes
    al  = r.*exp(1i*2*pi*rand(K,1));                 % complex amplitudes   
    x   = A*al;                                      % original signal
    Pn  = mean(abs(x(Mcal+1)).^2)*10^(-SNR/10);      % noise power
    eps0 = sqrt(0.5*Pn).*(randn(M,1)+1i*randn(M,1));  % complex Gaussian noise
    y   = x(Mcal+1) + eps0;                           % measured signal
    
    for Bit_dep = 1:length(B_all)
        B = B_all(Bit_dep);
        if B~=inf
                 y_real = [real(y);imag(y)];
                nbins = 2^B;   
                yy_max = 3*sqrt(K/2);
                yy_min = -yy_max; 
                alpha = (yy_max - yy_min)/(nbins);    
                tau_o = 1:nbins-1;
                tau = yy_min +tau_o.*alpha;                  
                yy = floor((y_real-yy_min)/alpha);
                index1 = find(y_real>=yy_max);
                yy(index1) = nbins-1;
                index2 = find(y_real<yy_min);
                yy(index2) = 0;     
                y_pre = yy_min + (yy+0.5)* alpha;
                y_pre_c = y_pre(1:end/2)+1j*y_pre(end/2+1:end);
        else
                 yy = y;
                 y_pre_c = y;
        end

        outVALSE_EP = VALSE_EP(yy, Mcal, 2, x, Iter_max, B, yy_min, alpha, method_EP);
%         outVALSE_EP.K(end)
%         outVALSE_EP.MSE(end)
        MSE_Iter(:,Bit_dep,mc) = outVALSE_EP.MSE;
        outVALSE= VALSE( y_pre_c, Mcal, 2, x, Iter_max,B );
%         outVALSE.K(end)
%         outVALSE.MSE(end)
        MSE_VALSE(:,Bit_dep,mc) = outVALSE.MSE;      
    end
end

MSE_Iter_mean = squeeze(mean(MSE_Iter,3));
MSE_VALSE_mean = squeeze(mean(MSE_VALSE,3));

Iter_max0 = 20;
Iter_index = 1:1:Iter_max0;
%% Figures
alw = 0.75;    % AxesLineWidth
fsz = 18;      % Fontsize
lw = 1.6;      % LineWidth
msz = 8;       % MarkerSize

h1 = figure(1);
set(h1,'position',[100,100,750,500])
set(gca, 'FontSize', fsz, 'LineWidth', alw);
plot(Iter_index,MSE_VALSE_mean(Iter_index,1),'-r^','LineWidth',lw, 'MarkerSize', msz);
hold on
plot(Iter_index,MSE_VALSE_mean(Iter_index,2),'-rv','LineWidth',lw, 'MarkerSize', msz);
hold on
plot(Iter_index,MSE_VALSE_mean(Iter_index,3),'-r<','LineWidth',lw, 'MarkerSize', msz);
hold on
plot(Iter_index,MSE_VALSE_mean(Iter_index,4),'-r>','LineWidth',lw, 'MarkerSize', msz);
plot(Iter_index,MSE_Iter_mean(Iter_index,1),'-.b+','LineWidth',lw, 'MarkerSize', msz);
hold on
plot(Iter_index,MSE_Iter_mean(Iter_index,2),'-.bx','LineWidth',lw, 'MarkerSize', msz);
hold on
plot(Iter_index,MSE_Iter_mean(Iter_index,3),'-.b.','LineWidth',lw, 'MarkerSize', msz+18);
hold on
plot(Iter_index,MSE_Iter_mean(Iter_index,4),'-.b*','LineWidth',lw, 'MarkerSize', msz);
hold on
hhh = legend('1 bit VALSE','2 bit VALSE','3 bit VALSE','\infty bit VALSE','1 bit VALSE-EP','2 bit VALSE-EP','3 bit VALSE-EP','\infty bit VALSE-EP');
set(hhh,'Fontsize',fsz-2);
xlabel('Iteration','Fontsize',fsz);
ylabel('${\rm NMSE}(\hat{\mathbf z})$ (dB)','interpreter','latex','Fontsize',fsz);
% ylabel('$20\log\frac{\|\hat{\mathbf z}-{\mathbf z}\|_2}{\|{\mathbf z}\|_2}$ (dB)','interpreter','latex','Fontsize',fsz);
xlim([1,20]);
% ylim([-24,-4])
set(gca,'xtick',Iter_index);
set(gca,'Fontsize',fsz);

if SNR==0
    pos = [12,-8];
    save('MSE_Itersnr0.mat','MSE_Iter_mean','MSE_Iter_mean','Iter_max')
elseif SNR==20
    pos = [8,-15];
    save('MSE_Itersnr20.mat','MSE_Iter_mean','MSE_Iter_mean','Iter_max')
elseif SNR==40
    pos = [10,-12];
%     save('MSE_Itersnr40.mat','MSE_Iter_mean','MSE_Iter_mean','Iter_max')
end
text(pos(1),pos(2),sprintf('SNR=%d dB',SNR),'Fontsize',fsz)   % 10 dB, pos = (1.5,-24) % 20 dB, pos = (8,-30) % 30 dB, pos = (3,-45)

