% Moment of Velocity (Hilbert Transform)
% Signal's DC componenrt is ignored

clc
clear
close all

% Person 1 (Healthy)
load('f1o02m');
s_1=(val(1,:));
s_1=s_1-mean(s_1);
r_1=floor(length(s_1)/6.25); %desired range
m_1=s_1(1:r_1);
s_1=s_1/(max(m_1)-min(m_1));

t=(0:length(s_1)-1);

m_1=s_1(1:r_1);
s_1=s_1-min(m_1);

s_1=s_1-mean(s_1);  % eliminating DC component


% % performing the HILBERT TRANSFORM
H_1 = hilbert(s_1);
HT_1=imag(H_1);
% % calculating the AMPLITUDE (ENVELOPE)
inst_amp_1 = sqrt([real(H_1).*real(H_1)]+[imag(H_1).*imag(H_1)]);
% % calculating the FREQUENCY
inst_freq_1 = gradient(atan(imag(H_1)./real(H_1)))/(1/(2*pi));
G_1=inst_freq_1(1:r_1);
inst_freq_1=inst_freq_1-(max(G_1)+min(G_1))/2;
G_1=inst_freq_1(1:r_1);
inst_freq_1=10* inst_freq_1/max(G_1);
% %calculating the MOMENT OF VELOCITY
mov_1= inst_freq_1.*(inst_amp_1.*inst_amp_1);


% Patient 2 (St. Petersburg Institute of Cardiological Technics 12-lead 
% Arrhythmia Database. Seventy-five half-hour recordings extracted from 32 
% Holter records from patients undergoing tests for coronary artery disease,
% with reference annotation files containing over 175,000 beat annotations
% in all.)
load('I01m');
s_2=(val(1,:));
s_2=s_2-mean(s_2);
r_2=floor(length(s_2)/6.25);
m_2=s_2(1:r_2);
s_2=s_2/(max(m_2)-min(m_2));

t=(0:length(s_2)-1);

m_2=s_2(1:r_2);
s_2=s_2-min(m_2);
s_2=s_2-mean(s_2);  % eliminating DC component


% % performing the HILBERT TRANSFORM
H_2 = hilbert(s_2);
HT_2=imag(H_2);
% % calculating the AMPLITUDE (ENVELOPE)
inst_amp_2 = sqrt([real(H_2).*real(H_2)]+[imag(H_2).*imag(H_2)]);
% % calculating the FREQUENCY
inst_freq_2 = gradient(atan(imag(H_2)./real(H_2)))/(1/(2*pi));
G_2=inst_freq_2(1:r_2);
inst_freq_2=inst_freq_2-(max(G_2)+min(G_2))/2;
G_2=inst_freq_2(1:r_2);
inst_freq_2=10* inst_freq_2/max(G_2);
% %calculating the MOMENT OF VELOCITY
mov_2= inst_freq_2.*(inst_amp_2.*inst_amp_2);




% Patient 3 (MIT-BIH Arrhythmia Database. This collection of 48 fully 
%annotated half-hour two-lead ECGs is available here in its entirety.)
load('100m');
s_3=(val(1,:));
s_3=s_3-mean(s_1);
r_3=floor(length(s_3)/6.25);
m_3=s_3(1:r_3);
s_3=s_3/(max(m_3)-min(m_3));

t=(0:length(s_3)-1);

m_3=s_3(1:r_3);
s_3=s_3-min(m_3);
s_3=s_3-mean(s_3);  % eliminating DC component

% % performing the HILBERT TRANSFORM
H_3 = hilbert(s_3);
HT_3=imag(H_3);
% % calculating the AMPLITUDE (ENVELOPE)
inst_amp_3 = sqrt([real(H_3).*real(H_3)]+[imag(H_3).*imag(H_3)]);
% % calculating the FREQUENCY
inst_freq_3 = gradient(atan(imag(H_3)./real(H_3)))/(1/(2*pi));
G_3=inst_freq_3(1:r_3);
inst_freq_3=inst_freq_3-(max(G_3)+min(G_3))/2;
G_3=inst_freq_3(1:r_3);
inst_freq_3=10* inst_freq_3/max(G_3);
% %calculating the MOMENT OF VELOCITY
mov_3= inst_freq_3.*(inst_amp_3.*inst_amp_3);


% Patient 4 ([Class 2; core] BIDMC Congestive Heart Failure Database. 
% Long-term ECGs (about 20 hours each) from 15 subjects with severe CHF 
% (NYHA class 3-4).)
load('chf14m');
s_4=(val(1,:));
s_4=s_4-mean(s_4);
r_4=floor(length(s_4)/6.25);
m_4=s_4(1:r_4);
s_4=s_4/(max(m_4)-min(m_4));

t=(0:length(s_4)-1);

m_4=s_4(1:r_4);
s_4=s_4-min(m_4);
s_4=s_4-mean(s_4);  % eliminating DC component


% % performing the HILBERT TRANSFORM
H_4 = hilbert(s_4);
HT_4=imag(H_4);
% % calculating the AMPLITUDE (ENVELOPE)
inst_amp_4 = sqrt([real(H_4).*real(H_4)]+[imag(H_4).*imag(H_4)]);
% % calculating the FREQUENCY
inst_freq_4 = gradient(atan(imag(H_4)./real(H_4)))/(1/(2*pi));
G_4=inst_freq_4(1:r_4);
inst_freq_4=inst_freq_4-(max(G_4)+min(G_4))/2;
G_4=inst_freq_4(1:r_4);
inst_freq_4=10* inst_freq_4/max(G_4);
% %calculating the MOMENT OF VELOCITY
mov_4= inst_freq_4.*(inst_amp_4.*inst_amp_4);



% % Patient 5 STAFF-III Database. The STAFF III database was acquired 
% during 1995–96 at Charleston Area Medical Center (WV, USA) where single 
% prolonged balloon inflation had been introduced to achieve optimal results 
% of percutaneous transluminal coronary angiography (PTCA) procedures, 
% replacing the typical series of brief inflations. The database consists 
% of standard 12-lead ECG recordings from 104 patients.

load('106em');
s_5=(val(6,:));
s_5=s_5-mean(s_5);
r_5=floor(length(s_5)/6.25);
m_5=s_5(1:r_5);
s_5=s_5/(max(m_5)-min(m_5));

t=(0:length(s_5)-1);

m_5=s_5(1:r_5);
s_5=s_5-min(m_5);
s_5=s_5-mean(s_5);  % eliminating DC component

% % performing the HILBERT TRANSFORM
H_5 = hilbert(s_5);
HT_5=imag(H_5);
% % calculating the AMPLITUDE (ENVELOPE)
inst_amp_5 = sqrt([real(H_5).*real(H_5)]+[imag(H_5).*imag(H_5)]);
% % calculating the FREQUENCY
inst_freq_5 = gradient(atan(imag(H_5)./real(H_5)))/(1/(2*pi));
G_5=inst_freq_5(1:r_5);
inst_freq_5=inst_freq_5-(max(G_5)+min(G_5))/2;
G_5=inst_freq_5(1:r_5);
inst_freq_5=10* inst_freq_5/max(G_5);
% %calculating the MOMENT OF VELOCITY
mov_5= inst_freq_5.*(inst_amp_5.*inst_amp_5);




% % Patient 6 PATIENT DATA:
% 
% PATIENT DATA:
% AGE/SEX:
%   41 year old male
% DIAGNOSIS/SURGERY:
%   Congestive heart failure
% PERTINENT HISTORY:
%   Severe mitral stenosis
%   Moderate mitral regurgitation
%   Aortic stenosis
%   s/p AVR and MVR 10 years earlier


load('mgh156m');
s_6=(val(1,1701:3500));
s_6=s_6-mean(s_6);

r_6=floor(length(s_6)/3.125); %time interval [1.6 sec]
m_6=s_6(1:r_6);
s_6=s_6/(max(m_6)-min(m_6));

t=(0:length(s_6)-1);

m_6=s_6(1:r_6);
s_6=s_6-min(m_6);
s_6=s_6-mean(s_6);  % eliminating DC component

% % performing the HILBERT TRANSFORM
H_6 = hilbert(s_6);
HT_6=imag(H_6);
% % calculating the AMPLITUDE (ENVELOPE)
inst_amp_6 = sqrt([real(H_6).*real(H_6)]+[imag(H_6).*imag(H_6)]);
% % calculating the FREQUENCY
inst_freq_6 = gradient(atan(imag(H_6)./real(H_6)))/(1/(2*pi));
G_6=inst_freq_6(1:r_6);
inst_freq_6=inst_freq_6-(max(G_6)+min(G_6))/2;
G_6=inst_freq_6(1:r_6);
inst_freq_6=10* inst_freq_6/max(G_6);
% %calculating the MOMENT OF VELOCITY
mov_6= inst_freq_6.*(inst_amp_6.*inst_amp_6);


% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%  FIGURE 1  %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(1);
% subplot(2,3,1), plot(t(1,1:r_1),s_1(1,1:r_1)),title('(a) Healthy Heartbeat', 'FontSize', 16)
% xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
% axis([0 r_1 -1 1])
% 
% subplot(2,3,2), plot(t(1,1:r_1),inst_amp_1(1,1:r_1)),title('(b) Instantaneous Amplitude', 'FontSize', 16)
% xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% % ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_1 0 1])
% 
% 
% subplot(2,3,3), plot(t(1,1:r_1),inst_freq_1(1,1:r_1)),title('(c) Instantaneous Frequency', 'FontSize', 16)
% xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Frequency in Hz', 'FontSize', 16) % y-axis label
% axis([0 r_1 -10 10])
% 
% subplot(2,3,4), plot(t(1,1:r_1),mov_1(1,1:r_1)),title('(d) Moment of Velocity', 'FontSize', 16)
% xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('V^{2}/s', 'FontSize', 16) % y-axis label
% xlim([0 r_1])
% ylim([-6 6])
% 
% 
% subplot(2,3,5), plot(s_1(1,1:r_1),HT_1(1,1:r_1)),title('(e) Signal Vs Hilbert T', 'FontSize', 16)
% hold on
% plot(0,0,'+');
% hold off
% 
% % xlabel('Time in seconds') % x-axis label
% % % ylabel('Voltage Amplitude') % y-axis label
% % axis([0 r_1 0 1])
% 
% 
% 
% subplot(2,3,6),
% 
% plot3(s_1(1,1:r_1),HT_1(1,1:r_1),t(1,1:r_1)),title('(f) Signal Vs Hilbert T (3D)', 'FontSize', 16)
% xlabel('Signal', 'FontSize', 16) % x-axis label
% ylabel('Hilbert Transform', 'FontSize', 16) % y-axis label
% zlabel('Time in seconds', 'FontSize', 16) % z-axis label
% zticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
% zticklabels({'0','0.5','1','1.5','2'})
% hold on
% plot3(0,0,t(1,1:r_1),'.')
% hold off
% x1 = 0;
% y1 = 0;
% z1 = r_1;
% txt = '\leftarrow-----x=y=0';
% text(x1,y1,z1,txt)
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%  FIGURE 2  %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(2);
% subplot(2,3,1), plot(t(1,1:r_2),s_2(1,1:r_2)),title('(a) Coronary Artery Disease Heartbeat', 'FontSize', 16)
% xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
% axis([0 r_2 -1 1])
% 
% 
% subplot(2,3,2), plot(t(1,1:r_2),inst_amp_2(1,1:r_2)),title('(b) Instantaneous Amplitude', 'FontSize', 16)
% xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% % ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_2 0 1])
% 
% subplot(2,3,3), plot(t(1,1:r_2),inst_freq_2(1,1:r_2)),title('(c) Instantaneous Frequency', 'FontSize', 16)
% xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Frequency in Hz', 'FontSize', 16) % y-axis label
% axis([0 r_2 -10 10])
% 
% subplot(2,3,4), plot(t(1,1:r_2),mov_2(1,1:r_2)),title('(d) Moment of Velocity', 'FontSize', 16)
% xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('V^{2}/s', 'FontSize', 16) % y-axis label
% xlim([0 r_2])
% ylim([-6 6])
% 
% 
% subplot(2,3,5), plot(s_2(1,1:r_2),HT_2(1,1:r_2)),title('(e) Signal Vs Hilbert T', 'FontSize', 16)
% hold on
% plot(0,0,'+');
% hold off
% 
% % xlabel('Time in seconds') % x-axis label
% % % ylabel('Voltage Amplitude') % y-axis label
% % axis([0 r_1 0 1])
% 
% 
% 
% subplot(2,3,6),
% 
% plot3(s_2(1,1:r_2),HT_2(1,1:r_2),t(1,1:r_2)),title('(f) Signal Vs Hilbert T (3D)', 'FontSize', 16)
% xlabel('Signal', 'FontSize', 16) % x-axis label
% ylabel('Hilbert Transform', 'FontSize', 16) % y-axis label
% zlabel('Time in seconds', 'FontSize', 16) % z-axis label
% zticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
% zticklabels({'0','0.5','1','1.5','2'})
% hold on
% plot3(0,0,t(1,1:r_2),'.')
% hold off
% x2 = 0;
% y2 = 0;
% z2 = r_2;
% txt = '\leftarrow-----x=y=0';
% text(x2,y2,z2,txt)
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%  FIGURE 3  %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(3);
% subplot(2,3,1), plot(t(1,1:r_3),s_3(1,1:r_3)),title('(a) Arrhythmia Heartbeat', 'FontSize', 16)
% xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_3 -1 1])
% 
% subplot(2,3,2), plot(t(1,1:r_3),inst_amp_3(1,1:r_3)),title('(b) Instantaneous Amplitude', 'FontSize', 16)
% xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% % ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_3 0 1])
% 
% subplot(2,3,3), plot(t(1,1:r_3),inst_freq_3(1,1:r_3)),title('(c) Instantaneous Frequency', 'FontSize', 16)
% xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% ylabel('Frequency in Hz') % y-axis label
% axis([0 r_3 -10 10])
% 
% subplot(2,3,4), plot(t(1,1:r_3),mov_3(1,1:r_3)),title('(d) Moment of Velocity', 'FontSize', 16)
% xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% ylabel('(V^{2}/s)') % y-axis label
% xlim([0 r_3])
% ylim([-6 6])
% 
% subplot(2,3,5), plot(s_3(1,1:r_3),HT_3(1,1:r_3)),title('(e) Signal Vs Hilbert T', 'FontSize', 16)
% hold on
% plot(0,0,'+');
% hold off
% 
% % xlabel('Time in seconds') % x-axis label
% % % ylabel('Voltage Amplitude') % y-axis label
% % axis([0 r_1 0 1])
% 
% 
% 
% subplot(2,3,6),
% 
% plot3(s_3(1,1:r_3),HT_3(1,1:r_3),t(1,1:r_3)),title('(f) Signal Vs Hilbert T (3D)', 'FontSize', 16)
% xlabel('Signal') % x-axis label
% ylabel('Hilbert Transform') % y-axis label
% zlabel('Time in seconds') % z-axis label
% zticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
% zticklabels({'0','0.5','1','1.5','2'})
% hold on
% plot3(0,0,t(1,1:r_3),'.')
% hold off
% x3 = 0;
% y3 = 0;
% z3 = r_3;
% txt = '\leftarrow-----x=y=0';
% text(x3,y3,z3,txt)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%  FIGURE 4  %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(4);
% subplot(2,3,1), plot(t(1,1:r_4),s_4(1,1:r_4)),title('(a) Congestive Heart Failure Signal', 'FontSize', 16)
% xticks([0 r_4/4 r_4/2 r_4*3/4 r_4])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_4 -1 1])
% 
% 
% subplot(2,3,2), plot(t(1,1:r_4),inst_amp_4(1,1:r_4)),title('(b) Instantaneous Amplitude', 'FontSize', 16)
% xticks([0 r_4/4 r_4/2 r_4*3/4 r_4])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% % ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_4 0 1])
% 
% subplot(2,3,3), plot(t(1,1:r_4),inst_freq_4(1,1:r_4)),title('(c) Instantaneous Frequency', 'FontSize', 16)
% xticks([0 r_4/4 r_4/2 r_4*3/4 r_4])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% ylabel('Frequency in HZ') % y-axis label
% axis([0 r_4 -10 10])
% 
% subplot(2,3,4), plot(t(1,1:r_4),mov_4(1,1:r_4)),title('(d) Moment of Velocity', 'FontSize', 16)
% xticks([0 r_4/4 r_4/2 r_4*3/4 r_4])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% ylabel('(V^{2}/s)') % y-axis label
% xlim([0 r_4])
% ylim([-6 6])
% 
% subplot(2,3,5), plot(s_4(1,1:r_4),HT_4(1,1:r_4)),title('(e) Signal Vs Hilbert T', 'FontSize', 16)
% hold on
% plot(0,0,'+');
% hold off
% 
% % xlabel('Time in seconds') % x-axis label
% % % ylabel('Voltage Amplitude') % y-axis label
% % axis([0 r_1 0 1])
% 
% 
% 
% subplot(2,3,6),
% 
% plot3(s_4(1,1:r_4),HT_4(1,1:r_4),t(1,1:r_4)),title('(f) Signal Vs Hilbert T (3D)', 'FontSize', 16)
% xlabel('Signal') % x-axis label
% ylabel('Hilbert Transform') % y-axis label
% zlabel('Time in seconds') % z-axis label
% zticks([0 r_4/4 r_4/2 r_4*3/4 r_4])
% zticklabels({'0','0.5','1','1.5','2'})
% hold on
% plot3(0,0,t(1,1:r_4),'.')
% hold off
% x4 = 0;
% y4 = 0;
% z4 = r_4;
% txt = '\leftarrow-----x=y=0';
% text(x4,y4,z4,txt)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%  FIGURE 5  %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(5);
% subplot(2,3,1), plot(t(1,1:r_5),s_5(1,1:r_5)),title('(a) Ventricular Tachycardia Signal', 'FontSize', 16)
% xticks([0 r_5/4 r_5/2 r_5*3/4 r_5])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_5 -1 1])
% 
% 
% subplot(2,3,2), plot(t(1,1:r_5),inst_amp_5(1,1:r_5)),title('(b) Instantaneous Amplitude', 'FontSize', 16)
% xticks([0 r_5/4 r_5/2 r_5*3/4 r_5])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% % ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_5 0 1])
% 
% subplot(2,3,3), plot(t(1,1:r_5),inst_freq_5(1,1:r_5)),title('(c) Instantaneous Frequency', 'FontSize', 16)
% xticks([0 r_5/4 r_5/2 r_5*3/4 r_5])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% ylabel('Frequency in Hz') % y-axis label
% axis([0 r_5 -10 10])
% 
% subplot(2,3,4), plot(t(1,1:r_5),mov_5(1,1:r_5)),title('(d) Moment of Velocity', 'FontSize', 16)
% xticks([0 r_5/4 r_5/2 r_5*3/4 r_5])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% ylabel('(V^{2}/s)') % y-axis label
% xlim([0 r_5])
% ylim([0 10])
% ylim([-6 6])
% 
% subplot(2,3,5), plot(s_5(1,1:r_5),HT_5(1,1:r_5)),title('(e) Signal Vs Hilbert T', 'FontSize', 16)
% hold on
% plot(0,0,'+');
% hold off
% 
% % xlabel('Time in seconds') % x-axis label
% % % ylabel('Voltage Amplitude') % y-axis label
% % axis([0 r_1 0 1])
% 
% 
% 
% subplot(2,3,6),
% 
% plot3(s_5(1,1:r_5),HT_5(1,1:r_5),t(1,1:r_5)),title('(f) Signal Vs Hilbert T (3D)', 'FontSize', 16)
% xlabel('Signal') % x-axis label
% ylabel('Hilbert Transform') % y-axis label
% zlabel('Time in seconds') % z-axis label
% zticks([0 r_5/4 r_5/2 r_5*3/4 r_5])
% zticklabels({'0','0.5','1','1.5','2'})
% hold on
% plot3(0,0,t(1,1:r_5),'.')
% hold off
% x5 = 0;
% y5 = 0;
% z5 = r_5;
% txt = '\leftarrow-----x=y=0';
% text(x5,y5,z5,txt)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%  FIGURE 6  %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(6);
% subplot(2,3,1), plot(t,s_6),title('(a) Mitral Valve Regurgitation Signal', 'FontSize', 16)
% xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_6 -1 1])
% 
% 
% subplot(2,3,2), plot(t,inst_amp_6),title('(b) Instantaneous Amplitude', 'FontSize', 16)
% xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% % ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_6 0 1])
% 
% subplot(2,3,3), plot(t,inst_freq_6),title('(c) Instantaneous Frequency', 'FontSize', 16)
% xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% ylabel('Frequency in Hz') % y-axis label
% axis([0 r_6 -10 10])
% 
% subplot(2,3,4), plot(t,mov_6),title('(d) Moment of Velocity', 'FontSize', 16)
% xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds') % x-axis label
% ylabel('(V^{2}/s)') % y-axis label
% xlim([1 r_6])
% ylim([-6 6])
% 
% subplot(2,3,5), plot(s_6(1,1:r_6),HT_6(1,1:r_6)),title('(e) Signal Vs Hilbert T', 'FontSize', 16)
% % xlim([1:r_6])
% % ylim([-6 6])
% hold on
% plot(0,0,'+');
% hold off
% 
% 
% % xlabel('Time in seconds') % x-axis label
% % % ylabel('Voltage Amplitude') % y-axis label
% % axis([0 r_1 0 1])
% 
% 
% 
% subplot(2,3,6),
% 
% plot3(s_6(1,1:r_6),HT_6(1,1:r_6),t(1,1:r_6)),title('(f) Signal Vs Hilbert T (3D)', 'FontSize', 16)
% xlabel('Signal') % x-axis label
% ylabel('Hilbert Transform') % y-axis label
% zlabel('Time in seconds') % z-axis label
% zticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
% zticklabels({'0','0.5','1','1.5','2'})
% hold on
% plot3(0,0,t(1,1:r_6),'.')
% hold off
% x6 = 0;
% y6 = 0;
% z6 = r_6;
% txt = '\leftarrow-----x=y=0';
% text(x6,y6,z6,txt)




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%  Figures 7  %%%%%%%%%%%%  Signals  %%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(7);
% 
% subplot(2,3,1), plot(t(1,1:r_1),s_1(1,1:r_1)),title('(a) Healthy Heartbeat', 'FontSize', 16)
% xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
% axis([0 r_1 -1 1])
% 
% 
% subplot(2,3,2), plot(t(1,1:r_2),s_2(1,1:r_2)),title('(b) Coronary Artery Disease', 'FontSize', 16)
% xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
% axis([0 r_2 -1 1])
% 
% 
% subplot(2,3,3), plot(t(1,1:r_3),s_3(1,1:r_3)),title('(c) Arrhythmia', 'FontSize', 16)
% xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
% axis([0 r_3 -1 1])
% 
% 
% subplot(2,3,4), plot(t(1,1:r_4),s_4(1,1:r_4)),title('(d) Congestive Heart Failure', 'FontSize', 16)
% xticks([0 r_4/4 r_4/2 r_4*3/4 r_4])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
% axis([0 r_4 -1 1])
% 
% 
% subplot(2,3,5), plot(t(1,1:r_5),s_5(1,1:r_5)),title('(e) Ventricular Tachycardia', 'FontSize', 16)
% xticks([0 r_5/4 r_5/2 r_5*3/4 r_5])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
% axis([0 r_5 -1 1])
% 
% 
% subplot(2,3,6), plot(t,s_6),title('(f) Mitral Valve Regurgitation', 'FontSize', 16)
% xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
% axis([0 r_6 -1 1])
% 
% 
% 
% h=gcf;
% set(h,'PaperPositionMode','auto');         
% set(h,'PaperOrientation','landscape');
% set(h,'Position',[50 50 1200 800]);
% print(gcf, '-dpdf', 'Signals.pdf')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%  Figures 8  %%%%%%%%%%%%  IA  %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(8);
% 
% subplot(2,3,1), plot(t(1,1:r_1),inst_amp_1(1,1:r_1)),title('(a)', 'FontSize', 16)
% xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Amplitude', 'FontSize', 16) % y-axis label
% % ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_1 0 1])
% 
% 
% subplot(2,3,2), plot(t(1,1:r_2),inst_amp_2(1,1:r_2)),title('(b)', 'FontSize', 16)
% xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Amplitude', 'FontSize', 16) % y-axis label
% % ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_2 0 1])
% 
% 
% subplot(2,3,3), plot(t(1,1:r_3),inst_amp_3(1,1:r_3)),title('(c)', 'FontSize', 16)
% xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Amplitude', 'FontSize', 16) % y-axis label
% % ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_3 0 1])
% 
% 
% subplot(2,3,4), plot(t(1,1:r_4),inst_amp_4(1,1:r_4)),title('(d)', 'FontSize', 16)
% xticks([0 r_4/4 r_4/2 r_4*3/4 r_4])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Amplitude', 'FontSize', 16) % y-axis label
% % ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_4 0 1])
% 
% subplot(2,3,5), plot(t(1,1:r_5),inst_amp_5(1,1:r_5)),title('(e)', 'FontSize', 16)
% xticks([0 r_5/4 r_5/2 r_5*3/4 r_5])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Amplitude', 'FontSize', 16) % y-axis label
% % ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_5 0 1])
% 
% subplot(2,3,6), plot(t,inst_amp_6),title('(f)', 'FontSize', 16)
% xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Amplitude', 'FontSize', 16) % y-axis label
% % ylabel('Voltage Amplitude') % y-axis label
% axis([0 r_6 0 1])
% 
% h=gcf;
% set(h,'PaperPositionMode','auto');         
% set(h,'PaperOrientation','landscape');
% set(h,'Position',[50 50 1200 800]);
% print(gcf, '-dpdf', 'ia.pdf')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%  Figures 9  %%%%%%%%%%%%  IF  %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(9);
% 
% subplot(2,3,1), plot(t(1,1:r_1),inst_freq_1(1,1:r_1)),title('(a)', 'FontSize', 16)
% xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Frequency in Hz', 'FontSize', 16) % y-axis label
% axis([0 r_1 -10 10])
% 
% 
% subplot(2,3,2), plot(t(1,1:r_2),inst_freq_2(1,1:r_2)),title('(b)', 'FontSize', 16)
% xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Frequency in Hz', 'FontSize', 16) % y-axis label
% axis([0 r_2 -10 10])
% 
% subplot(2,3,3), plot(t(1,1:r_3),inst_freq_3(1,1:r_3)),title('(c)', 'FontSize', 16)
% xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Frequency in Hz', 'FontSize', 16) % y-axis label
% axis([0 r_3 -10 10])
% 
% 
% 
% subplot(2,3,4), plot(t(1,1:r_4),inst_freq_4(1,1:r_4)),title('(d)', 'FontSize', 16)
% xticks([0 r_4/4 r_4/2 r_4*3/4 r_4])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Frequency in HZ', 'FontSize', 16) % y-axis label
% axis([0 r_4 -10 10])
% 
% subplot(2,3,5), plot(t(1,1:r_5),inst_freq_5(1,1:r_5)),title('(e)', 'FontSize', 16)
% xticks([0 r_5/4 r_5/2 r_5*3/4 r_5])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Frequency in Hz', 'FontSize', 16) % y-axis label
% axis([0 r_5 -10 10])
% 
% subplot(2,3,6), plot(t,inst_freq_6),title('(f)', 'FontSize', 16)
% xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('Frequency in Hz', 'FontSize', 16) % y-axis label
% axis([0 r_6 -10 10])
% 
% h=gcf;
% set(h,'PaperPositionMode','auto');         
% set(h,'PaperOrientation','landscape');
% set(h,'Position',[50 50 1200 800]);
% print(gcf, '-dpdf', 'if.pdf')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%  Figures 10  %%%%%%%%%%%%  HT  %%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(10);
% 
% subplot(2,3,1),
% 
% plot3(s_1(1,1:r_1),HT_1(1,1:r_1),t(1,1:r_1),'LineWidth',1.5),title('(a)', 'FontSize', 16)
% xlabel('Real', 'FontSize', 16) % x-axis label
% ylabel('Imaginary', 'FontSize', 16) % y-axis label
% zlabel('Time (S)', 'FontSize', 16) % z-axis label
% zticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
% zticklabels({'0','0.5','1','1.5','2'})
% hold on
% plot3(0,0,t(1,1:r_1),'.')
% hold off
% x1 = 0;
% y1 = 0;
% z1 = r_1;
% txt = '\leftarrow-----x=y=0';
% text(x1,y1,z1,txt)
% 
% 
% 
% subplot(2,3,2),
% 
% plot3(s_2(1,1:r_2),HT_2(1,1:r_2),t(1,1:r_2),'LineWidth',1.5),title('(b)', 'FontSize', 16)
% xlabel('Real', 'FontSize', 16) % x-axis label
% ylabel('Imaginary', 'FontSize', 16) % y-axis label
% zlabel('Time (S)', 'FontSize', 16) % z-axis label
% zticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
% zticklabels({'0','0.5','1','1.5','2'})
% hold on
% plot3(0,0,t(1,1:r_2),'.')
% hold off
% x2 = 0;
% y2 = 0;
% z2 = r_2;
% txt = '\leftarrow-----x=y=0';
% text(x2,y2,z2,txt)
% 
% 
% subplot(2,3,3),
% 
% plot3(s_3(1,1:r_3),HT_3(1,1:r_3),t(1,1:r_3),'LineWidth',1.5),title('(c)', 'FontSize', 16)
% xlabel('Real', 'FontSize', 16) % x-axis label
% ylabel('Imaginary', 'FontSize', 16) % y-axis label
% zlabel('Time (S)', 'FontSize', 16) % z-axis label
% zticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
% zticklabels({'0','0.5','1','1.5','2'})
% hold on
% plot3(0,0,t(1,1:r_3),'.')
% hold off
% x3 = 0;
% y3 = 0;
% z3 = r_3;
% txt = '\leftarrow-----x=y=0';
% text(x3,y3,z3,txt)
% 
% 
% 
% subplot(2,3,4),
% 
% plot3(s_4(1,1:r_4),HT_4(1,1:r_4),t(1,1:r_4),'LineWidth',1.5),title('(d)', 'FontSize', 16)
% xlabel('Real', 'FontSize', 16) % x-axis label
% ylabel('Imaginary', 'FontSize', 16) % y-axis label
% zlabel('Time (S)', 'FontSize', 16) % z-axis label
% zticks([0 r_4/4 r_4/2 r_4*3/4 r_4])
% zticklabels({'0','0.5','1','1.5','2'})
% hold on
% plot3(0,0,t(1,1:r_4),'.')
% hold off
% x4 = 0;
% y4 = 0;
% z4 = r_4;
% txt = '\leftarrow-----x=y=0';
% text(x4,y4,z4,txt)
% 
% 
% 
% subplot(2,3,5),
% 
% plot3(s_5(1,1:r_5),HT_5(1,1:r_5),t(1,1:r_5),'LineWidth',1.5),title('(e)', 'FontSize', 16)
% xlabel('Real', 'FontSize', 16) % x-axis label
% ylabel('Imaginary', 'FontSize', 16) % y-axis label
% zlabel('Time (S)', 'FontSize', 16) % z-axis label
% zticks([0 r_5/4 r_5/2 r_5*3/4 r_5])
% zticklabels({'0','0.5','1','1.5','2'})
% hold on
% plot3(0,0,t(1,1:r_5),'.')
% hold off
% x5 = 0;
% y5 = 0;
% z5 = r_5;
% txt = '\leftarrow-----x=y=0';
% text(x5,y5,z5,txt)
% 
% 
% 
% subplot(2,3,6),
% 
% plot3(s_6(1,1:r_6),HT_6(1,1:r_6),t(1,1:r_6),'LineWidth',1.5),title('(f)', 'FontSize', 16)
% xlabel('Real', 'FontSize', 16) % x-axis label
% ylabel('Imaginary', 'FontSize', 16) % y-axis label
% zlabel('Time (S)', 'FontSize', 16) % z-axis label
% zticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
% zticklabels({'0','0.5','1','1.5','2'})
% hold on
% plot3(0,0,t(1,1:r_6),'.')
% hold off
% x6 = 0;
% y6 = 0;
% z6 = r_6;
% txt = '\leftarrow-----x=y=0';
% text(x6,y6,z6,txt)
% 
% 
% h=gcf;
% set(h,'PaperPositionMode','auto');         
% set(h,'PaperOrientation','landscape');
% set(h,'Position',[50 50 1200 800]);
% print(gcf, '-dpdf', 'HT.pdf')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%  Figures 11  %%%%%%%%%%%% MoV  %%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(11);
% 
% subplot(2,3,1), plot(t(1,1:r_1),mov_1(1,1:r_1)),title('(a)', 'FontSize', 16)
% xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('V^{2}/s', 'FontSize', 16) % y-axis label
% xlim([0 r_1])
% ylim([-6 6])
% 
% subplot(2,3,2), plot(t(1,1:r_2),mov_2(1,1:r_2)),title('(b)', 'FontSize', 16)
% xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('V^{2}/s', 'FontSize', 16) % y-axis label
% xlim([0 r_2])
% ylim([-6 6])
% 
% subplot(2,3,3), plot(t(1,1:r_3),mov_3(1,1:r_3)),title('(c)', 'FontSize', 16)
% xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('V^{2}/s', 'FontSize', 16) % y-axis label
% xlim([0 r_3])
% ylim([-6 6])
% 
% subplot(2,3,4), plot(t(1,1:r_4),mov_4(1,1:r_4)),title('(d)', 'FontSize', 16)
% xticks([0 r_4/4 r_4/2 r_4*3/4 r_4])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('V^{2}/s', 'FontSize', 16) % y-axis label
% xlim([0 r_4])
% ylim([-6 6])
% 
% subplot(2,3,5), plot(t(1,1:r_5),mov_5(1,1:r_5)),title('(e)', 'FontSize', 16)
% xticks([0 r_5/4 r_5/2 r_5*3/4 r_5])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('V^{2}/s', 'FontSize', 16) % y-axis label
% xlim([0 r_5])
% ylim([0 10])
% ylim([-6 6])
% 
% subplot(2,3,6), plot(t,mov_6),title('(f)', 'FontSize', 16)
% xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
% xticklabels({'0','0.5','1','1.5','2'})
% xlabel('Time in seconds', 'FontSize', 16) % x-axis label
% ylabel('V^{2}/s', 'FontSize', 16) % y-axis label
% xlim([1 r_6])
% ylim([-6 6])
% 
% h=gcf;
% set(h,'PaperPositionMode','auto');         
% set(h,'PaperOrientation','landscape');
% set(h,'Position',[50 50 1200 800]);
% print(gcf, '-dpdf', 'mov.pdf')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Figures 12  %%%%%%%%%%%%  Signals & IF  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(12);

fig = figure(12);
left_color = [0.85 0.325 0.0980];
right_color = [0 0.4470 0.741];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);




subplot(2,2,1), 

yyaxis left
plot(t(1,1:r_1),s_1(1,1:r_1),'LineWidth',1.5),title('(a)', 'FontSize', 16)
xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_1 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t(1,1:r_1),inst_freq_1(1,1:r_1),'LineWidth',1),title('(a)', 'FontSize', 16)
xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Instantaneous Frequency (Hz)', 'FontSize', 16) % y-axis label
axis([0 r_1 -10 10])
ylim([-30 10])
yticks([-10 -5 0 5 10])




subplot(2,2,2), 

yyaxis left
plot(t(1,1:r_2),s_2(1,1:r_2),'LineWidth',1.5),title('(b)', 'FontSize', 16)
xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_2 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t(1,1:r_2),inst_freq_2(1,1:r_2),'LineWidth',1),title('(b)', 'FontSize', 16)
xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Instantaneous Frequency (Hz)', 'FontSize', 16) % y-axis label
axis([0 r_2 -10 10])
ylim([-30 10])
yticks([-10 -5 0 5 10])





subplot(2,2,3), 
yyaxis left
plot(t(1,1:r_3),s_3(1,1:r_3),'LineWidth',1.5),title('(c)', 'FontSize', 16)
xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_3 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t(1,1:r_3),inst_freq_3(1,1:r_3),'LineWidth',1),title('(c)', 'FontSize', 16)
xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Instantaneous Frequency (Hz)', 'FontSize', 16) % y-axis label
axis([0 r_3 -10 10])
ylim([-30 10])
yticks([-10 -5 0 5 10])






subplot(2,2,4),
yyaxis left
plot(t,s_6,'LineWidth',1.5),title('(d)', 'FontSize', 16)
xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_6 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t,inst_freq_6,'LineWidth',1),title('(d)', 'FontSize', 16)
xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Instantaneous Frequency (Hz)', 'FontSize', 16) % y-axis label
axis([0 r_6 -10 10])
ylim([-30 10])
yticks([-10 -5 0 5 10])




h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'Position',[50 50 1200 800]);
print(gcf, '-dpdf', 'SigVsIF.pdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Figures 12  %%%%%%%%%%%%  Signals & IF  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(13);

fig = figure(13);
left_color = [0.85 0.325 0.0980];
right_color = [0 0.4470 0.741];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);




subplot(2,2,1), 

yyaxis left
plot(t(1,1:r_1),s_1(1,1:r_1),'LineWidth',1.5),title('(a)', 'FontSize', 16)
xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_1 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t(1,1:r_1),mov_1(1,1:r_1),'LineWidth',1),title('(a)', 'FontSize', 16)
xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Moment of Velocity', 'FontSize', 16) % y-axis label
axis([0 r_1 -10 10])
ylim([-25 8])
yticks([-5 -2.5 0 2.5 5])




subplot(2,2,2), 

yyaxis left
plot(t(1,1:r_2),s_2(1,1:r_2),'LineWidth',1.5),title('(b)', 'FontSize', 16)
xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_2 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t(1,1:r_2),mov_2(1,1:r_2),'LineWidth',1),title('(b)', 'FontSize', 16)
xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Moment of Velocity', 'FontSize', 16) % y-axis label
axis([0 r_2 -10 10])
ylim([-25 8])
yticks([-5 -2.5 0 2.5 5])





subplot(2,2,3), 
yyaxis left
plot(t(1,1:r_3),s_3(1,1:r_3),'LineWidth',1.5),title('(c)', 'FontSize', 16)
xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_3 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t(1,1:r_3),mov_3(1,1:r_3),'LineWidth',1),title('(c)', 'FontSize', 16)
xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Moment of Velocity', 'FontSize', 16) % y-axis label
axis([0 r_3 -10 10])
ylim([-25 8])
yticks([-5 -2.5 0 2.5 5])






subplot(2,2,4),
yyaxis left
plot(t,s_6,'LineWidth',1.5),title('(d)', 'FontSize', 16)
xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_6 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t,mov_6,'LineWidth',1),title('(d)', 'FontSize', 16)
xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Moment of Velocity', 'FontSize', 16) % y-axis label
axis([0 r_6 -10 10])
ylim([-25 8])
yticks([-5 -2.5 0 2.5 5])




h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'Position',[50 50 1200 800]);
print(gcf, '-dpdf', 'SigVsMoV.pdf')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%         Noise        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Person 1 (Healthy)
s_1=awgn(s_1,20) ;  % adding some white Gussian noise
r_1=floor(length(s_1)/6.25); %desired range
m_1=s_1(1:r_1);
s_1=s_1/(max(m_1)-min(m_1));

t=(0:length(s_1)-1);

m_1=s_1(1:r_1);
s_1=s_1-min(m_1);

s_1=s_1-mean(s_1);  % eliminating DC component


% % performing the HILBERT TRANSFORM
H_1 = hilbert(s_1);
HT_1=imag(H_1);
% % calculating the AMPLITUDE (ENVELOPE)
inst_amp_1 = sqrt([real(H_1).*real(H_1)]+[imag(H_1).*imag(H_1)]);
% % calculating the FREQUENCY
inst_freq_1 = gradient(atan(imag(H_1)./real(H_1)))/(1/(2*pi));
G_1=inst_freq_1(1:r_1);
inst_freq_1=inst_freq_1-(max(G_1)+min(G_1))/2;
G_1=inst_freq_1(1:r_1);
inst_freq_1=10* inst_freq_1/max(G_1);
% %calculating the MOMENT OF VELOCITY
mov_1= inst_freq_1.*(inst_amp_1.*inst_amp_1);


% Patient 2 (St. Petersburg Institute of Cardiological Technics 12-lead 
% Arrhythmia Database. Seventy-five half-hour recordings extracted from 32 
% Holter records from patients undergoing tests for coronary artery disease,
% with reference annotation files containing over 175,000 beat annotations
% in all.)

s_2=awgn(s_2,20) ;  % adding some white Gussian noise
r_2=floor(length(s_2)/6.25);
m_2=s_2(1:r_2);
s_2=s_2/(max(m_2)-min(m_2));

t=(0:length(s_2)-1);

m_2=s_2(1:r_2);
s_2=s_2-min(m_2);
s_2=s_2-mean(s_2);  % eliminating DC component


% % performing the HILBERT TRANSFORM
H_2 = hilbert(s_2);
HT_2=imag(H_2);
% % calculating the AMPLITUDE (ENVELOPE)
inst_amp_2 = sqrt([real(H_2).*real(H_2)]+[imag(H_2).*imag(H_2)]);
% % calculating the FREQUENCY
inst_freq_2 = gradient(atan(imag(H_2)./real(H_2)))/(1/(2*pi));
G_2=inst_freq_2(1:r_2);
inst_freq_2=inst_freq_2-(max(G_2)+min(G_2))/2;
G_2=inst_freq_2(1:r_2);
inst_freq_2=10* inst_freq_2/max(G_2);
% %calculating the MOMENT OF VELOCITY
mov_2= inst_freq_2.*(inst_amp_2.*inst_amp_2);




% Patient 3 (MIT-BIH Arrhythmia Database. This collection of 48 fully 
%annotated half-hour two-lead ECGs is available here in its entirety.
s_3=awgn(s_3,20) ;  % adding some white Gussian noise
r_3=floor(length(s_3)/6.25);
m_3=s_3(1:r_3);
s_3=s_3/(max(m_3)-min(m_3));

t=(0:length(s_3)-1);

m_3=s_3(1:r_3);
s_3=s_3-min(m_3);
s_3=s_3-mean(s_3);  % eliminating DC component

% % performing the HILBERT TRANSFORM
H_3 = hilbert(s_3);
HT_3=imag(H_3);
% % calculating the AMPLITUDE (ENVELOPE)
inst_amp_3 = sqrt([real(H_3).*real(H_3)]+[imag(H_3).*imag(H_3)]);
% % calculating the FREQUENCY
inst_freq_3 = gradient(atan(imag(H_3)./real(H_3)))/(1/(2*pi));
G_3=inst_freq_3(1:r_3);
inst_freq_3=inst_freq_3-(max(G_3)+min(G_3))/2;
G_3=inst_freq_3(1:r_3);
inst_freq_3=10* inst_freq_3/max(G_3);
% %calculating the MOMENT OF VELOCITY
mov_3= inst_freq_3.*(inst_amp_3.*inst_amp_3);


% Patient 4 ([Class 2; core] BIDMC Congestive Heart Failure Database. 
% Long-term ECGs (about 20 hours each) from 15 subjects with severe CHF 
% (NYHA class 3-4).)
s_4=awgn(s_4,20) ;  % adding some white Gussian noise
r_4=floor(length(s_4)/6.25);
m_4=s_4(1:r_4);
s_4=s_4/(max(m_4)-min(m_4));

t=(0:length(s_4)-1);

m_4=s_4(1:r_4);
s_4=s_4-min(m_4);
s_4=s_4-mean(s_4);  % eliminating DC component


% % performing the HILBERT TRANSFORM
H_4 = hilbert(s_4);
HT_4=imag(H_4);
% % calculating the AMPLITUDE (ENVELOPE)
inst_amp_4 = sqrt([real(H_4).*real(H_4)]+[imag(H_4).*imag(H_4)]);
% % calculating the FREQUENCY
inst_freq_4 = gradient(atan(imag(H_4)./real(H_4)))/(1/(2*pi));
G_4=inst_freq_4(1:r_4);
inst_freq_4=inst_freq_4-(max(G_4)+min(G_4))/2;
G_4=inst_freq_4(1:r_4);
inst_freq_4=10* inst_freq_4/max(G_4);
% %calculating the MOMENT OF VELOCITY
mov_4= inst_freq_4.*(inst_amp_4.*inst_amp_4);



% % Patient 5 STAFF-III Database. The STAFF III database was acquired 
% during 1995–96 at Charleston Area Medical Center (WV, USA) where single 
% prolonged balloon inflation had been introduced to achieve optimal results 
% of percutaneous transluminal coronary angiography (PTCA) procedures, 
% replacing the typical series of brief inflations. The database consists 
% of standard 12-lead ECG recordings from 104 patients.
s_5=awgn(s_5,20) ;  % adding some white Gussian noise
r_5=floor(length(s_5)/6.25);
m_5=s_5(1:r_5);
s_5=s_5/(max(m_5)-min(m_5));

t=(0:length(s_5)-1);

m_5=s_5(1:r_5);
s_5=s_5-min(m_5);
s_5=s_5-mean(s_5);  % eliminating DC component

% % performing the HILBERT TRANSFORM
H_5 = hilbert(s_5);
HT_5=imag(H_5);
% % calculating the AMPLITUDE (ENVELOPE)
inst_amp_5 = sqrt([real(H_5).*real(H_5)]+[imag(H_5).*imag(H_5)]);
% % calculating the FREQUENCY
inst_freq_5 = gradient(atan(imag(H_5)./real(H_5)))/(1/(2*pi));
G_5=inst_freq_5(1:r_5);
inst_freq_5=inst_freq_5-(max(G_5)+min(G_5))/2;
G_5=inst_freq_5(1:r_5);
inst_freq_5=10* inst_freq_5/max(G_5);
% %calculating the MOMENT OF VELOCITY
mov_5= inst_freq_5.*(inst_amp_5.*inst_amp_5);




% % Patient 6 PATIENT DATA:
% 
% PATIENT DATA:
% AGE/SEX:
%   41 year old male
% DIAGNOSIS/SURGERY:
%   Congestive heart failure
% PERTINENT HISTORY:
%   Severe mitral stenosis
%   Moderate mitral regurgitation
%   Aortic stenosis
%   s/p AVR and MVR 10 years earlier
s_6=awgn(s_6,20) ;  % adding some white Gussian noise
r_6=floor(length(s_6)/3.125); %time interval [1.6 sec]
m_6=s_6(1:r_6);
s_6=s_6/(max(m_6)-min(m_6));

t=(0:length(s_6)-1);

m_6=s_6(1:r_6);
s_6=s_6-min(m_6);
s_6=s_6-mean(s_6);  % eliminating DC component

% % performing the HILBERT TRANSFORM
H_6 = hilbert(s_6);
HT_6=imag(H_6);
% % calculating the AMPLITUDE (ENVELOPE)
inst_amp_6 = sqrt([real(H_6).*real(H_6)]+[imag(H_6).*imag(H_6)]);
% % calculating the FREQUENCY
inst_freq_6 = gradient(atan(imag(H_6)./real(H_6)))/(1/(2*pi));
G_6=inst_freq_6(1:r_6);
inst_freq_6=inst_freq_6-(max(G_6)+min(G_6))/2;
G_6=inst_freq_6(1:r_6);
inst_freq_6=10* inst_freq_6/max(G_6);
% %calculating the MOMENT OF VELOCITY
mov_6= inst_freq_6.*(inst_amp_6.*inst_amp_6);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   Plot Noisy  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Noise %%%%%%%%%%%  Figures 14  %%%%%%%%%%%%  Signals & IF  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(14);

fig = figure(14);
left_color = [0.85 0.325 0.0980];
right_color = [0 0.4470 0.741];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);




subplot(2,2,1), 

yyaxis left
plot(t(1,1:r_1),s_1(1,1:r_1),'LineWidth',1.5),title('(a)', 'FontSize', 16)
xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_1 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t(1,1:r_1),inst_freq_1(1,1:r_1),'LineWidth',1),title('(a)', 'FontSize', 16)
xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Instantaneous Frequency (Hz)', 'FontSize', 16) % y-axis label
axis([0 r_1 -10 10])
ylim([-30 10])
yticks([-10 -5 0 5 10])




subplot(2,2,2), 

yyaxis left
plot(t(1,1:r_2),s_2(1,1:r_2),'LineWidth',1.5),title('(b)', 'FontSize', 16)
xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_2 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t(1,1:r_2),inst_freq_2(1,1:r_2),'LineWidth',1),title('(b)', 'FontSize', 16)
xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Instantaneous Frequency (Hz)', 'FontSize', 16) % y-axis label
axis([0 r_2 -10 10])
ylim([-30 10])
yticks([-10 -5 0 5 10])





subplot(2,2,3), 
yyaxis left
plot(t(1,1:r_3),s_3(1,1:r_3),'LineWidth',1.5),title('(c)', 'FontSize', 16)
xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_3 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t(1,1:r_3),inst_freq_3(1,1:r_3),'LineWidth',1),title('(c)', 'FontSize', 16)
xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Instantaneous Frequency (Hz)', 'FontSize', 16) % y-axis label
axis([0 r_3 -10 10])
ylim([-30 10])
yticks([-10 -5 0 5 10])






subplot(2,2,4),
yyaxis left
plot(t,s_6,'LineWidth',1.5),title('(d)', 'FontSize', 16)
xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_6 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t,inst_freq_6,'LineWidth',1),title('(d)', 'FontSize', 16)
xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Instantaneous Frequency (Hz)', 'FontSize', 16) % y-axis label
axis([0 r_6 -10 10])
ylim([-30 10])
yticks([-10 -5 0 5 10])




h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'Position',[50 50 1200 800]);
print(gcf, '-dpdf', 'SigVsIFN.pdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Noise %%%%%%%%%%  Figures 15  %%%%%%%%%%%%  Signals & MOV  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(15);

fig = figure(15);
left_color = [0.85 0.325 0.0980];
right_color = [0 0.4470 0.741];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);




subplot(2,2,1), 

yyaxis left
plot(t(1,1:r_1),s_1(1,1:r_1),'LineWidth',1.5),title('(a)', 'FontSize', 16)
xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_1 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t(1,1:r_1),mov_1(1,1:r_1),'LineWidth',1),title('(a)', 'FontSize', 16)
xticks([0 r_1/4 r_1/2 r_1*3/4 r_1])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Moment of Velocity', 'FontSize', 16) % y-axis label
axis([0 r_1 -10 10])
ylim([-25 8])
yticks([-5 -2.5 0 2.5 5])




subplot(2,2,2), 

yyaxis left
plot(t(1,1:r_2),s_2(1,1:r_2),'LineWidth',1.5),title('(b)', 'FontSize', 16)
xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_2 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t(1,1:r_2),mov_2(1,1:r_2),'LineWidth',1),title('(b)', 'FontSize', 16)
xticks([0 r_2/4 r_2/2 r_2*3/4 r_2])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Moment of Velocity', 'FontSize', 16) % y-axis label
axis([0 r_2 -10 10])
ylim([-25 8])
yticks([-5 -2.5 0 2.5 5])





subplot(2,2,3), 
yyaxis left
plot(t(1,1:r_3),s_3(1,1:r_3),'LineWidth',1.5),title('(c)', 'FontSize', 16)
xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_3 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t(1,1:r_3),mov_3(1,1:r_3),'LineWidth',1),title('(c)', 'FontSize', 16)
xticks([0 r_3/4 r_3/2 r_3*3/4 r_3])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Moment of Velocity', 'FontSize', 16) % y-axis label
axis([0 r_3 -10 10])
ylim([-25 8])
yticks([-5 -2.5 0 2.5 5])






subplot(2,2,4),
yyaxis left
plot(t,s_6,'LineWidth',1.5),title('(d)', 'FontSize', 16)
xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Voltage Amplitude', 'FontSize', 16) % y-axis label
axis([0 r_6 -1 1])
ylim([-1 3])
yticks([-1 -.5 0 .5 1])

yyaxis right
plot(t,mov_6,'LineWidth',1),title('(d)', 'FontSize', 16)
xticks([0 r_6/4 r_6/2 r_6*3/4 r_6])
xticklabels({'0','0.5','1','1.5','2'})
xlabel('Time (seconds)', 'FontSize', 16) % x-axis label
ylabel('Moment of Velocity', 'FontSize', 16) % y-axis label
axis([0 r_6 -10 10])
ylim([-25 8])
yticks([-5 -2.5 0 2.5 5])




h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'Position',[50 50 1200 800]);
print(gcf, '-dpdf', 'SigVsMoVN.pdf')
