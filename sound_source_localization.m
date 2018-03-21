%% clearing and closing commands
clc
clear all 
close all
%% Assuming the source location,sampling frequency,number of microphones,distance between each microphone
c=343;% speed of sound in air
fs=10000;%sampling frequncy
source=[8 6];%source location
N_mics=5;% Number of microphones
d=2;% distance between each microphone
%% setting up Microphone positions
for i=0:N_mics-1;
    mics(i+1,:)=[d*i 0];
end
%% Estimating the time delay
for i=1:length(mics)
    X=[source(1,1) source(1,2);mics(i,1) mics(i,2)];
    mic_source_distance(i)=pdist(X,'euclidean');
    mics_delay_time(i)=(mic_source_distance(i)-mic_source_distance(1))/c;%Time delay at each microphone in seconds 
    mics_delay_samples(i)=mics_delay_time(i)*fs;%Time delay at each microphone in samples
end
%% Delaying the signal using sinc and filter functions
N=10000;
X0=randn(N,1);% generating random point sound source signal
filter_length=400;% filter length
figure;
for i=1:N_mics-1
    w=sinc([-filter_length/2:filter_length/2]-mics_delay_samples(1,i+1));%sinc filter used for delaying the signal 
    Xdelayed(:,i)=filter(w,1,X0); % delayed signal
    plot(w);
    hold on
end
%% plotting impulse response of filter
grid on
hold off
title('Impulse Response of the filter');
xlabel('Samples(n)');
ylabel('Impulse Response');
legend('mic1','mic2','mic3','mic4')
%% Estimating the Time delay using the LMS algorithm
w1=zeros(filter_length,N_mics-1);
mu=0.001; %step size
for i=filter_length:N
    X1=X0(i:-1:(i-filter_length+1));
    e=Xdelayed(i,:)-X1.'*w1; %error
    w1=w1+mu*X1*e; %updating the filter coefficients
end
%% plotting the Impulse Response of filter coefficients obtained in LMS algorithm
figure;
plot(w1) 
title('Impulse Response of filter coefficients using LMS algorithm');
xlabel('Samples(n)');
ylabel('Impulse Response');
legend('mic1','mic2','mic3','mic4')
grid on
%% Estimating the Time Delay in seconds and in samples
w_fft=fft(w1);%fft of filter coefficient for finding the phase of the filter
phase=unwrap(angle(w_fft));%phase of the filter
p=0:(pi/filter_length):(pi-(pi/filter_length));%normalizing the phase
t=[p;ones(1,filter_length)].';
for i=1:N_mics-1;
    phase_slope=t\phase(:,i);%slope of phase of filter
    k(i,1)=phase_slope(1,:);
    T(i,1)=-((filter_length/2)+k(i,1)./2);%Entire delay in samples
end
mics_delay_time1=T/fs;%Entire delay in seconds
t1=p/pi;
%% plotting the phase response of the filter
figure;
plot(t1,phase)
title('Phase response of the filter');
xlabel('w/pi');
ylabel('Phase Response');
legend('mic1','mic2','mic3','mic4')
grid on
%% Source Localization using Steepest Descent algorithm
syms x y %creating symbolic variables
G=0;
dx=0;
dy=0;
for i=1:N_mics-1
    G1(i)=(((sqrt((x-(i*d)).^2 + (y).^2) - sqrt((x).^2 + (y).^2))) - (mics_delay_time1(i,1)*c)).^2;     
    G=G+G1(i);
end
G2=matlabFunction(G);
xdiff=diff(G2,x);
ydiff=diff(G2,y);
Fx=matlabFunction(xdiff);
Fy=matlabFunction(ydiff);
initiate_source=[1;1];%initiating the source location
mu1=0.01;%step size
estimated_source(:,1)=initiate_source-((mu1).*[Fx(initiate_source(1,1),initiate_source(2,1));Fy(initiate_source(1,1),initiate_source(2,1))]);
for i=2:5000
estimated_source(:,i)=estimated_source(:,i-1)-((mu1).*[Fx(estimated_source(1,i-1),estimated_source(2,i-1));Fy(estimated_source(1,i-1),estimated_source(2,i-1))]);
end
disp(estimated_source(:,i));%displaying the estimated source
%% plotting the microphone location,source location and estimated source location  
figure;
plot(source(1,1),source(1,2),'*',mics(:,1),mics(:,2),'+')
axis([-5 20 -5 15])
hold on;
plot(estimated_source(1,i),estimated_source(2,i),'o')
title('2-D coordinate system')
legend('source location','mic location','Estimated source')
xlabel('X----->');
ylabel('Y----->');
grid  on
%% plotting the Calculated Time delay and Estimated Time delay
mics_delay_time1=[0 mics_delay_time1.'];
y=0;
figure;
plot(mics_delay_time,y,'*')
hold on
plot(mics_delay_time1,y,'o')
title('Calculated time delay and Estimated time delay')
xlabel('X----->');
ylabel('Y----->');
legend('Calculated time delay','Estimated time delay')
grid on


