
%% exercice 1.1 

Ts=0.2; %sec
tfc = tf([0, 3], [0.5, 0.4, 1]);% transfer fonction continuous
tfd=c2d(tfc, Ts, 'zoh');

%step fcn
simin.time = (0:Ts:10)';
%simin.signals.values =[zeros(1/Ts,1); 0.7*ones(fix((10-1)/(Ts+1)),1)];
simin.signals.values = 0.7*(simin.time >= 1);
simulation_time=10;

%simulation
simout=sim('compexx',simulation_time);

%plot step response
opt = stepDataOptions('StepAmplitude',0.7);
[true_step, t_step] = step(tfd, 10, opt);% without the noise

figure(1);
hold on

plot(simout.simout.time,simout.simout.signals.values);
plot(simin.time,simin.signals.values);%step
plot([0;t_step+1],[0;true_step]);
legend('Real Output with noise', 'Input Step', 'Noiseless step response');
title('Plot the step response of the model');
xlabel('Time(s)');
ylabel('Data');
xlim([0,10]);

hold off

%impulse fcn
simin.time = (0:Ts:10)';
%simin.signals.values =[zeros(1/Ts,1);0.7*ones(1/Ts,1); zeros(fix((10-2)/(Ts+1)),1)];
simin.signals.values =0.7*(simin.time == 1)
simulation_time=10;

%simulation
simout=sim('compexx',simulation_time);

%plot impulse response
%opt = stepDataOptions('ImpulseAmplitude',0.7);
[true_impulse, t_impulse] = impulse(tfd);% without the noise

figure(2);
hold on

plot(simout.simout.time,simout.simout.signals.values);
plot(simin.time,simin.signals.values);%step
plot([0;t_impulse+1],0.7*[0;true_impulse]);
legend('Real Output with noise', 'Input Impulse', 'Noiseless impulse response');
title('Plot the impulse response of the model');
xlabel('Time(s)');
ylabel('Data');
xlim([0,10]);

hold off


%% exercice 1.2: Auto Correlation of a PRBS signal
u=prbs(7,5);

[R,h] = intcor(u,u);

% PLOT 
figure(3);
hold on
stem(h,R);
title('Autocorrelation function of a PRBS signal');
xlabel('h');
ylabel('R(u,u)')
hold off

%% exercice 1.3:

Ts= 0.2; %sec
N=50/Ts;
time_vector=[0:Ts:(N-1)*Ts];

random_signal= rand(N, 1);

toeplitz_row = zeros(1,length(random_signal));
toeplitz_row(1) = random_signal(1);
U = toeplitz(random_signal, toeplitz_row);

simin.signals.values = random_signal;
simin.time= time_vector;
simout=sim('compexx',50);

output=simout.simout.signals.values;

%finite impulse response


U_p = U(:,1:200); %k=200
f_impulse_response = (U_p'*U_p) \ U_p' * output(1:end-1);

%finite impulse response with bias
tfc = tf([0, 3], [0.5, 0.4, 1]);% transfer fonction continuous
tfd=c2d(tfc, Ts, 'zoh');
[true_impulse, impulse_t] = impulse(tfd, time_vector(200)); %to compare with the values given bu lambda


%finite impulse response with bias

%eye= eye(200); %matrice identite


lambda=linspace(1,200,400);
erreur_theta_truei=zeros(length(lambda),1);

for i=1:length(lambda)-1
    impulse_bias= ((U_p'*U_p)+lambda(i)*eye(200)) \ (U_p' * output(1:end-1));
    erreur_theta_truei(i)=sqrt(sum((impulse_bias-true_impulse*Ts).^2));
end

%doable in exercices because we have the true impulse response but in
%theory we don't have it so we need to find lambda by trial and error

[~, I] = min(erreur_theta_truei);
lambda= lambda(I);
impulse_bias= (U_p'*U_p+lambda*eye(200)) \ (U_p' * output(1:end-1));


figure(4)
hold on
plot(time_vector(1:length(f_impulse_response)),f_impulse_response);
plot(time_vector(1:length(impulse_bias)),impulse_bias);
plot(impulse_t, true_impulse*Ts); %true impulse
legend('Finite impulse response by deconvolution (K=200)','Finite impulse response by regularization', 'True impulse response');
xlabel('Time(s)');
ylabel('Module of the different inpulse responses ');
hold off
%


%% Exercice 1.4
Ts=0.2;
K=200;

u_14=prbs(7,2);
N=length(u_14);

time_vector=[0:Ts:(N-1)*Ts];


tfc = tf([0, 3], [0.5, 0.4, 1]);% transfer fonction continuous
tfd=c2d(tfc, Ts, 'zoh');
[true_impulse, impulse_t] = impulse(tfd, time_vector(length(u_14))); %to compare with the values given bu lambda




simin.signals.values = u_14;
simin.time= time_vector;

simout=sim('compexx',(N-1)*Ts);

output=simout.simout.signals.values;


%impulse response of the intercorrelation

[Ruu,~] = intcor(u_14,u_14);
[Ryu,~] = intcor(output,u_14);

Ruu_k = Ruu(1:K);%reduced
Ryu_k = Ryu(1:K);

R_toeplitz = toeplitz(Ruu_k, Ruu_k); 
impulse_intercor=R_toeplitz\Ryu_k';

%impulse response of the xcorrelation
XRuu=xcorr(u_14,u_14);
XRuy=xcorr(output,u_14);

XRuu_k = XRuu(N:N+K-1);%reduced
XRuy_k = XRuy(N:N+K-1);

XR_toeplitz = toeplitz(XRuu_k, XRuu_k); 
impulse_xcorr=XR_toeplitz\XRuy_k;


%plot impulse responses 
figure(5);
hold on

plot(time_vector(1:length(impulse_intercor)),impulse_intercor);
plot(time_vector(1:length(impulse_xcorr)),impulse_xcorr);
plot(impulse_t, true_impulse*Ts); %true impulse
legend('impulse_intercor','impulse_xcorr', 'True impulse response');
xlabel('Time(s)');
ylabel('Module of the different impulse responses ');

hold off



