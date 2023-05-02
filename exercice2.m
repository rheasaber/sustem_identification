
%% 2.1.1
[y,u,Ts]=GetExperimentData('logs.bin');
N = length(u);
%time_vector=[0:Ts:N*Ts];
time_vector=linspace(0,N*Ts,N);
m=200;

%phi
z = zeros(1,m); %All previous inputs are assumed to be equal to zero up ..
Phi = [0;u(1:end-1)]; 
Phi = toeplitz(Phi,z);

Theta_LS= inv(Phi'*Phi)*(Phi'*y);

%Compute the predicited output
y_predict = Phi*Theta_LS;

%Loss fonction
J= sum((y - y_predict).^2);

%covariance
var=J/(N-m);
covar=var*inv(Phi'*Phi);

%plot y and y_predict
figure(1)
hold on
plot(time_vector,y);
plot(time_vector,y_predict);
legend('y measured','y predictes');
xlabel('Time(s)');
ylabel('Different outputs');
hold off

%plot FIR(theta)
errorbar(Theta, 2*sqrt(diag(covar)));

%% 2.1.2 
[y,u,Ts]=GetExperimentData('logs.bin');
N = length(u);
%time_vector=[0:Ts:N*Ts];
time_vector=linspace(0,N*(1/100),N);
m=200;

%phi
z = zeros(1,m); %All previous inputs are assumed to be equal to zero up ..
Phi = [-y(2:end-1),-y(1:end-1),u(2:end-1),u(1:end-1)]; 
y=y(3:end);

Theta= (Phi'\Phi)*(Phi'*y);

%Compute the predicited output
y_predict = Phi*Theta_LS;

%Loss fonction
J= sum((y - y_predict).^2);


%plot y and y_predict
figure(2)
hold on
plot(time_vector,y);
plot(time_vector,y_predict);
legend('y measured','y predictes');
xlabel('Time(s)');
ylabel('Different outputs');
hold off

%ym
Bo=Theta(3:end);
Ao=Theta(1:2);
sys=tf([0 Bo'],[1 Ao'],Ts);
ym=lsim(sys,u,timevector);

%plot y and y_predict
figure(2)
hold on
plot(time_vector,y);
plot(time_vector,ym);
legend('y measured','ym');
xlabel('Time(s)');
ylabel('Different outputs');
hold off

error == sum((y - ym).^2);