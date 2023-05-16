
%% 2.1.1 FIR model identification
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
legend('y measured','y predicted');
xlabel('Time(s)');
ylabel('Different outputs');
hold off

%plot FIR(theta)
errorbar(Theta_LS, 2*sqrt(diag(covar)));

%% 2.1.2 ARX model identification

[y,u,Ts]=GetExperimentData('logs.bin');
N = length(u);
M = length(y);
%time_vector=[0:Ts:N*Ts];
time_vector=linspace(0,N*(1/100),N-2);
m=200;

%phi
z = zeros(1,m); %All previous inputs are assumed to be equal to zero up ..
Phi = [-y(2:M-1),-y(1:M-2),u(2:M-1),u(1:M-2)]; 
yp=y(3:end);

Theta= inv(Phi'*Phi)*(Phi'*yp);

%Compute the predicited output
y_predict = Phi*Theta;

%Loss fonction
J= sum((yp - y_predict).^2);


%plot y and y_predict
figure(2)
hold on
plot(time_vector,yp);
plot(time_vector,y_predict);
legend('y','y predict');
xlabel('Time(s)');
ylabel('Different outputs');
hold off

%ym
time=0:Ts:Ts*(N-1);
Bo=Theta(3:end);
Ao=Theta(1:2);
sys=tf([0 Bo'],[1 Ao'],Ts);
ym=lsim(sys,u,time);

%plot y and y_predict
figure(3)
hold on
plot(time_vector,yp);
plot(time,ym);
legend('y ','ym');
xlabel('Time(s)');
ylabel('Different outputs');
hold off

error == sum((yp - ym).^2);

%% 2.1.3 state- space model identification

r=20;
N = length(u);
%Y and U
[y,u,Ts]=GetExperimentData('logs.bin');
Y = buffer(y, r, r-1, 'nodelay');
U = buffer(u, r, r-1, 'nodelay');
Y = Y(:, 1:end-1);
U = U(:, 1:end-1);

n=length(U);

U_p=eye(n) - U'*inv(U*U')*U;
Q=Y*U_p;

% compute the singular values using svd
SV=svd(Q);
%number of states
plot(SV, 'o')

n=2; % according to the graphe try to see with 2,3 or 4
Or=Q(:,1:n);

C=Or(1,:);
A=Or(1:r-1,:)\Or(2:r,:);

%estimate D annd B  using LS
D=0;

q = tf('z',Ts);
F = C*inv(q*eye(size(A))-A);

Uf = [];  % Initialize an empty matrix to store the simulation outputs
time=0:Ts:Ts*(N-1);
for i = 1:2
    F_i = F(i);  % Retrieve the i-th element of vector F
    output = lsim(F_i, u, time);  % Perform linear simulation
    
    Uf = [Uf, output];  % Append the output to Uf matrix
end


Phi=Uf';
B = (Phi*Phi')\(Phi*y);

%compute ym
SS = ss(A,B,C,0,Ts);
sys = tf(SS);
y_predict=lsim(sys,u,time);

%plot y and y_predict
figure(4)
hold on
plot(time,y);
plot(time,ym);
legend('y ','ym');
xlabel('Time(s)');
ylabel('Different outputs');
hold off



%% Exercice 2.2.1 Order estimation
load("ASSdata.mat")
Ts=0.01;
data_objct=iddata(y,u,Ts);
y_o= detrend(y);
u_o = detrend(u);

data_objct_mean=iddata(y_o,u_o,Ts);
Y = data_objct_mean.y;
U = data_objct_mean.u;
N=length(U);

Fu = fft(U);
Fy = fft(Y);

%ARX
nk=1;
loss =[];
lossFcn = zeros(10);
for i=1:10
    model_arx = arx(data_objct_mean, [i i 1]);  % Identify ARX model of specific order
    loss = [loss model_arx.EstimationInfo.LossFcn];  
end
%order seems to be 2

%Plot loss function ARX
figure(5)
close all
plot(loss)
xlabel('Order Number')
ylabel('Loss Value')

%ARMAX
min =3;
max=8;
for i=min:max
    model_armax= armax(data_objct_mean,[i i i 1]);
    subplot(max-min+1,1,i-min+1)
    h = iopzplot(model_armax);
    showConfidence(h,2);
end

%estimation of time delay
delta= 4;%=order
model_armax= armax(data_objct_mean,[delta delta delta 1]);
errorbar(model_armax.b,model_armax.db)
xlabel('Coefficient of B')
ylabel('Magnitude')

%1<= nb<=delta-nk+1   thus nb<= delta=4 but how much?


%Compute the loss function for variable nb
min = 1;
max = 7;
lossb = [];
for i=min:max
    model = arx(data_objct_mean,[4 i 1]);
    lossb = [lossb model.EstimationInfo.LossFcn];  
end

%Compute the loss function for variable na
min = 1;
max = 7;
lossa = [];
for i=min:max
    model = arx(data_objct_mean,[i 4 1]);
    lossa = [lossa model.EstimationInfo.LossFcn];  
end
close all

figure(11)
plot(lossb)
xlabel('Order Number for n_{b}')
ylabel('Loss Value')

%Looking at the plots it seems that nb is around 3 
close all
figure(12)
plot(lossa)
xlabel('Order Number for n_{a}')
ylabel('Loss Value')


%Looking at the plots it seems that na is around 2

NA = 1:1:10;
NB = 1:1:10;
NK = 1:5;
NN = struc(NA,NB,NK);
V = arxstruc(data_objct_mean,data_objct_mean,NN);
selstruc(V)
