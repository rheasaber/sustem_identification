
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

error2 = sum((yp - y_predict).^2);
error2

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

ym = ym(3:end);
error = sum((yp - ym).^2);
error

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
figure(3)
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
plot(time,y_predict);
legend('y ','y_{predict}');
xlabel('Time(s)');
ylabel('Different outputs');
hold off



error = sum((y_predict - y).^2);
error



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
p=8;

%bode diagram
Fu = fft(U);
Fy = fft(Y);
g = Fy(1:p:end)./Fu(1:p:end);
n=N/p;
omega=2*pi/Ts;
freq_vect = 0 : omega*(1/n) : (N-1)*omega*(1/n);
frequency_model= frd(g(1:(n-1)/2),freq_vect(1:(n-1)/2),Ts);

figure(5);
hold on
bode(frequency_model, freq_vect);
hold off

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
figure(6)

plot(loss)
xlabel('Order Number')
ylabel('Loss Value')

% ARMAX
min =1;
max=6;
for i=min:max
    %figure()
    model_armax= armax(data_objct_mean,[i i i 1]);
   
    %h = iopzplot(model_armax);

    showConfidence(h,2);

end % delta = 4 seems to be the max value for wich there is no pole /zero cancellation for fixedpoles and zeros and (thus delta_0=4=order)

% estimation of time delay
delta= 4;%=order
model_armax= armax(data_objct_mean,[delta delta delta 1]);
figure();
errorbar(model_armax.b,model_armax.db)
xlabel('Coefficient of B')
ylabel('Magnitude')

%nk=1 because bo=0
%1<= nb<=delta-nk+1   thus nb<= delta=4 but how much?-->loss fcn
%na=n=4
%nb=m-d+1  =m

%Compute the loss function for variable nb --> seems to be equal to 3
min = 4;
max = 7;
lossb = [];
for i=min:max
    model = arx(data_objct_mean,[4 i 1]);
    lossb = [lossb model.EstimationInfo.LossFcn];  
end

figure(11)
hold on
plot(lossb)
xlabel('Order Number for n_{b}')
ylabel('Loss Value')
hold off

%

na = 1:1:10;
nb = 1:1:10;
nk = 1:5;
nn = struc(na,nb,nn);
V = arxstruc(data_objct_mean,data_objct_mean,nn)
selstruc(V)

%% 2.2.2 Parametric identification

load("ASSdata.mat");
Ts=0.01;
data = iddata(y,u,Ts);
data = detrend(data);

N = length(u);
data_id = data(1:N/2);
data_val = data(N/2+1:end);

na=4;
nc=na;
nd=na;
nf=na;
nx=na;
nb=3;
nk=1;

model_arx = arx(data_id, [na, nb, nk]);
model_iv4 = iv4(data_id, [na, nb, nk]);
model_armax = armax(data_id, [na, nb, nc, nk]);
model_oe = oe(data_id, [nb, nf, nk]);
model_bj = bj(data_id, [nb, nc, nd, nf, nk]);
model_n4sid = n4sid(data_id, nx);

%time domain
figure
compare(data_val, model_arx, model_armax, model_oe, model_bj, model_iv4, model_n4sid);

%frequency domain
U = fft(u);
Y = fft(y);
p=8;
n=N/p; %=125

g = Y(1:p:end)./U(1:p:end);
omega=2*pi/Ts;
freq_vect = 0 : omega*(1/n) : (N-1)*omega*(1/n);
G= frd(g(1:(n-1)/2),freq_vect(1:(n-1)/2),Ts);

figure
compare(G, model_arx, model_armax, model_oe, model_bj, model_iv4, model_n4sid);
