clc
clear all
% the wavelength of the ligh.
lambda = 0.55e-6;
% total number of frequency components.
Nf=80; 
% number of layers.
Nh=40;
% number of frequency-points to be used for minimization.
Np=Nh+30;
% spatial resolution or sampling interval of the speckle image on the analysis plane. 
Dx=0.009;
% the distance below the ground surface to which the analysis plane is shifted.
hsg=1000;
hh=zeros(1,Nh);
Cn=zeros(1,Nh);
ff=zeros(1,Nf);
% Simulated Cn2 h=3 km to be recovered following minimization.
Cn(7)=0.6;
% Simulated Cn2 h=4 km to be recovered following minimization.
Cn(9)=0.5;  
%The simulated experimental power spectrum.
W=zeros(Nh,Nf);
W_exp=zeros(1,Nf);
%%%%%%%%%%%%%%%%%%%%%%%%%
%the thickness of the turbulent layer
for j=1:Nh
         hh(j)=(j-1)*500;        
end

% The spatial frequency
Df=1/(Nf*Dx);
for kk=1:Nf
        ff(kk)=(kk-1)*Df;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The simulated experimental power spectrum (W_exp)

for i=1 : Nh
    h=hh(i)+hsg;  
    for j=1:Nf
        %f=ff(j);
        f=abs(j-Nf/2-1)*Df;
        W(i,j)=(sin(3.14*lambda*h*f^2))^2;%*f^(-11.0/3.0);
        W(i,Nf/2+1)=0.0;       
    end
end

W=W./max(max(W)) + 0.0*rand;

for mf=1:Nf 
 
    for kh=1:length(hh)
        W_exp(mf)=W_exp(mf)+Cn(kh)*W(kh,mf);
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Objectif Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   fun=@(x)objf(x,hh,ff,hsg,W_exp,Np); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x0: This variable is initialized as an array of zeros and is likely used as the initial guess or starting point for optimization.
x0 = zeros(1,Nh);
%lb: This variable is initialized as an array of zeros and likely represents the lower bounds for the optimization problem.
lb = zeros(1,Nh); 
%ub: This variable is initialized as an array of twos (2s) and likely represents the upper bounds for the optimization problem.
ub=2.*ones(1,Nh);
%A: This variable is initialized as an array of ones and is used in the optimization problem, as part of linear constraints.
A=ones(1,Nh); 
%b: This variable is calculated as the sum of the upper bounds plus 1.
b=sum(ub)+1; 
%x: This variable stores the solution obtained from the optimization process.
%fval: This variable stores the value of the objective function (the result of the optimization).
[x fval]= fmincon(fun,x0,A,b,[],[],lb,ub);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%b: This variable is generated as random values in an array, likely for use as noise.
b=rand(1,Nf);
%snr0: This variable calculates the Signal-to-Noise Ratio (SNR) by taking the maximum value of W_exp divided by the standard deviation of b.
snr0=max(W_exp)/std(b);
%noise: This variable represents the noise component, calculated as 1.5 times b divided by snr0.
noise=1.5*b./snr0;
W_exp=W_exp+noise;
SNR=std(W_exp)/std(noise);
SNR2=mean(W_exp)/std(noise);
SNR3=5*10*log10(SNR2);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Objectif Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   fun=@(x)objf(x,hh,ff,hsg,W_exp,Np); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = zeros(1,Nh); % %%% start point
lb = zeros(1,Nh); %%%%% min of Cn2
ub=3.*ones(1,Nh);
A=ones(1,Nh); 
b=sum(ub)+1; 
[x1 fval]= fmincon(fun,x0,A,b,[],[],lb,ub);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure (1)
x=x.*10^(-16);
x1=x1.*10^(-16);
Cn=Cn.*10^(-16);
hh=hh/1000;
plot(hh, Cn, 'r-x', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on
plot(hh,x, 'b-+', 'LineWidth', 1.5)
plot(hh,x1,'g-s', 'LineWidth', 1.5,'MarkerSize', 8)
xlim([0 20])
legend('C_n^2(h): Simulated', ['C_n^2(h): Retrieved, SNR= inf'], ['C_n^2(h): Retrieved, SNR=', num2str(floor(SNR3))])
xlabel('h(km)', 'FontSize', 14);
ylabel('C_n^2(h) (m^{-2/3})', 'FontSize', 14);
title(strcat('SSS- C_n^2(h); Nh=40, \delta h=0.5km') ,'FontSize', 14);




