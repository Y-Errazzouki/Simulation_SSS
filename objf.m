function fun=objf(x,hh,ff, hsg,W1,Np)

lambda = 0.55e-6;
Nh0=length(hh);
Nf0=length(ff);
W_theo=zeros(Nh0,Nf0);
Df=ff(2)-ff(1);
% The simulated theorical power spectrum.
for i=1 : Nh0
    h=hh(i)+hsg;  
    for j=1:Nf0
        %f=ff(j);
        f=abs(j-Nf0/2-1)*Df;
        W_theo(i,j)=(sin(3.14*lambda*h*f^2))^2;%*f^(-11.0/3.0);
        W_theo(i,Nf0/2+1)=0.0;       
    end
end
W_theo=W_theo./max(max(W_theo));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Objective Function for Minimization (fun)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fun=0.0;
for m=1:Np 
    s1=0.0;
    for k=1:length(hh)
        s1=s1+x(k)*W_theo(k,m);
    end
fun=fun+abs(s1-W1(m));
end
