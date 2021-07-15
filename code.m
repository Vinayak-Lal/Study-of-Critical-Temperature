%Computational Assignment
%Vinayak Lal 190969
%Chloroform Peng Robinson

%Declaring constants

eps=1-sqrt(2);
sigma=1+sqrt(2);
gamma=0.077802;
psi=0.45724;
Tc=536.4;
Pc=5.470e6;
Vc=0.239e-3;
R=8.314413;
b=gamma*R*(Tc/Pc);
omega=0.218;
kappa=0.37464+1.54226*omega-0.26992*omega*omega;
SP=zeros(1,6);%To store Saturation pressures
SV=zeros(1,11);%Corrsponding Vg and Vl
%taking guess values from Inspection of isotherms
Guess=[5012990,5109610,5190340,5255060,5318070,5403310];
%%Plotting isotherms and points on the saturation curve
for T=531.4:535.4
    i=T-530.4;%counter variable
    alpha=(1+kappa*(1-sqrt(T/Tc)))^2;
    a=alpha*psi*R*R*Tc*Tc/Pc;
    hold on
    %%plotting P vs V using Peng Robinson EOS for 1 mole of chloroform
    syms V;
    fplot((R*T/(V-b))-(a/((V+(eps*b))*(V+(sigma*b)))),[0.0001 0.2]);
    
    xlim([0.0001 0.00055])
    ylim([4.8e6 5.5e6])
    %%Now we have to approximate Psat at temperature T so as to get the
    %%points on the saturation curve
    %Antoine coefficients
    A=13.7324;
    B= 2548.74 ;
    C=218.552;
    %ln(Psat(in Kpa))=A-(B/(T(in C)+C))
    %taking PSat from ANtoine eqn as initial guess
    
    %Psat=(exp(A-(B/(T-273.15+C))))*10^3;
    Psat=Guess(i);
    
    %MUl=chemical potential of liquid,MUg=chemical potential of gas,At
    %equilibrium both are the same and it is enough to take Departure function
    %of Gibbs free energy
    
     %finding values of Vl and Vg for Psat
     %coefficient of cubic equation of v
     a0=1;
     
     a1=b*(sigma+eps-1)-R*T/Psat;
     a2=b*b*(sigma*eps-(sigma+eps))-R*T*b*(1/Psat)*(sigma+eps)+a/Psat;
     a3=-((eps*sigma*b*b*b)+R*T*(1/Psat)*eps*sigma*b*b+a*b/Psat);
     k=[a0 a1 a2 a3];
     v=roots(k);
     v=sort(v);
     Vl=v(1);
     Vg=v(3);
     T
     
     
     
     DepGl=departureliquid(R,T,Psat,Vl,kappa,Tc,alpha,Pc,a,b,eps,sigma);
     DepGg=departuregas(R,T,Psat,Vg,kappa,Tc,alpha,Pc,a,b,eps,sigma);
     ini=abs(DepGg-DepGl);
     if(abs(DepGg-DepGl)<=1e-3)
         SP(i)=Psat;
         SV(2*i-1)=Vl;
         SV(2*i)=Vg;
         Psat
         continue;
     end
     %checking if we should add to this approximate value or not
     Psat=Psat+2000;
     a0=1;
     
     a1=b*(sigma+eps-1)-R*T/Psat;
     a2=b*b*(sigma*eps-(sigma+eps))-R*T*b*(1/Psat)*(sigma+eps)+a/Psat;
     a3=-((eps*sigma*b*b*b)+R*T*(1/Psat)*eps*sigma*b*b+a*b/Psat);
     
     k=[a0 a1 a2 a3];
     Vtemp=roots(k);
     Vtemp=sort(Vtemp);
     Vltemp=Vtemp(1);
     Vgtemp=Vtemp(3);
     DepGltemp=departureliquid(R,T,Psat,Vltemp,kappa,Tc,alpha,Pc,a,b,eps,sigma);
     DepGgtemp=departuregas(R,T,Psat,Vgtemp,kappa,Tc,alpha,Pc,a,b,eps,sigma);

     
     
     if(abs(DepGltemp-DepGgtemp)>ini)
        
         
         while(abs(DepGltemp-DepGgtemp)>=1e-3)
             
             Psat=Psat-10;
             a0=1;
     
             a1=b*(sigma+eps-1)-R*T/Psat;
             a2=b*b*(sigma*eps-(sigma+eps))-R*T*b*(1/Psat)*(sigma+eps)+a/Psat;
             a3=-((eps*sigma*b*b*b)+R*T*(1/Psat)*eps*sigma*b*b+a*b/Psat);
             k=[a0 a1 a2 a3];
             Vtemp=roots(k);
             Vtemp=sort(Vtemp);
             Vltemp=Vtemp(1);
             Vgtemp=Vtemp(3);
             DepGltemp=departureliquid(R,T,Psat,Vltemp,kappa,Tc,alpha,Pc,a,b,eps,sigma);
             DepGgtemp=departuregas(R,T,Psat,Vgtemp,kappa,Tc,alpha,Pc,a,b,eps,sigma);
             
         end
         %%we have obtained values of Psat and corresponding Vl and Vg
         SP(i)=Psat;
         Psat
         SV(2*i-1)=Vltemp;
         SV(2*i)=Vgtemp;
         continue;
         
     end
     %We have to move in same direction
     while(abs(DepGltemp-DepGgtemp)>=1e-2)
             
             Psat=Psat+1000;
             a0=1;
     
             a1=b*(sigma+eps-1)-R*T/Psat;
             a2=b*b*(sigma*eps-(sigma+eps))-R*T*b*(1/Psat)*(sigma+eps)+a/Psat;
             a3=-((eps*sigma*b*b*b)+R*T*(1/Psat)*eps*sigma*b*b+a*b/Psat);
             k=[a0 a1 a2 a3];
             Vtemp=roots(k);
             Vtemp=sort(Vtemp);
             Vltemp=Vtemp(1);
             Vgtemp=Vtemp(3);
             DepGltemp=departureliquid(R,T,Psat,Vltemp,kappa,Tc,alpha,Pc,a,b,eps,sigma);
             DepGgtemp=departuregas(R,T,Psat,Vgtemp,kappa,Tc,alpha,Pc,a,b,eps,sigma);
             
             
     end
     SP(i)=Psat;
     Psat
     SV(2*i-1)=(Vltemp);
     SV(2*i)=(Vgtemp);
     
     

     
     
end
%Let us also plot for Tc
T=536.4;
syms V;
    fplot((R*T/(V-b))-(a/((V+(eps*b))*(V+(sigma*b)))),[0.0001 0.2]);


X=zeros(1,11);
for i=1:5
    X(2*i)=SP(i);
    X(2*i-1)=SP(i);
SV(11)=0.00025;
X(11)=5446910;


% plotting the saved values using a cubic spline
y=X;
x=SV;

xx = 0.0001:.000001:0.0010;
yy = spline(x,y,xx);
plot(x,y,'o',xx,yy);
hold on

     
end
function depart=departureliquid(R,T,Psat,Vl,kappa,Tc,alpha,Pc,a,b,eps,sigma)


Zl=Psat*Vl/(R*T);
Pr=Psat/Pc;
Tr=T/Tc;
dl=1/Vl;

A1=0.457*Pr*alpha/Tr;
B1=0.0778*Pr/Tc;
depart=-log(1-dl*b)-(a/(b*R*T))*(1/(sigma-eps))*log((1+sigma*dl*b)/(1+eps*dl*b))+Zl-1-log(Zl);


end


function depart2=departuregas(R,T,Psat,Vg,kappa,Tc,alpha,Pc,a,b,eps,sigma)


Zg=Psat*Vg/(R*T);
dg=1/Vg;
Tr=T/Tc;


depart2=-log(1-dg*b)-(a/(b*R*T))*(1/(sigma-eps))*log((1+sigma*dg*b)/(1+eps*dg*b))+Zg-1-log(Zg);


end




