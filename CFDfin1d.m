%%%%%%%%%%%% p  1      %%%%%%%%%%%%%%%%%%% 
clc
clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------parameters
l = 1;
node = 100;
h = 250;
k = 12;

r = 0.006;
m=900;
To=100-25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=zeros(1,node);
T_exact=zeros(1,node);
D=zeros(1,node-1);
A=zeros(1,node-2);
C=zeros(1,node-2);
TDMA=zeros(node-1,node-1);
bb=zeros(1,node-1);
dd=zeros(1,node-1);
L=zeros(node-1,node-1);
U=zeros(node-1,node-1);
TT=zeros(m,node-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = l/(node-1);

x(1) = 0;
for i = 1:node-1
    
    x(i+1) = x(i)+dx;
end
c1 = 2*h/(k*r);
n=node-1;

c2 = (c1)^0.5 ; 
nx=length(x);
for ix = 1 : nx
    T_exact(ix) = (75)*(cosh(c2*(x(ix)-l))+(h/(c2*k))*sinh(c2*(l-x(ix))))/(cosh(c2*l)+(h/(c2*k))*sinh(c2*l)) ; 
end
for i=1:n
    D(i)=-(2+c1*dx^2);
     D(n)=-(2+c1*dx^2+2*dx*h/k);
end
for i=1:n-1
    A(i)=1;
    C(i)=1;
     C(n-1)=2;
end
for i=1:n
    TDMA(i,i)=D(i);
    if i<n
    TDMA(i,i+1)=A(i);
    end
    if i<n
    TDMA(i+1,i)=C(i);
    end
end

%--------------------------------------------
d=D;
a=[0,A];
b=[C,0];

rh=zeros(n,1);
rh(1)=-To;

Tprimary=zeros(n,1);
for i=1:n
  Tprimary(i)=0;
end
rhs=rh+Tprimary;
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   bb(1) =0;
dd(1) = d(1);
for i = 2 : n
	bb(i) = b(i) / dd(i-1);
	dd(i) = d(i) - bb(i) * a(i-1);
end
%--------------------------------------------
for i=1:n
    L(i,i)=1;
   if i<n
    L(i+1,i)=bb(i+1);
   end
end
for i=1:n
    U(i,i)=dd(i);
  if i<n
    U(i,i+1)=a(i);
  end
end

%--------------------------------------------
for i=1:m
y=L\rhs;
T=U\y;
rhs=T+rh;
 TT(i,:)=T;
end
%--------------------------------------------

T1=zeros(m,node);
T1(:,1)=To;
T1(:,2:end)=TT;



%--------------------------------------------
figure(1)
plot(x,T1(1,:),x,T_exact)
xlabel('x(m)') ; 
ylabel('Tetha [{\circC}]')
legend('T','T_e_x_a_c_t')
grid on 

figure(2)
Tn=[0.8254239,0.00077051,0.0075578,0.0071318];
Te=0.0071419;
Ter=abs(Tn-Te);
h_x=[0.1,0.01,0.001,0.0001];
loglog(h_x,Ter)
xlabel('h_x') ; 
ylabel('TR')

figure(3)
Ters=[0.8254239-0.0077051,0.0077051-0.0075578,0.0075578-0.0071318];
h_xs=[0.1,0.01,0.001];
semilogy(h_xs,Ters)
xlabel('h') ; 
ylabel('successive error')
