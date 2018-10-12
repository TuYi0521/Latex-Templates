clc;
clear all;
close all;
% format long e
%%
n=3;
h=1./(n+1);
h4=h.^4;
x1 = linspace ( 0, 1.0, n+2);
y1 = linspace ( 0, 1.0, n+2);

x = x1(2:(n+1));
y = y1(2:(n+1));

v = zeros(n.^2,1);
u = v;
f = v;
Fxx=f;
Fyy=f;

A = zeros(n.^2,n.^2);
C = A;

syms a b;
F(a,b)= sin(pi.*a) * (pi.^2 * b * (exp(b) - exp(1)) - (2*exp(b) + b * exp(b))); 
Faa(a,b)=diff(F(a,b),a,2);
Fbb(a,b)=diff(F(a,b),b,2);

for j=1:n
    for i=1:n
        u(i+(j-1)*n) = sin(pi * x(i)) * y(j) * (exp(y(j)) - exp(1));
        f(i+(j-1)*n) = sin(pi.*x(i)) * (pi.^2 * y(j) * (exp(y(j)) - exp(1)) - (2*exp(y(j)) + y(j) * exp(y(j)))); 
%         f(i+(j-1)*n) = -((-pi.*pi.*sin(pi.*x(i)).*y(j).*(exp(y(j))-exp(1)))+(sin(pi.*x(i)).*exp(y(j)).*(2+y(j))));
    end   
end

for j=1:n
    for i=1:n
        Fxx(i+(j-1)*n)=Faa(x(i),y(j));
        Fyy(i+(j-1)*n)=Fbb(x(i),y(j));
    end   
end

B2 = diag(8*ones(n,1))+ diag(ones(n-1,1),1)+ diag(ones(n-1,1),-1);
I = eye(n,n);

for i=1:n
    for j=1:n
        for k=1:n
            C(n*(i-1)+j,n*(i-1)+k)=B2(j,k);
        end
    end
end

for i=1:(n-1)
    for j=1:n
        for k=1:n
            C(n*(i-1)+j,n*i+k)=I(j,k);
        end
    end
end

for i=1:(n-1)
    for j=1:n
        for k=1:n
            C(n*i+j,n*(i-1)+k)=I(j,k);
        end
    end
end


B1 = diag((-20)*ones(n,1))+ diag(4*ones(n-1,1),1)+ diag(4*ones(n-1,1),-1);
B3 = diag(4*ones(n,1))+ diag(ones(n-1,1),1)+ diag(ones(n-1,1),-1);
for i=1:n
    for j=1:n
        for k=1:n
            A(n*(i-1)+j,n*(i-1)+k)=B1(j,k);
        end
    end
end

for i=1:(n-1)
    for j=1:n
        for k=1:n
            A(n*(i-1)+j,n*i+k)=B3(j,k);
        end
    end
end

for i=1:(n-1)
    for j=1:n
        for k=1:n
            A(n*i+j,n*(i-1)+k)=B3(j,k);
        end
    end
end




A1=A/(6*(h^2));
C1=C/12;

A1_inv=inv(A1);

d = (h^2/12)*(Fxx+Fyy)+f;
v=-A1_inv*d;
v;

Eh=(max(abs(u-v)));

error1 = -A1*u;
error2 =d;
error=max(abs(error1-error2));

V=zeros(n+2,n+2);
U=zeros(n+2,n+2);
for i=1:n
    for j=1:n
       V(i+1,j+1)=v((i-1)*n+j);
       U(i+1,j+1)=u((i-1)*n+j);
    end
end

s=sprintf('V(x,y), n=%d',n);

figure;
mesh ( x1, y1, V);
xlabel ( '<--x-->' );
ylabel ( '<--y-->');
title (s);

s=sprintf('U(x,y), n=%d',n);
figure;
mesh ( x1, y1, U);
xlabel ( '<--x-->' );
ylabel ( '<--y-->');
title (s);

