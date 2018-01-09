function howard_example_4
% howard_example_4: MATLAB script that solves the PDE
% stored in deglinsys.m, psysbc.m, and degsysinit.m
m = 0;
nx = 100;
x = linspace(-25,25,100);
t = linspace(0,2,20);
% use helper class to improve performance of ODE evaluations
fc = FuncCalculator(@pdegwave, nx);
pdeFunc = @(x,t,u,DuDx) deglinsys(x,t,u,DuDx, fc);
sol = pde1d(m,pdeFunc,@degsysinit,@psysbc,x,t);
u1 = sol(:,:,1);
u2 = sol(:,:,2);
flag = 1;
while flag==1
answer = input('Finished iteration. View plot (y/n)','s');
if isequal(answer,'y')
figure;
hold on;
fig1=plot(x,u1(1,:),'erasemode','xor');
axis([min(x) max(x) -1 1]);
fig2=plot(x,u2(2,:),'r','erasemode','xor');
for k=2:length(t)
set(fig1,'ydata',u1(k,:));
set(fig2,'ydata',u2(k,:));
pause(.5);
end
else
flag=0;
end
end
end

function xprime = degode(t,x);
% DEGODE: Stores an ode for a standing wave
% solution to the p-system.
xprime=[-2*(x(1)+2)-x(2); -(x(1)+2)-2*x(2)-(x(1)^3+8)];
end

function u1bar=pdegwave(x)
%P DEGWAVE: Function that takes input x and returns
% the vector value of a degenerate wave.
% in degode.m
small = .000001;
if x <= -20
u1bar = -2;
u2bar = 0;
else
tspan = [-20 x];
%Introduce small perturbation from initial point
x0 = [-2+small,-small];
odeFunc = @(x, t) degode(t,x);
x=lsode(odeFunc, x0, tspan);
u1bar = x(end,1);
u2bar = x(end,2);
end
end

function [pl,ql,pr,qr]=psysbc(xl,ul,xr,ur,t)
% PSYSBC: Boundary conditions for the linearized
% p-system.
pl=[ul(1);ul(2)];
ql=[0;0];
pr=[ur(1);ur(2)];
qr=[0;0];
end

function value = degsysinit(x);
% DEGSYSINIT: Contains initial condition for linearized
% p-system.
value = [exp(-(x-5)^2);exp(-(x+5)^2)];
end

function [c,b,s] = deglinsys(x,t,u,DuDx, fc)
% DEGLINSYS: MATLAB function that contains the coeficents for
% a system of two PDE in time and one space dimension.
c = [1; 1];
b = [1; 1] .* DuDx + [2*u(1)+u(2);u(1)+2*u(2)+3*fc.getVal(x)^2*u(1)];
s = [0;0];
end

