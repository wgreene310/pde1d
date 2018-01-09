function howard_example_1
%PDE1: MATLAB script M-le that solves and plots
%solutions to the PDE stored in eqn1.m
m = 0;
%NOTE: m=0 speci?es no symmetry in the problem. Taking
%m=1 species cylindrical symmetry, while m=2 speci?es
%spherical symmetry.
%
%De?ne the solution mesh
x = linspace(0,1,20);
t = linspace(0,2,10);
%Solve the PDE
u = pde1d(m,@eqn1,@initial1,@bc1,x,t);
%Plot solution
surf(x,t,u);
title('Surface plot of solution.');
xlabel('Distance x');
ylabel('Time t');
figure;
fig = plot(x,u(1,:),'erasemode','xor');
for k=2:length(t)
pause(.5);
set(fig,'xdata',x,'ydata',u(k,:));
drawnow;
end
end

function [c,b,s] = eqn1(x,t,u,DuDx)
%EQN1: MATLAB function M-?le that species
%a PDE in time and one space dimension.
c = 1;
b = DuDx;
s = 0;
end

function [pl,ql,pr,qr] = bc1(xl,ul,xr,ur,t)
%BC1: MATLAB function M-le that speci?es boundary conditions
%for a PDE in time and one space dimension.
pl = ul;
ql = 0;
pr = ur-1;
qr = 0;
end

function value = initial1(x)
%INITIAL1: MATLAB function M-le that species the initial condition
%for a PDE in time and one space dimension.
value = 2*x/(1+x^2);
end

