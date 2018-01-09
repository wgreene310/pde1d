function howard_example_3
% howard_example_3: MATLAB script that solves the PDE
% stored in eqn2.m, bc2.m, and initial2.m
m = 0;
x = linspace(0,1,10);
t = linspace(0,1,10);
sol = pde1d(m,@eqn2,@initial2,@bc2,x,t);
u1 = sol(:,:,1);
u2 = sol(:,:,2);
figure
subplot(2,1,1)
surf(x,t,u1);
title('u1(x,t)');
xlabel('Distance x');
ylabel('Time t');
subplot(2,1,2)
surf(x,t,u2);
title('u2(x,t)');
xlabel('Distance x');
ylabel('Time t');
end

function [c,b,s] = eqn2(x,t,u,DuDx)
% EQN2: MATLAB function that contains the coeficents for
% a system of two PDE in time and one space dimension.
c = [1; 1];
b = [1; 1] .* DuDx;
s = [u(1)*(1-u(1)-u(2)); u(2)*(1-u(1)-u(2))];
end

function [pl,ql,pr,qr] = bc2(xl,ul,xr,ur,t)
% BC2: MATLAB function that defines boundary conditions
% for a system of two PDE in time and one space dimension.
pl = [0; ul(2)];
ql = [1; 0];
pr = [ur(1)-1; 0];
qr = [0; 1];
end

function value = initial2(x);
% INITIAL2: MATLAB function that defines initial conditions
% for a system of two PDE in time and one space variable.
value = [x^2; x*(x-2)];
end