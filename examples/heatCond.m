## Copyright (C) 2016 bgreene
## 
## Author: bgreene
## Created: 2016-01-06

function heatCond

x = linspace(0,1,10);
t = linspace(0,.05,40);

pdeFunc = @(x,t,u,DuDx) heatpde(x,t,u,DuDx);
icFunc = @(x) heatic(x);
bcFunc = @(xl,ul,xr,ur,t) heatbc(xl,ul,xr,ur,t);

m=0;
u = pde1d(m, pdeFunc,icFunc,bcFunc,x,t);

figure; plot(t, u(:,end)); grid on;
xlabel('Time'); ylabel('Temperature');
title('Temperature at right end as a function of time');
figure; plot(x, u(end,:)); grid on;
xlabel('x'); ylabel('Temperature');
title('Temperature along the length at final time');

end

function [c,f,s] = heatpde(x,t,u,DuDx)
c = 1;
f = 10*DuDx;
s = 0;
end

function u0 = heatic(x)
u0 = 0;
end

% --------------------------------------------------------------
function [pl,ql,pr,qr] = heatbc(xl,ul,xr,ur,t)
pl = ul-10;
ql = 0;
pr = 0;
qr = 1;
end
