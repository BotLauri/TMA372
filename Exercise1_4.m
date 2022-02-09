% b) cG(1) numerical methods for the IVP,
% {u'(t) + a*u(t) = f(t), 0 < t <= T.
% {u(0) = u0. 
clear
tic

% Parameters. 
T = 1; 
n = 50;
h = T/n;
u0 = 1; 
u = zeros(1, n);
u(1) = u0;

t = linspace(0, T, n);
for i = 2:n
    u(i) = (((t(i)^3 - t(i-1)^3)/3 - (u(i-1)*(2*h - 1)))/(1 + 2*h));
end

% c) Schemes for the same problem.
% u'(t) = -a*u(t), 0 < t <= T,
% u(0) = u0 = 1.

% Parameters. 
u0 = 1; T = 1; a = 10;
N = 9;
k = T/N;
uExp = zeros(1, N); uImp = zeros(1, N); uCN = zeros(1, N);
uExp(1) = u0; uImp(1) = u0; uCN(1) = u0;

t = linspace(0, T, N);
% Explicit/forward Euler. 
% Stable oscillation for N = a/2.
% Overshoots if N < a.
for n = 2:N
    uExp(n) = uExp(n-1) + k*(-a)*uExp(n-1);
end

% Implicit/backward Euler.
% y_{n+1} = y_n + k*f(y_{n+1}) -> y_{n+1} = y_n/(1 + a*k).
for n = 2:N
    uImp(n) = uImp(n-1)/(1 + a*k);
end

% Crank-Nicolson. 
% y_{n+1} = y_n + k/2*[f(y_{n+1}) + f(y_n)] -> 
% y_{n+1} = y_n*(1 - a*k/2)/(1 + a*k/2).
for n = 2:N
    uCN(n) = uCN(n-1)*(1 - a*k/2)/(1 + a*k/2);
end

toc

%% b) Plot of the solution and the cG(1) approximation. 
ts = linspace(0, T, n);
syms t
hold on
fplot(31/32*exp(-4*t) + t^2/4 - t/8 + 1/32, [0 T])
plot(ts, u)
title('Plot of the exact solution, u(t), and the approximation cG(1), u_n(t).')
legend('u(t)', 'u_n(t)')
xlabel('t')
hold off

%% c) Plot of the solution and the cG(1) approximation. 
ts = linspace(0, T, N);
syms t
hold on
fplot(exp(-a*t), [0 T])
plot(ts, uExp)
plot(ts, uImp)
plot(ts, uCN)
title('Plot of the exact solution, u(t), and the approximation cG(1), u_n(t).')
legend('u = e^{-at}', 'u_{Exp.}(t)', 'u_{Imp.}(t)', 'u_{CN}(t)')
xlabel('t')
hold off
