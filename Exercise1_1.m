% cG(1) finite element approximation of the two-point BVP,
% {-u''(x) = f(x), 0 < x < 1.
% {u(0) = u(1) = 0. 
clear
tic

% Parameters. 
f = @fcn;
Th = [0, 1/6, 1/3, 1/2, 2/3, 5/6, 1];
%Th = [0, 1/4, 1/2, 3/4, 1];
syms x
v = HatFunction(Th, x);

% Calculation of stiffness matrix and load vector. 
A = zeros(length(Th)-2, length(Th)-2);
vDiff = diff(v);
for i = 1:length(A)
    for j = 1:length(A)
        A(i,j) = int(vDiff(i)*vDiff(j), 0, Th(end));
    end
end
b = zeros(length(Th)-2, 1);
for i = 1:length(b)
    b(i) = int(f(x)*v(i), 0, Th(end));
end

% Finally solve the matrix equation.
xi = linsolve(A, b);
un = xi'.*v;
disp(un)
toc

%% Plot of the hat functions. 
hold on
for i = 1:(length(Th) - 2)
    fplot(v(i), [0 1])
end
title('Plot of the hat functions.')
legend('\phi_1(x)', '\phi_2(x)', '\phi_3(x)', '\phi_4(x)', '\phi_5(x)')
xlabel('x')
hold off

%% Plot of the solution and the cG(1) approximation. 
hold on
fplot(x - x^3, [0 1])
fplot(sum(un), [0 1])
title('Plot of the exact solution, u(x), and the approximation cG(1), u_n(x).')
legend('u(x)', 'u_n(x)')
xlabel('x')
hold off

%% Functions.
function f = fcn(x)
    f = 6*x;
end

function phi = HatFunction(Th, x)
    phi = [];
    for i = 1:(length(Th) - 2)
        phi = [phi, triangularPulse(Th(i), (Th(i)+Th(i+2))/2, Th(i+2), x)];
    end
end
 