% cG(2) finite element approximation of the two-point BVP,
% {-u''(x) = f(x), 0 < x < 1.
% {u(0) = u(1) = 0. 
clear
tic

% Parameters. 
f = @fcn;
Th = [0, 1/6, 1/3, 1/2, 2/3, 5/6, 1];
%Th = linspace(0, 1, 11); % Should be odd.
syms x
lambda = LagrangePolynomial(Th, x);

% Calculation of stiffness matrix and load vector. 
A = zeros(length(Th)-2, length(Th)-2);
lambdaDiff = diff(lambda);
for i = 1:length(A)
    for j = 1:length(A)
        A(i,j) = int(lambdaDiff(i+1)*lambdaDiff(j+1), 0, 1);
    end
end
b = zeros(length(Th)-2, 1);
for i = 1:length(b)
    b(i) = int(f(x)*lambda(i+1), 0, 1);
end

% Finally solve the matrix equation.
xi = linsolve(A, b);
un = xi'.*lambda(2:(length(Th)-1));
disp(lambda)
disp(xi)

toc

%% Plot of the Lagrange polynomials. 
hold on
for i = 1:length(Th)
    fplot(lambda(i), [0 1])
end
title('Plot of the Lagrange polynomials.')
legend('\lambda_0(x)', '\lambda_{1/2}(x)', '\lambda_1(x)', '\lambda_{3/2}(x)', '\lambda_2(x)', '\lambda_{5/2}(x)', '\lambda_3(x)')
xlabel('x')
hold off

%% Plot of the solution and the cG(2) approximation. 
hold on
fplot(x - x^3, [0 1], 'LineStyle', '-')
fplot(sum(un), [0 1], 'LineStyle', '--')
title('Plot of the exact solution, u(x), and the approximation cG(2), u_n(x).')
legend('u(x)', 'u_n(x)')
xlabel('x')
hold off

%% Functions.
function f = fcn(x)
    f = 6*x;
end

function lambda = LagrangePolynomial(Th, x)
    lambda = [];
    q = length(Th);
    h = Th(1+2) - Th(1);
    for i = 1:q
        if (i == 1) % Special case.
            int = 2*(Th(i+1) - x)*(Th(i+2) - x)/(Th(i+2) - Th(i))^2;
            temp = piecewise(x < Th(i+2), int, x > Th(i+2), 0);
        elseif (i == q) % Special case.
            int = 2*(x - Th(i-2))*(x - Th(i-1))/(Th(i) - Th(i-2))^2;
            temp = piecewise(x > Th(i-2), int, x < Th(i-2), 0);
        elseif (mod(i, 2) ~= 0) % Whole number. 
            int1 = 2*(x - Th(i-2))*(x - Th(i-1))/(Th(i) - Th(i-2))^2;
            int2 = 2*(Th(i+1) - x)*(Th(i+2) - x)/(Th(i+2) - Th(i))^2;
            temp = piecewise(Th(i-2) < x < Th(i), int1, Th(i) < x < Th(i+2), int2, x < Th(i-2), 0, x > Th(i+2), 0);
        else % Half number. 
            int = 4*(x - Th(i-1))*(Th(i+1) - x)/(2*(Th(i+1) - Th(i)))^2;
            temp = piecewise(Th(i-1) < x < Th(i+1), int, x < Th(i-1), 0, x > Th(i+1), 0);
        end
        lambda = [lambda, temp];
    end
end
