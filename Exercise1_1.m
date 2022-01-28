% cG(1) finite element approximation of the two-point BVP,
% {-u''(x) = f(x), 0 < x < 1.
% {u(0) = u(1) = 0. 
clear
tic

% Parameters. 
f = @fcn;
Th = [0, 1/4, 1/2, 3/4, 1];
syms x
v = hatFunction(Th, x);

% Plot of the hat functions. 
hold on
for i = 1:(length(Th) - 2)
    fplot(v(i), [0 1])
end
hold off

% Calculation of stiffness matrix and load vector. 
A = zeros(length(Th) - 2, length(Th) - 2);
for i = 1:length(A)
    for j = 1:length(A)
        if abs(i - j) < 2
            if (i == j)
                h1 = Th(i+1) - Th(i);
                h2 = Th(i+2) - Th(i+1);
                A(i,j) = 1/h1 + 1/h2;
            else
                h2 = Th(i+1) - Th(i);
                A(i,j) = -1/h2;
            end
        end
    end
end
b = zeros(length(Th) - 2, 1);
for i = 2:(length(b)+1)
    hi = Th(i) - Th(i-1);
    firstTerm = int(f(x)*(x - Th(i-1))/hi, Th(i-1), Th(i));
    hiPlusOne = Th(i+1) - Th(i);
    secondTerm = int(f(x)*(Th(i+1) - x)/hiPlusOne, Th(i), Th(i+1));
    b(i-1) = firstTerm + secondTerm;
end

% Finally solve the matrix equation.
xi = linsolve(A, b);
un = xi'.*v;
disp(un)

toc

function f = fcn(x)
    f = 6*x;
end

function phi = hatFunction(Th, x)
    phi = [];
    for i = 1:(length(Th) - 2)
        phi = [phi, triangularPulse(Th(i), (Th(i)+Th(i+2))/2, Th(i+2), x)];
    end
end
 