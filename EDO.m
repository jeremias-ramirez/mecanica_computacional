syms x T DT D2T % Defino variables simbólicas
% Matlab entiende a D como derivada primera y D2 como derivada segunda en
% lenguaje simbólico

% Ecuacion diferencial con v=2, k=1, c=3 y G=4
eq = '2*DT - 2*D2T + 3*T - 4 = 0'; 

condD = '1*DT(0)+2*T(0)=2*3'; % Condicion Robin con k=1, h=1 y Tinf=3
condN = '-1*DT(0)=1'; % Condicion Neumann con k=1 y q=1
condR = 'T(1)=10'; % Condicion Dirichlet

T = dsolve(eq,condD,condR,x); % T es la solución simbólica

xdis = 0:0.01:1;
sol = subs(T,x,xdis); % Evalúo T en los puntos xdis que defino.

plot(xdis,sol)

