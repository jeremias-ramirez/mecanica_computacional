function [ht,u1,u2,u3] = Ejemplo_NoEstacionario(hx)
%  Calculamos por diferencias finitas la solucion de una PDE
%       (du/dt)-d^2u/dx^2 = 100*x; 0 < x < 1;
%       u(0,t) = 1; u(1,t) = 0;
%       u(x,0) = 0;
   u_left = 1;   %  Izquierda BC
   u_right = 0;  %  Derecha BC
   x = [0:hx:1]';            %  puntos de la malla, incluyendo la frontera.
   n = length(x);

%  En base al numero de Fourier calculamos un delta t correcto
   %fd=0.5; %factor de disminucion de ht, si fd=1 entonces ht es el maximo posible
   fd=1;
   ht =(0.5*(hx^2))*fd;

%   Definimos un numero maximo de iteraciones temporales nt.
%   Por lo tanto T maximo = ht*nt, si llega a su estado estacionario antes,
%   informamos en que tiempo fue, tomando una tolerancia definida por tol.
   tol = 0.0001; 
   nt = 100; T = nt*ht;
   u = zeros(n,nt); 
   u(1,:) = u_left; %  Izquierda BC
   u(n,:) = u_right; %  Derecha BC
   
   disp('Metodo de integracion explicito');
   no_est = 1; % bandera para saber si llego a un estado estacionario
   a = ht/(hx^2);
   b = 1-((2*ht)/(hx^2));
   figure(1);
   for t=2:nt
       for i=2:(n-1)
           u(i,t)=a*u(i-1,t-1)+b*u(i,t-1)+a*u(i+1,t-1)+(100*x(i)*ht);
       end
       figure(1)
       plot(x,u(:,t),'*r--');
       xlabel('eje x')
       ylabel('u(x)')
       titulo = ['Num. de puntos:  ', num2str(n)];
       title(titulo)
       legend('Método Explicito')
       pause(0.15);
       if (norm(u(:,t)-u(:,t-1),2)/norm(u(:,t),2) < tol)
           disp ('La solucion ha llegado a un estado estacionario');
           fprintf('Los pasos realizados fueron: %d\n',t)
           fprintf('El tiempo transcurrido fue de: %f segundos\n',t*ht)
           no_est = 0;
           break;
       end
   end
   if (no_est)
       disp ('La solucion no llego a un estado estacionario, se requieren mas pasos de simulacion');
   end
   u1=u(:,t);
   
   disp ('Presione una tecla para continuar o espere');
   pause(0.7)
   
   disp('Metodo de integracion implicito');
   fa=2; % factor de aumento de ht
   ht=ht*fa; % aumento el paso de tiempo para acelerar la convergencia 
             % no hay restricciones sobre ht
   no_est = 1; % bandera para saber si llego a un estado estacionario
   u = zeros(n,nt);
   u_ap = zeros(n-2,1);
   K = zeros(n-2);
   Q = zeros (n-2,1);
   A = (1/ht)+(2/(hx^2));
   fila = [(-1/(hx^2)) A (-1/(hx^2))];
   K(1,1:2) = fila(2:3); Q(1) = 100*x(2);
   for i=2:(n-3)
        K(i,i-1:i+1) = fila;
        Q(i) = 100*x(i+1);
   end
   K(n-2,n-3:n-2) = fila(1:2); Q(n-2) = 100*x(n-1);
   for t=2:nt
       b = Q + (u_ap./ht);
       if (t ~= 2)
           b(1) = b(1) + (u_left/(hx^2));
           b(n-2) = b(n-2) + (u_right/(hx^2)); 
       end
       u_ap = K\b;
       u(2:n-1,t) = u_ap;
       u(1,t) = 1;
       u(n,t) = 0;
       figure(2)
       plot(x,u(:,t),'*g--');
       xlabel('eje x')
       ylabel('u(x)')
       titulo = ['Num. de puntos:  ', num2str(n)];
       title(titulo)
       legend('Método Implicito')  
       pause(0.15);
       if (norm(u(:,t)-u(:,t-1),2)/norm(u(:,t),2) < tol)
           disp ('La solucion ha llegado a un estado estacionario');
           fprintf('Los pasos realizados fueron: %d\n',t)
           fprintf('El tiempo transcurrido fue de: %f segundos\n',t*ht)
           no_est=0;
           break;
       end
   end
   if (no_est)
       disp ('La solucion no llego a un estado estacionario, se requieren mas pasos de simulacion');
   end
   u2=u(:,t);
   
   disp ('Presione una tecla para continuar o espere');
   pause(0.7)
   
   disp('Metodo de integracion de Crank Nicholson'); 
   ht=ht*fa; % aumento el paso de tiempo para acelerar la convergenicia
   no_est = 1; % bandera para saber si llego a un estado estacionario
   u = zeros(n,nt);
   u_ap = zeros(n-2,1);
   K = zeros(n-2);
   Q = zeros(n-2,1);
   A = (1/ht)+(1/(hx^2));
   fila = [(-1/(2*(hx^2))) A (-1/(2*(hx^2)))];
   K(1,1:2) = fila(2:3); Q(1) = 100*x(2); 
   for i=2:(n-3)
        K(i,i-1:i+1) = fila;
        Q(i) = 100*x(i+1);
   end
   K(n-2,n-3:n-2) = fila(1:2); Q(n-2) = 100*x(n-1);
   b = zeros (n-2,1);
   u_left = 1;   %  Izquierda BC
   u_right = 0;  %  Derecha BC
   for t=2:nt
       for i=2:(n-1)
           b(i-1) = Q(i-1) + ((1/(2*(hx^2)))*u(i-1,t-1)) + (((1/ht)-(1/(hx^2)))*u(i,t-1)) + ((1/(2*(hx^2)))*u(i+1,t-1));
       end
       if (t ~= 2)
           b(1) = b(1) + ((1/(2*(hx^2)))*u_left);
           b(n-2) = b(n-2) + ((1/(2*(hx^2)))*u_right);
       end
       u_ap = K\b;
       u(2:n-1,t) = u_ap;
       u(1,t) = 1;
       u(n,t) = 0;
       figure(3)
       plot(x,u(:,t),'*b--');
       xlabel('eje x')
       ylabel('u(x)')
       titulo = ['Num. de puntos:  ', num2str(n)];
       title(titulo)
       legend('Crank Nicholson')
       pause(0.15);
       if (norm(u(:,t)-u(:,t-1),2)/norm(u(:,t),2) < tol)
           disp ('La solucion ha llegado a un estado estacionario');
           fprintf('Los pasos realizados fueron: %d\n',t)
           fprintf('El tiempo transcurrido fue de: %f segundos\n',t*ht)
           no_est=0;
           break;
       end
   end
   if (no_est)
       disp ('La solucion no llego a un estado estacionario, se requieren mas pasos de simulacion');
   end
   u3=u(:,t);
   
   disp ('Presione una tecla para continuar o espere');
   pause(0.7)
   
   %analitica
   u_an=(-50/3)*(x.^3)+(47/3)*x+1;
   figure(4);
   plot(x,u1,'+r--',x,u2,'*g-',x,u3,'b--',x,u_an,'oy-')
   xlabel('eje x')
   ylabel('u(x)')
   titulo = ['Num. de puntos:  ', num2str(n)];
   title(titulo)
   legend('Metodo Explicito','Metodo implicito','Crank Nicholson','Analitica')
end

