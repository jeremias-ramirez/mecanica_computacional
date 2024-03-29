E = 70e9;
A = 0.0003;
l = 3;
P = 50e3;
a = A*E / l;
b = (A * E) / (sqrt(2)*2*l);
c = (A * E) / (2*l);
m = [b+c, b, -b, -c; b, b+a, -b, 0; -b, -b, b+a, 0; -c, 0, 0, b+c];
f = [0; 0; -P; -P]
u=m\f

mr = [0 0 0 -b;0 0 -a -b; -b -b b 0; -b, -b, a+b, 0; 0, -a, 0, b; -c, 0, 0, b+c];
fv = mr * u

f1 = a * [0 -1 0 1] * [0 0 u(2) u(3)]'
f2 = a * [-1 0 1 0] * [0 0 0 u(1)]'

f3 = a * [-1 0 1 0] * [u(2) u(3) 0 u(4)]'

d = (A * E) / (2*l);
f4 = d * [1 -1 -1 1] * [0 u(1) u(2) u(3)]'

f5 = d *  [-1 -1 1 1] * [0 0 0 u(4)]'

f6 = a * [0 -1 0 1] * [0 u(1) 0 u(4)]'

fx3 = f3  - cos(135*pi/180) * f4
fy3 = +P + f1 + sin(135*pi/180) * f4

fx4 = -f3 + f3 +f5 * cos(pi/4) - f5 * cos(pi/4)
fy4 = 2*P +f6 +2 * f5 * sin(pi/4)



