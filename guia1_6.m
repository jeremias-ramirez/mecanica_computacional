## Copyright (C) 2019 jeremias
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {Function File} {@var{retval} =} guia1_6 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

source "funaux.m"

function [T, T1, T2] = guia1_6 (N,p)
	xl = 0;
	xr = 1;
	T0 = 0;
	TN = 1;
	
	dx = (xr-xl) / N;

	a1 = 1  / dx;
	b1 = -2 / dx + p ;
	c1 = 1 /dx - p;

	m1 = getMatrixTriang(N,a1,b1,c1);
	
	a2 = 1  / dx + p/2;
	b2 = -2 / dx ;
	c2 = 1 /dx - p/2;

	m2 = getMatrixTriang(N,a2,b2,c2);


	contor1 = zeros(N-2, 1);
	contor1(1, 1) = -a1*T0;
	contor1(end, 1) = -c1*TN;
	
	contor2 = zeros(N-2, 1);
	contor2(1, 1) = -a2*T0;
	contor2(end, 1) = -c2*TN;


	theta = @(x,p) (1- exp(p .* x)) ./ (1-exp(p));
	
	x = linspace(xl,xr,N)';
	
	T = theta(x,p);
	
	T1 = zeros(N,1);
	T1(1,1) = T0;
	T1(N,1) = TN;

	T2 = zeros(N,1);
	T2(1,1) = T0;
	T2(N,1) = TN;
	
	T1(2:N-1,1) = m1 \ contor1;
	T2(2:N-1,1) = m2 \ contor2;

end


iter = 10:10:100;
e1 = zeros(size(iter),1);
e2 = zeros(size(iter),1);


for i = 10:10:100
	[T, T1, T2] = guia1_6(i,20);
	e1(i/10) = norm(T1-T,2);
	e2(i/10) = norm(T2-T,2);
end
plot(iter,log(e1),"b",iter,log(e2),"r")
eR = abs(e1-e2) ./ e1

