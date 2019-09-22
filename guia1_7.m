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
## @deftypefn {Function File} {@var{retval} =} guia1_7 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: jeremias <jeremias@debian-jere>
## Created: 2019-09-10
source "funaux.m"
source "difFinitas.m"

function [x, T1] = guia1_7 (dx, k, v)
	xl = 0;
	xr = 1;
	T0 = 1;
	TN = 0;
	
	dx = (2 * k) / v 

	x = [xl:dx:xr]';

	N = length(x)	

	a1 = k  / dx^2 + v / (2 * dx);
	b1 = -2 * k / dx^2 ;
	c1 = k  / dx^2 - v / (2 * dx);

	m1 = getMatrixTriang(N,a1,b1,c1)
	
	contor1 = zeros(N-2, 1);
	contor1(1, 1) = -a1*T0;
	contor1(end, 1) = -c1*TN;
	
	
	T1 = zeros(N,1);
	T1(1,1) = T0;
	T1(N,1) = TN;
	T1(2:N-1,1) = m1 \ contor1;
	figure(4)
	plot(x,T1)
	
end

model = struct();
model.k = 1;
model.v = 1;
model.c = 0;
model.rhoCp = 0;

cb = [[1,1,-1];[1,0,-1]];

xnode = [0:0.1:1]';

model.G = 0 .* xnode;
et = 0;

T1 = difFinitas(xnode, model, cb, et);
figure(1)
plot(xnode,T1)


model.k = 0.01;
T2 = difFinitas(xnode, model, cb, et);
figure(2)
plot(xnode,T2)


