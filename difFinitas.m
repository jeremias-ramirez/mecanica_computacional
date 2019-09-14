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
## @deftypefn {Function File} {@var{T} = } difFinitas (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: jeremias <jeremias@debian-jere>
## Created: 2019-09-11
source "funaux.m"

function [T] = difFinitas (xnode, model, cb, et)
	N = length(xnode);

	dx = xnode(2,1) - xnode(1,1);
	k = model.k;
	v = model.v;
	c = model.c;
	G = model.G;

	a = 1 + (v * dx) / (2 * k);
	b = -(2 + c * dx^2 / k);
	c = 1 - (v * dx) / (2 * k);
	
	G = (-dx^2 / k) * G;

	if length(G) == 1
		G = ones(N,1) * G;
	end

	G(1,1) = 0;
	G(N,1) = 0;
		

	m = getMatrixN(N, a, b, c);
	[m, b] = getCondicionBordes(m, dx, model , cb);
	T = m \ (G + b);
	figure(1)
	plot(xnode, T);

	if et == 0
		T = m \ (G + b);
		return	
	end
	TI = zeros(N,1);
	figure(2)
	showNoEstacionario(m, b, G, TI, et, dx, model, xnode)

	#fd = 1 ;
	#fa = 2 ;
	#dt = 0.5 * dx^2  / k * fd;
	#rhoCp = model.rhoCp;
	#
	#dt_rhoCp = dt/rhoCp;
	#
	#T = zeros(N,1);
	#I = diag(ones(N,1));

	#m(2:N-1,:) = dt_rhoCp * m(2:N-1,:) + I(2:N-1,:);
	#G = -dt_rhoCp * G;
	#T =  G + m * T + b;
	#
	#tol = 0.0000001;
	#e = 1;
	#iter = 0;

	#figure(2)
	#while ( e > tol || iter < 1000)	
	#	TNew =  G + m * T;

	#	e = norm(TNew - T, 2);
	#	T = TNew;
	#	plot(xnode, T)
	#	pause(0.05)
		
end

