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
	
	[M, F] = getSystem(N, dx, model, cb);


	if et == 0
		T = M \ F;
		return	
	end
	
	TI = zeros(N,1);
	if cb(1,1) == 1
		TI(1,1) = cb(1,2);
	end

	if cb(2,1) == 1
		TI(N,1) = cb(2,2);
	end

	figure(2)
	
	%TI = 100 - 50 * xnode;

	showNoEstacionario(M, F, TI, et, dx, model, xnode)

end

