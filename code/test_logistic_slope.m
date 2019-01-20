% Copyright 2013, M. Bogdan, E. van den Berg, W. Su, and E.J. Candes

% This file is part of SLOPE Toolbox version 1.0.
%
%    The SLOPE Toolbox is free software: you can redistribute it
%    and/or  modify it under the terms of the GNU General Public License
%    as published by the Free Software Foundation, either version 3 of
%    the License, or (at your option) any later version.
%
%    The SLOPE Toolbox is distributed in the hope that it will
%    be useful, but WITHOUT ANY WARRANTY; without even the implied
%    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%    See the GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with the SLOPE Toolbox. If not, see
%    <http://www.gnu.org/licenses/>.

% n = 10;
% p = 1000;
% 
% X = [-1.5 + 1.5*randn([n/2, p]); 2.6 + .8*randn([n/2,p])];
% y = [ones(n/2,1); -ones(n/2,1)];
% 
% lambda = .1;
% 
% lambda = linspace(1,1e-3,p);

X      = load('data/testData_A.txt');
y      = load('data/testData_b.txt');
lambda = load('data/testData_lambda.txt');

% lambda = median(lambda);
% lambda = lambda / 100;

y_med = median(y);
ind = y>y_med;
y(ind) = 1;
y(~ind) = -1;

[w,w0,info] = logistic_slope(X,y,lambda);

