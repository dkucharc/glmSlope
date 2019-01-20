function [w,w0,info] = logistic_slope(X,y,lambda,options)
% logistic_slope   Sorted L1 parameter estimation solver for logistic
% regression
%
% [w,info] = logistic_slope(X,y,lambda,options) solves the Slope problem
%
%       Minimize_{w}  sum(log(1+exp(-diag(y)*X*w)) + sum_i (lambda_i * |w|_[i])
%
% where X is an n x p data matrix, y \in {-1,1}^n contains binary
% responses, and |w|_[i] denotes the i-th largest entry in |w|. The entries
% in lambda must be nonnegative and in non-increasing order. When lambda is a
% scalar, the above formulation is reduced to the Lasso:
%
%       Minimize_{w} sum(log(1+exp(-diag(y)*X*w)) + lambda * ||w||_1.
%
% The options parameter is a structure with the following optional fields
% with [default value]: 
%
%    .iterations    Maximum number of iterations                  [10,000]
%    .verbosity     0 = nothing, 1 = major, 2 = every                  [1]
%    .fid           File identifier for output                [1 = stdout]
%    .optimIter     Iterations between optimality-condition checks     [1]
%    .tolInfeas     Maximum allowed dual infeasibility              [1e-6]
%    .tolRelGap     Stopping criterion for relative primal-dual gap [1e-6]
%    .xInit         Initial value of x                        [zeros(n,1)]
%
% The info output structure contains the following fields
%
%    .runtime       Runtime
%    .objPrimal     Primal objective
%    .objDual       Dual objective (possibly for infeasible dual point)
%    .infeas        Dual infeasibility
%    .status        Status: 1 = optimal, 2 = iterations
%

% Copyright 2016, S. Lee

% This file used Adlas.m as a template, from the SLOPE Toolbox version 1.0.
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

% -------------------------------------------------------------
% Start timer
% -------------------------------------------------------------
t0 = tic();


% -------------------------------------------------------------
% Parse parameters
% -------------------------------------------------------------
if (nargin <  4), options = struct(); end;

iterations = getDefaultField(options,'iterations',10000); %10000
verbosity  = getDefaultField(options,'verbosity',1);
fid        = getDefaultField(options,'fid',1);
optimIter  = getDefaultField(options,'optimIter',1);
tolInfeas  = getDefaultField(options,'tolInfeas',1e-6);
tolRelGap  = getDefaultField(options,'tolRelGap',1e-6);
wInit      = getDefaultField(options,'wInit',[]);

% Ensure that lambda is non-increasing
if ((length(lambda) > 1) && any(lambda(2:end) > lambda(1:end-1)))
   error('Lambda must be non-increasing.');
end
if (lambda(end) < 0)
   error('Lambda must be nonnegative');
elseif (lambda(1) == 0)
   error('Lambda must have at least one nonnegative entry.');
end


% -------------------------------------------------------------
% Initialize
% -------------------------------------------------------------

% Get problem dimension
p = size(X,2);

% Get initial lower bound on the Lipschitz constant
s = RandStream('mt19937ar','Seed',0);
x = randn(s,p,1); x = x / norm(x,2);
x = X'*(X*x);
L = norm(x,2)/4;

% Constants for exit status
STATUS_RUNNING    = 0;
STATUS_OPTIMAL    = 1;
STATUS_ITERATIONS = 2;
STATUS_MSG = {'Optimal','Iteration limit reached'};

% Initialize parameters and iterates
if (isempty(wInit)), wInit = zeros(p,1); end;

t       = 1;
eta     = 2;
lambda  = lambda(:);
y       = y(:);
Y       = diag(y);
YX      = Y*X;
w       = wInit;
w0      = 0;
v       = w;
iter    = 0;
status  = STATUS_RUNNING;

% Deal with Lasso case
modeLasso = (numel(lambda) == 1);
if (modeLasso)
   proxFunction = @(v1,v2) proxL1(v1,v2);
else
   proxFunction = @(v1,v2) proxSortedL1(v1,v2);
end

if (verbosity > 0)
   fprintf(fid,'%5s   %9s    %9s   %9s\n','Iter', 'Gap','Infeas.', 'Rel. gap');
end

% -------------------------------------------------------------
% Main loop
% -------------------------------------------------------------
while (true)

   % Compute the gradient at f(y)
   r = 1./(1+exp(YX*v));
   g = -YX'*r;       % full grad [-YX'*r; -y'*r]
   log1mr = log(1-r);
   f = -sum(log1mr);
   
   % Increment iteration count
   iter = iter + 1;

   % Check optimality conditions
   if ((mod(iter,optimIter) == 0))
      % Compute 'dual', check infeasibility and gap
      if (modeLasso)
         infeas = max(abs(g) - lambda, 0);

         objPrimal =  f + lambda*norm(v,1);
         objDual   = (r-1)'*log1mr - r'*log(r);
      else
         vs = sort(abs(v),'descend');
         gs = sort(abs(g),'descend');
         
         % Compute primal and dual objective
         infeas = max(max(cumsum(gs-lambda),0));
         
         objPrimal =  f + lambda'*vs;
         objDual = (r-1)'*log1mr - r'*log(r);
      end
      
      % Format string
      if (verbosity > 0)
         str = sprintf('   %9.2e  %9.2e   %9.2e', objPrimal - objDual, infeas/lambda(1), abs(objPrimal - objDual) / max(1,objPrimal));
      end
      
      % Check primal-dual gap
      if ((abs(objPrimal - objDual)/max(1,objPrimal) < tolRelGap) && ...
         (infeas < tolInfeas * lambda(1)))
         status = STATUS_OPTIMAL;
      end
   else
      str = '';
   end

   if (verbosity > 0)
       if ((verbosity == 2) || ...
               ((verbosity == 1) && (mod(iter,optimIter) == 0)))
           fprintf(fid,'%5d  %s\n', iter,str);
       end
   end
   
   % Stopping criteria
   if (status == 0)
      if (iter >= iterations)
         status = STATUS_ITERATIONS;
      end
   end
   
   if (status ~= 0)
      if (verbosity > 0)
         fprintf(fid,'Exiting with status %d -- %s\n', status, STATUS_MSG{status});
      end
      break;
   end
   

   % Keep copies of previous values
   wPrev  = w;
   fPrev  = f;
   tPrev  = t;
   
   % Lipschitz search
   while (true)
      % Compute prox mapping
      w = proxFunction(v - (1/L)*g, lambda/L);
      d = w - v;
      
      r = 1./(1+exp(YX*w));
      f = -sum(log(1-r));
      q  = fPrev + d'*g + (L/2)*(d'*d);
      
      if (q >= f*(1-1e-12))
         break;
      else
         L = L * eta;
      end
   end

   % Update
   t = (1 + sqrt(1 + 4*t^2)) / 2;
   v = w + ((tPrev - 1) / t) * (w - wPrev);
end

% Solution
w = v;

% Information structure
info = struct();
if (nargout > 1)
   info.runtime   = toc(t0);
   info.objPrimal = objPrimal;
   info.objDual   = objDual;
   info.infeas    = infeas;
   info.status    = status;
   info.L         = L;
end

end 


% ------------------------------------------------------------------------
function opt = getDefaultField(data,field,default)
% ------------------------------------------------------------------------
   if isfield(data,field)
      opt = data.(field);
   else
      opt = default;
   end
end


% ------------------------------------------------------------------------
function x = proxL1(y,lambda)
% ------------------------------------------------------------------------
   % Normalization
   y    = y(:);
   sgn  = sign(y);
   x    = sgn .* max(abs(y) - lambda,0);
end
