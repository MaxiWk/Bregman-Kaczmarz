% find argmin g(t)=f(z-t*a)+t*y where grad(f) is shrinkage with c
% we have g'(t)= m*t-b continuous, increasing and piecewise linear, i.e. slope m and intercept -b vary on different intervals
% find t with 0=g'(t) <=> t=b/m in the appropriate interval
% if g'(0)<0 (>0) then we must have t>0 (<0)
% it is assumed that x=shrinkage(z,c)
% Authors: Frank Schöpfer, Maximilian Winkler
% Original code from https://ieeexplore.ieee.org/abstract/document/7025269/?casa_token=gLoOkpKVjnsAAAAA:8MY-JZXarAfZlcTOnpgyz9u0cm-0OSfLnHLu1AzlJPRnVDOKnT3KPy_NT7JcZo5nELK4TBs
function [vars,t]=linesearch_shrinkage(vars,a,y,c,b)

    x = vars.x;
    z = vars.xstar;

    % first handle trivial cases (hyperplane empty or full space; differentiable objective)
    
    if norm(a) < 1e-11
        if abs(y) < 1e-11
            t = 0; 
        else
             disp('Error in linesearch_shrinkage: a=0, b!=0 -> empty hyperplane!')
        end
             return
    end
    
   % start intervall contains initial value t0=0, i.e. initial b=-g'(0)
   if nargin < 6
        b = a'*x-y;
   end
    
   if abs(c) < 1e-11  % if g is differentiable, done already with explicit solution
        t = b/norm(a)^2;
        s = 1;
   else 
       
        % nontrivial cases       
    
        s = sign(b);
        if s == 0 % <=> g'(0)=0 and we are already done
            t = 0;
            return;
        elseif s == -1
            a = -a; % to ensure that we look for t>0
            b = -b;
        end

        % from now on only entries are needed where a~=0 
        I = (a~=0);
        n = numel(I); 
        z_old = z;
        z = z(I);
        a_old = a;
        a = a(I);
        a_square = a.^2; % will be needed for updates of m and b
        csa = c*sign(a);

        % compute kinks and sort them

        r = (z+csa)./a; % right kinks
        l = (z-csa)./a; % left kinks
        m = sum(a_square(l>0)) + sum(a_square(r<0)); % initial slope m
        % keep track of relevant indices
        I = 1:n; 
        % only kinks >= 0 may be crossed
        Ir = I(r>=0);
        Il = I(l>=0);
        r = r(Ir);
        l = l(Il);
        nl = numel(l); % number of left kinks
        % collect all kinks
        kink = [l;r];
        [kink,ind] = sort(kink); % now the successive intervals are [kink(k),kink(k+1)]
        % keep track of original indices
        I = [Il,Ir];
        I = I(ind);
        N = numel(kink);

        k = 0;
        while (k < N) && (m*kink(k+1) < b) % <=> g'(kink(k+1)) < 0
            k = k+1;
            if ind(k) <= nl % t crosses some left kink
                b = b-kink(k)*a_square(I(k));
                m = m-a_square(I(k));
            else % t crosses some right kink
                b = b+kink(k)*a_square(I(k));
                m = m+a_square(I(k));
            end    
        end
        if (m==0) && (N > 0)
            t = s*kink(max(k,1));
        else
            t = s*b/m;
        end
        
        z = z_old;
        a = a_old;
        
    end
    
    z = z-t*s*a;
    vars.xstar = z;
    vars.x = min(max(0,z-c),z+c); % x=shrinkage(z,c)
    
end

