function [U,f,g,stp,info,nfev] = cvsrch_alm(U,tA,N,R,normA,Lambda,C,dims,f,g,stp,s)
%CVSRCH   More-Thuente line search from MINPACK.
%
%   Re-factor of API for use in the Poblano Toolbox, 
%   Daniel M. Dunlavy, August 2008, March 2009
%    
%   Translation of minpack subroutine cvsrch
%   Dianne O'Leary   July 1991


% Additions from D. Dunlavy
      x = cell2vec(U,N,R,dims);

      % Initialize params
      xtol = 1e-15;
      ftol = 1e-4;
      gtol = 1e-2;
      stpmin = 1e-15;
      stpmax = 1e15;
      maxfev = 20;
      
      % move up in case of early termination (before this was initialized
      nfev = 0;
      
% Start of D. O'Leary translation      
      p5 = .5;
      p66 = .66;
      xtrapf = 4;
      info = 0;
      infoc = 1;
%
%     Check the input parameters for errors.
%
      if (stp <= 0.0 || ftol < 0.0 ||  ...
          gtol < 0.0 || xtol < 0.0 || stpmin < 0.0  ...
          || stpmax < stpmin || maxfev <= 0) 
         return
      end
%
%     Compute the initial gradient in the search direction
%     and check that s is a descent direction.
%
      dginit = g'*s;
      if (dginit >= 0.0) 
          return
      end
%
%     Initialize local variables.
%
      brackt = 0;
      stage1 = 1;
      % moved up to initialize before any potential return
      % nfev = 0;
      finit = f;
      dgtest = ftol*dginit;
      width = stpmax - stpmin;
      width1 = 2*width;
      wa = x;
%
%     The variables stx, fx, dgx contain the values of the step, 
%     function, and directional derivative at the best step.
%     The variables sty, fy, dgy contain the value of the step,
%     function, and derivative at the other endpoint of
%     the interval of uncertainty.
%     The variables stp, f, dg contain the values of the step,
%     function, and derivative at the current step.
%
      stx = 0.0;
      fx = finit;
      dgx = dginit;
      sty = 0.0;
      fy = finit;
      dgy = dginit;
%
%     Start of iteration.
%
   while (1)   
%
%        Set the minimum and maximum steps to correspond
%        to the present interval of uncertainty.
%
         if (brackt) 
            stmin = min(stx,sty);
            stmax = max(stx,sty);
         else
            stmin = stx;
            stmax = stp + xtrapf*(stp - stx);
         end 
%
%        Force the step to be within the bounds stpmax and stpmin.
%
         stp = max(stp,stpmin);
         stp = min(stp,stpmax);
%
%        If an unusual termination is to occur then let 
%        stp be the lowest point obtained so far.
%
         if ((brackt && (stp <= stmin || stp >= stmax)) ...
            || nfev >= maxfev-1 || infoc == 0 ...
            || (brackt && stmax-stmin <= xtol*stmax)) 
            stp = stx;
         end
%
%        Evaluate the function and gradient at stp
%        and compute the directional derivative.
%
         x = wa + stp * s;
         U = vec2cell(x,N,R,dims);
         [f,g] = fun_alm(U,tA,N,R,normA,Lambda,C,dims);
         nfev = nfev + 1;
         dg = g' * s;
         ftest1 = finit + stp*dgtest;
%
%        Test for convergence.
%
         if ((brackt && (stp <= stmin || stp >= stmax)) || infoc == 0) 
                  info = 6;
         end
         if (stp == stpmax && f <= ftest1 && dg <= dgtest) 
                  info = 5;
         end
         if (stp == stpmin && (f > ftest1 || dg >= dgtest)) 
                  info = 4;
         end
         if (nfev >= maxfev) 
                  info = 3;
         end
         if (brackt && stmax-stmin <= xtol*stmax) 
                  info = 2;
         end
         if (f <= ftest1 && abs(dg) <= gtol*(-dginit)) 
                  info = 1;
         end
%
%        Check for termination.
%
         if (info ~= 0) 
                  return
         end
%
%        In the first stage we seek a step for which the modified
%        function has a nonpositive value and nonnegative derivative.
%
         if (stage1 && f <= ftest1 && dg >= min(ftol,gtol)*dginit) 
                stage1 = 0;
         end
%
%        A modified function is used to predict the step only if
%        we have not obtained a step for which the modified
%        function has a nonpositive function value and nonnegative 
%        derivative, and if a lower function value has been  
%        obtained but the decrease is not sufficient.
%
         if (stage1 && f <= fx && f > ftest1) 
%
%           Define the modified function and derivative values.
%
            fm = f - stp*dgtest;
            fxm = fx - stx*dgtest;
            fym = fy - sty*dgtest;
            dgm = dg - dgtest;
            dgxm = dgx - dgtest;
            dgym = dgy - dgtest;
% 
%           Call cstep to update the interval of uncertainty 
%           and to compute the new step.
%
            [stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm,brackt,infoc] ...
             = cstep(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm, ...
                     brackt,stmin,stmax);
          
%
%           Reset the function and gradient values for f.
%
            fx = fxm + stx*dgtest;
            fy = fym + sty*dgtest;
            dgx = dgxm + dgtest;
            dgy = dgym + dgtest;
         else
%
%           Call cstep to update the interval of uncertainty 
%           and to compute the new step.
%
            [stx,fx,dgx,sty,fy,dgy,stp,f,dg,brackt,infoc] ...
             = cstep(stx,fx,dgx,sty,fy,dgy,stp,f,dg, ...
                     brackt,stmin,stmax);
          
         end
%
%        Force a sufficient decrease in the size of the
%        interval of uncertainty.
%
         if (brackt) 
            if (abs(sty-stx) >= p66*width1) 
              stp = stx + p5*(sty - stx);
            end
            width1 = width;
            width = abs(sty-stx);
         end
%
%        End of iteration.
%
   end
%
%     Last card of subroutine cvsrch.
%


