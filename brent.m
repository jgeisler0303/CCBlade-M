function [xs, iter]= brent(f, ab, t)
% copied from https://de.wikipedia.org/wiki/Brent-Verfahren#Algorithmus_von_Brent_f%C3%BCr_Matlab

if ~exist('t', 'var')
    t= 1e-6;
end

a= ab(1);
b= ab(2);

fa= f(a); fb= f(b);

if fa*fb>0
    error('f(a) and f(b) must have different sign');
end

c= a; fc= fa;

c= a; fc= fa; d= b-a; e= d;

iter= 0;
maxiter= 1000;

while iter<maxiter
    iter= iter+1;
    
    if fb*fc>0
        c= a; fc= fa; d= b-a; e= d;
    end
    
    if abs(fc)<abs(fb)
        a= b; b= c; c= a;
        fa= fb; fb= fc; fc= fa;
    end

    tol= 2*eps*abs(b)+t; m= (c-b)/2;

    if (abs(m)>tol) && (abs(fb)>0)
    
        if (abs(e)<tol) || (abs(fa)<=abs(fb))
            d= m; e= m;
        else
            s= fb/fa;
            if a==c
                p= 2*m*s; q= 1-s;
            else
                q= fa/fc; r= fb/fc;
                p= s*(2*m*q*(q-r)-(b-a)*(r-1));
                q= (q-1)*(r-1)*(s-1);
            end
            if p>0
                q= -q;
            else
                p= -p;
            end
            s= e; e= d;
            if ( 2*p<3*m*q-abs(tol*q) ) && (p<abs(s*q/2))
                d= p/q;
            else
                d= m; e= m;
            end
        end

        a= b; fa= fb;

        if abs(d)>tol
            b= b+d;
        else
            if m>0
                b= b+tol;
            else
                b= b-tol;
            end
        end
    else
        break;
    end
        
    fb= f(b);
end

if iter>=maxiter
    error('Maximum number of iterations reached without convergence');
end
        
xs= b;
