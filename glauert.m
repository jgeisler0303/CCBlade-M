function a= glauert(a)
% Glauert correction
acorr= 0.4;

k= -a/(a-1);

if k < -acorr/(acorr-1) 
    a= k/(1+k);
else 
    if k==0
        warning('k==0');
    end

    k_= 1/k;
    k__= (k_*(1-2*acorr)+2)^2 + 4*(k_*acorr^2-1);
    if k__<-0.01
        warning('k__= %f', k__);
    end
    if k__<0
        k__= eps;
    end
    a= 0.5 * ( 2 + k_ * (1 - 2*acorr) - sqrt( k__ ) );
end
    
