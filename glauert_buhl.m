function a= glauert_buhl(a)
k= -a/(a-1);
F= 1;

if k <= 2/3  % momentum state
    a= k/(1+k);
else 
    % Glauert(Buhl) correction
    g1= 2*F*k - (10/9-F);
    g2= 2*F*k - (4/3-F)*F;
    g3= 2*F*k - (25/9-2*F);

    if abs(g3) < 1e-6_dp  % avoid singularity
        a= 1 - 1/2/sqrt(g2);
    else
        a= (g1 - sqrt(g2)) / g3;
    end
end
