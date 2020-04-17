function [residual, a, ap, cl, cd]= phiResidual(phi, data, lambda, pitch, node)

[cl, cd]= clcd(phi, data, lambda, pitch, node);

sphi= sin(phi);
cphi= cos(phi);

if data.AIDrag
    cn= cl*cphi;
else
    cn= cl*cphi + cd*sphi;
end

if data.TIDrag
    ct= cl*sphi;
else
    ct= cl*sphi - cd*cphi;
end

Ftip= 1;
if data.TipLoss
    factortip= data.B/2*(data.R(end) - data.R(node))/(data.R(node)*abs(sphi));
    Ftip= 2/pi*acos(exp(-factortip)) + eps;
end

Fhub= 1;
if data.HubLoss
    % factorhub= data.B/2*(data.R(node) - data.R(1))/(Rhub*sphi); error in original code
    factorhub= data.B/2*(data.R(node) - data.R(1))/(data.R(node)*abs(sphi));
    Fhub= 2/pi*acos(exp(-factorhub)) + eps;
end

F= Ftip * Fhub;

k= data.sigma_p(node)*cn/4/F/sphi/sphi;
kp= data.sigma_p(node)*ct/4/F/sphi/cphi;

% compute axial induction factor
if phi > 0 % momentum/empirical
    % update axial induction factor
    if k < -data.acorr/(data.acorr-1)  % momentum state
        % if k <= 2/3  % momentum state
        a= k/(1+k);
    else 
        % Glauert(Buhl) correction
        %          g1= 2*F*k - (10/9-F);
    %          g2= 2*F*k - (4/3-F)*F;
%          g3= 2*F*k - (25/9-2*F);
%  
%          if abs(g3) < 1e-6_dp  % avoid singularity
%              a= 1 - 1/2/sqrt(g2);
%          else
%              a= (g1 - sqrt(g2)) / g3;
%          end

        % Glauert correction
        if k==0
            warning('k==0');
        end

        k_= 1/k;
        k__= (k_*(1-2*data.acorr)+2)^2 + 4*(k_*data.acorr^2-1);
        if k__<-0.01
            warning('k__= %f', k__);
        end
        if k__<0
            k__= eps;
        end
        a= 0.5 * ( 2 + k_ * (1 - 2*data.acorr) - sqrt( k__ ) );
        k= -(a-1)/a;
    end
else  % propeller brake region (a and ap not directly used but update anyway)
    if k > 1
        a= k/(k-1);
    else
        a= 0;  % dummy value
    end
end
        
% compute tangential induction factor
ap= kp/(1-kp);
if ~data.TanInd
    ap= 0;
    kp= 0;
end

% error function
lambda_r= lambda * data.R(node)/data.R(end);
if phi > 0 % momentum/empirical
    residual= sphi/(1+eps-a) - cphi/lambda_r*(1-kp);
else  % propeller brake region
    residual= sphi*(1-k) - cphi/lambda_r*(1-kp);
end

    
