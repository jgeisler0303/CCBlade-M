function result= CCBlade(data, lambda, pitch, v_wind)
% ------ BEM solution method see (Ning, doi:10.1002/we.1636) ------

data.sigma_p= data.B/2/pi*data.chord./data.R;
v_rot= lambda*v_wind * data.R/data.R(end);
omega= lambda*v_wind/data.R(end);

phi= zeros(size(data.R));
cl= zeros(size(data.R));
cd= zeros(size(data.R));
v_res= zeros(size(data.R));
a= zeros(size(data.R));
ap= zeros(size(data.R));
rho= data.rho;


epsilon= 1e-6;

for node= 1:length(data.R)
    phi_lower= epsilon;
    phi_upper= pi/2;

    if phiResidual(phi_lower, data, lambda, pitch, node)*phiResidual(phi_upper, data, lambda, pitch, node) > 0  % an uncommon but possible case
        if phiResidual(-pi/4, data, lambda, pitch, node) < 0 && phiResidual(epsilon, data, lambda, pitch, node) > 0
            phi_lower= -pi/4;
            phi_upper= -epsilon;
        else
            phi_lower= pi/2;
            phi_upper= pi - epsilon;
        end
    end    
            
    if phi_upper==-epsilon && phiResidual(-epsilon, data, lambda, pitch, node)<0
        phi(node)= -epsilon;
        iter= -1;
    else
        [phi(node), iter]= brent(@(phi)phiResidual(phi, data, lambda, pitch, node), [phi_lower, phi_upper], data.IndToler);
%          [phi(node), ~, ~, out]= fzero(@(phi)phiResidual(phi, data, lambda, pitch, node), [phi_lower, phi_upper], struct('TolX', data.IndToler)); iter= out.iterations;
    end
    
    % fprintf('node: %d, iter: %d\n', node, iter);

    [~, a(node), ap(node), cl(node), cd(node)]= phiResidual(phi(node), data, lambda, pitch, node);
    
    v_res(node)= sqrt((v_wind.*(1-a(node))).^2 + (v_rot(node).*(1+ap(node))).^2);
end

result.a= a;
result.ap= ap;
result.phi= phi;
result.cl= cl;
result.cd= cd;
result.v_res= v_res;

fl= rho/2*data.chord .* v_res.^2 .* cl;
fd= rho/2*data.chord .* v_res.^2 .* cd;
fax= fl.*cos(phi) + fd.*sin(phi);
ftan= -fd.*cos(phi) + fl.*sin(phi);


dR= diff(data.R);
Fl_sect= (fl(1:end-1)+fl(2:end))/2 .* dR;
Fd_sect= (fd(1:end-1)+fd(2:end))/2 .* dR;

Fl= 0.5*[Fl_sect(1); Fl_sect(1:end-1)+Fl_sect(2:end); Fl_sect(end)];
Fd= 0.5*[Fd_sect(1); Fd_sect(1:end-1)+Fd_sect(2:end); Fd_sect(end)];

result.Fl= Fl;
result.Fd= Fd;

Fax_sect= (fax(1:end-1)+fax(2:end))/2 .* dR;;
Ftan_sect= (ftan(1:end-1)+ftan(2:end))/2 .* dR;

Fax= Fl.*cos(phi) + Fd.*sin(phi);
Ftan= -Fd.*cos(phi) + Fl.*sin(phi);
result.Fax= Fax;
result.Ftan= Ftan;

Fflap= Fax.*cos(pitch) + Ftan.*sin(pitch);
Fedge= -Fax.*sin(pitch) + Ftan.*cos(pitch);
result.Fflap= Fflap;
result.Fedge= Fedge;

Max_sect= int_torque(fax, data.R);
result.Max=  0.5*[Max_sect(1); Max_sect(1:end-1)+Max_sect(2:end); Max_sect(end)];

Mtan_sect= int_torque(ftan, data.R);
Mtan= 0.5*[Mtan_sect(1); Mtan_sect(1:end-1)+Mtan_sect(2:end); Mtan_sect(end)];

result.Mtan= Mtan;

Prot= 3*sum(Mtan)*omega;
ProtVec= 3*Mtan*omega;
Fwind= rho/2 * pi*data.R(end)^2 * v_wind^2;
Pwind= Fwind*v_wind;

result.cp= Prot / Pwind;
result.ct= 3*sum(Fax) / Fwind;

result.cp_i= ProtVec / Pwind;
result.ct_i= 3*Fax / Fwind;



function [residual, a, ap, cl, cd]= phiResidual(phi, data, lambda, pitch, node)

alpha= phi-pitch-data.twist(node);
if alpha>pi || alpha<-pi
    warning('alpha not in [-pi pi]')
end

idx_af= data.airfoil_idx(node);
if alpha>data.AirFoil(idx_af).alpha(end)
    cl= data.AirFoil(idx_af).cl(end);
    cd= data.AirFoil(idx_af).cd(end);
elseif alpha<data.AirFoil(idx_af).alpha(1)
    cl= data.AirFoil(idx_af).cl(1);
    cd= data.AirFoil(idx_af).cd(1);
else
    cl= interp1(data.AirFoil(idx_af).alpha, data.AirFoil(idx_af).cl, alpha);
    cd= interp1(data.AirFoil(idx_af).alpha, data.AirFoil(idx_af).cd, alpha);
end


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

function i= int_torque(a, r)
a0= a(1:end-1);
a1= a(2:end);
dR= diff(r);
R0= r(1:end-1);

i= ((2*a1+a0).*dR.^2 + (3*R0.*a1 + 3*R0.*a0).*dR)/6;
