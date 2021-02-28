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

Fax_sect= (fax(1:end-1)+fax(2:end))/2 .* dR;
Ftan_sect= (ftan(1:end-1)+ftan(2:end))/2 .* dR;

Fax= Fl.*cos(phi) + Fd.*sin(phi);
Ftan= -Fd.*cos(phi) + Fl.*sin(phi);
result.Fax= Fax;
result.Ftan= Ftan;

FBx= Fax.*cos(pitch) + Ftan.*sin(pitch);
FBy= Fax.*sin(pitch) - Ftan.*cos(pitch);
result.FBx= FBx;
result.FBy= FBy;

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
result.cs= sum(Ftan) / Fwind;

result.cp_i= ProtVec / Pwind;
result.ct_i= 3*Fax / Fwind;
result.cs_i= Ftan / Fwind;
result.cbx_i= FBx / Fwind;
result.cby_i= FBy / Fwind;
