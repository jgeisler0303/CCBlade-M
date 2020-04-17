function [cl, cd]= clcd(phi, data, lambda, pitch, node)

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

