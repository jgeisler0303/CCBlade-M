function AeroFields= calcAeroFields(data, LambdaVec, ThetaVec)
fprintf('Calculating aero fields. Please be patient.\n');

if ~exist('LambdaVec', 'var')
    LambdaVec= 1:0.5:13;
end
if ~exist('ThetaVec', 'var')
    ThetaVec= 0:0.5:45;
end

AeroFields.lambda= LambdaVec;
AeroFields.theta= ThetaVec;
AeroFields.R= data.R;

n_lam= length(LambdaVec);
n_th= length(ThetaVec);
n_R= length(data.R);

AeroFields.cp= nan(n_th, n_lam);
AeroFields.ct= nan(n_th, n_lam);

for idxTheta= 1:n_th
    theta= ThetaVec(idxTheta);
    for idxLambda= 1:n_lam
        lambda= LambdaVec(idxLambda);
        v_wind= interp1([0 4 10 100], [40 14 5 0], lambda);
        omega= lambda*v_wind/data.R(end);
    
        result= CCBlade_mex(data, lambda, theta/180*pi, v_wind);
    
        AeroFields.cp(idxTheta, idxLambda)= result.cp;
        AeroFields.ct(idxTheta, idxLambda)= result.ct;        
    end
end
