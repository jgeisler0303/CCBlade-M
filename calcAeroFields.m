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
AeroFields.cs= nan(n_th, n_lam);
AeroFields.cpi= nan(n_th, n_lam, n_R);
AeroFields.cti= nan(n_th, n_lam, n_R);
AeroFields.csi= nan(n_th, n_lam, n_R);
AeroFields.cbxi= nan(n_th, n_lam, n_R);
AeroFields.cbyi= nan(n_th, n_lam, n_R);


for idxTheta= 1:n_th
    theta= ThetaVec(idxTheta);
    for idxLambda= 1:n_lam
        lambda= LambdaVec(idxLambda);
        v_wind= interp1([0 4 10 100], [40 14 5 0], lambda);
        omega= lambda*v_wind/data.R(end);
    
        result= CCBlade_mex(data, lambda, theta/180*pi, v_wind);
    
        AeroFields.cp(idxTheta, idxLambda)= result.cp;
        AeroFields.ct(idxTheta, idxLambda)= result.ct;        
        AeroFields.cs(idxTheta, idxLambda)= result.cs;
        AeroFields.cpi(idxTheta, idxLambda, :)= result.cp_i;
        AeroFields.cti(idxTheta, idxLambda, :)= result.ct_i;
        AeroFields.csi(idxTheta, idxLambda, :)= result.cs_i;
        AeroFields.cbxi(idxTheta, idxLambda, :)= result.cbx_i;
        AeroFields.cbyi(idxTheta, idxLambda, :)= result.cby_i;
    end
end

AeroFields.cm= AeroFields.cp ./ repmat(AeroFields.lambda(:)', length(AeroFields.theta), 1);
AeroFields.cmi= AeroFields.cpi ./ repmat(AeroFields.lambda(:)', length(AeroFields.theta), 1, length(AeroFields.R));

AeroFields.lambda_offset= AeroFields.lambda(1);
AeroFields.theta_offset= AeroFields.theta(1);
AeroFields.lambda_scale= 1/diff(AeroFields.lambda(1:2));
AeroFields.theta_scale= 1/diff(AeroFields.theta(1:2));

if isfield(data, 'ModalShapes')
    AeroFields.ModalShapes= data.ModalShapes;
    
    for i= 1:length(AeroFields.ModalShapes)
        for idxTheta= 1:n_th
            for idxLambda= 1:n_lam
                AeroFields.(sprintf('cb%di', i))(idxTheta, idxLambda, :)= sum([squeeze(AeroFields.cbxi(idxTheta, idxLambda, :)) squeeze(AeroFields.cbyi(idxTheta, idxLambda, :))] .* AeroFields.ModalShapes{i}, 2);
                AeroFields.(sprintf('cb%d', i))(idxTheta, idxLambda)= sum(AeroFields.(sprintf('cb%di', i))(idxTheta, idxLambda, :));
            end
        end
    end
    
    for i= 1:length(AeroFields.ModalShapes)
        AeroFields= modalDeriv(AeroFields, i, 'cm');
        AeroFields= modalDeriv(AeroFields, i, 'ct');
        AeroFields= modalDeriv(AeroFields, i, 'cs');
        for j= 1:length(AeroFields.ModalShapes)
            AeroFields= modalDeriv(AeroFields, i, sprintf('cb%d', j));
        end
    end
end


function AeroFields= modalDeriv(AeroFields, n_shape, sensor)

sensor_i= [sensor 'i'];

cXi_= AeroFields.(sensor_i);
ni= size(cXi_, 3);
dcXi_dlam= zeros(size(cXi_));
for i= 1:ni
    dcXi_dlam(:, :, i)= gradient(cXi_(:, :, i), AeroFields.lambda, AeroFields.theta);
end

nth= size(dcXi_dlam, 1);
nlam= size(dcXi_dlam, 2);
nseg= size(dcXi_dlam, 3);

LAM= repmat(AeroFields.lambda(:)', nth, 1, nseg);
R= AeroFields.R(:);
RR= repmat(permute(R(:), [3 2 1]), nth, nlam, 1);
COS= repmat(cosd(AeroFields.theta(:)), 1, nlam, nseg);
SIN= repmat(sind(AeroFields.theta(:)), 1, nlam, nseg);

% vx = axial positive down wind
dcXi_dvx_v= dcXi_dlam .* -LAM; % dcXi_dvx_v= dcXi_dvx * v;

% vy = tangential positiv in direction of positive rotation
dcXi_dom_v= dcXi_dlam * R(end); % dcXi_dom_v= dcXi_dom * v;
dcXi_dvy_v= dcXi_dom_v ./ RR;

% verified in cx_qe_qf.wxmx
ModalShapeX= repmat(permute(AeroFields.ModalShapes{n_shape}(:, 1), [3 2 1]), nth, nlam, 1);
ModalShapeY= repmat(permute(AeroFields.ModalShapes{n_shape}(:, 2), [3 2 1]), nth, nlam, 1);

dcXi_dvbx_v= (-dcXi_dvx_v .* COS + dcXi_dvy_v .* SIN - 2*COS.*AeroFields.(sensor_i));
dcXi_dvby_v= (-dcXi_dvx_v .* SIN - dcXi_dvy_v .* COS - 2*SIN.*AeroFields.(sensor_i));

AeroFields.(sprintf('d%s_dvb%d_v', sensor, n_shape))= sum(dcXi_dvbx_v.*ModalShapeX, 3) + sum(dcXi_dvby_v.*ModalShapeY, 3);

AeroFields.(sprintf('d%si_dvx_v', sensor))= dcXi_dvx_v;
AeroFields.(sprintf('d%si_dvy_v', sensor))= dcXi_dvy_v;
AeroFields.(sprintf('d%si_dom_v', sensor))= dcXi_dom_v;
AeroFields.(sprintf('d%si_dvbx_v', sensor))= dcXi_dvbx_v;
AeroFields.(sprintf('d%si_dvby_v', sensor))= dcXi_dvby_v;

