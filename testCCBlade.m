old_dir= pwd;
cd('../matlab-toolbox/FAST2MATLAB')
addpath(pwd)
cd('../Utilities')
addpath(pwd)
cd(old_dir)


%%
fst_file= 'AeroData/5MW_Land_DLL_WTurb.fst';

[fst_dir, fname]= fileparts(fst_file);
base_path= fullfile(fst_dir, fname);

%%
fstDataOut = FAST2Matlab(fst_file);

EDFile= strrep(GetFASTPar(fstDataOut, 'EDFile'), '"', '');
edDataOut = FAST2Matlab(fullfile(fst_dir, EDFile));

%%
AeroFile= strrep(GetFASTPar(fstDataOut, 'AeroFile'), '"', '');
adDataOut = FAST2Matlab(fullfile(fst_dir, AeroFile));

ADBldFile= strrep(GetFASTPar(adDataOut, 'ADBlFile(1)'), '"', '');
adbldDataOut = FAST2Matlab(fullfile(fst_dir, ADBldFile));

data= [];
data.R= adbldDataOut.BldNodes(:, find(strcmp(adbldDataOut.BldNodesHdr, 'BlSpn'))) + GetFASTPar(edDataOut, 'HubRad');
data.chord= adbldDataOut.BldNodes(:, find(strcmp(adbldDataOut.BldNodesHdr, 'BlChord')));
data.twist= adbldDataOut.BldNodes(:, find(strcmp(adbldDataOut.BldNodesHdr, 'BlTwist')))/180*pi;
data.airfoil_idx= adbldDataOut.BldNodes(:, find(strcmp(adbldDataOut.BldNodesHdr, 'BlAFID')));
data.rho= GetFASTPar(adDataOut, 'AirDens');
data.TipLoss= double(strcmpi(GetFASTPar(adDataOut, 'TipLoss'), 'true'));
data.HubLoss= double(strcmpi(GetFASTPar(adDataOut, 'HubLoss'), 'true'));
data.TanInd= double(strcmpi(GetFASTPar(adDataOut, 'TanInd'), 'true'));
data.AIDrag= double(strcmpi(GetFASTPar(adDataOut, 'AIDrag'), 'true'));
data.TIDrag= double(strcmpi(GetFASTPar(adDataOut, 'TIDrag'), 'true'));
data.SkewMod= GetFASTPar(adDataOut, 'SkewMod');
data.SkewModFactor= GetFASTPar(adDataOut, 'SkewModFactor');
data.AeroFile= AeroFile;
data.B= GetFASTPar(edDataOut, 'NumBl');
data.IndToler= GetFASTPar(adDataOut, 'IndToler');
if strcmpi(strrep(data.IndToler, '"', ''), 'default')
    data.IndToler= 1e-6;
end

idx_alpha= GetFASTPar(adDataOut, 'InCol_Alfa');
idx_cl= GetFASTPar(adDataOut, 'InCol_Cl');
idx_cd= GetFASTPar(adDataOut, 'InCol_Cd');
% idx_cm= GetFASTPar(adDataOut, 'InCol_Cm');
for i= 1:length(adDataOut.FoilNm)
    AirFoil= FAST2Matlab(fullfile(fst_dir, strrep(adDataOut.FoilNm{i}, '"', '')));
    data.AirFoil(i).alpha= AirFoil.AFCoeff(:, idx_alpha)/180*pi;
    data.AirFoil(i).cl= AirFoil.AFCoeff(:, idx_cl);
    data.AirFoil(i).cd= AirFoil.AFCoeff(:, idx_cd);
end

data.acorr= 0.3;

%%
result= CCBlade(data, 7.55, 0, 10);

%%
makeCCBlade_mex
AeroFields= calcAeroFields(data);

cp= AeroFields.cp;
cp(cp<0)=nan;
surf(AeroFields.lambda, AeroFields.theta, cp)

hold on
plot3(7.55, 0, result.cp, 'ro')
view(61, 60)
