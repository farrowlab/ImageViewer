H.filepath = [];
H.ZStack = [];
H.zproj = [];
H.roi = [];
H.Max = [];
H.Min = [];
H.CM = colormap(CubeHelix(256,0.5,-1.5,1.2,1.0));    % Color map for figures.
load('roicolours.mat');
H.roiCM = roicolours;
H.npools = 10;
if isfield(H,'poolobj') && isempty(H.poolobj)   
    parpool(H.npools,'nocreate'); 
end
H.poolobj = gcp('nocreate');