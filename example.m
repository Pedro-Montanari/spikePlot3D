clear
close all
clc
%% Variable definition
% Create dummy sample points
ptsS = 10; 

% Define radius of the polar domains
r0 = 3;

% Define the discretization of the polar domains and of the structured 
% domain. To obtain a good enough discretization leave as below
gridS = 10*r0; 
gridD = ptsS*gridS; 

% Generate dummy x coordinates of sample points
Xmin = 0;
Xmax = 10;
X(:,1) = Xmin + (Xmax-Xmin).*rand(ptsS,1);

% Generate dummy y coordinates sample points
Ymin = 0;
Ymax = 20;
Y(:,1) = Ymin + (Ymax-Ymin).*rand(ptsS,1);

% Generate dummy values of sample points
Vmin = 50;
Vmax = 300;
V(:,1) = Vmin + (Vmax-Vmin).*rand(ptsS,1);

S = [X,Y,V];
[Xdomain,Ydomain,V] = spikePlot3D(S,r0,gridS,gridD);

% Plot
figure;
colormap("turbo");
surf(Xdomain,Ydomain,V,"EdgeColor","none","FaceColor","interp");
view(2);
colorbar;
title("Without domain size specified");
xlim([min(Xdomain,[],"all") max(Xdomain,[],"all")]);
ylim([min(Ydomain,[],"all") max(Ydomain,[],"all")]);
%%
% % Optional: Define the structured domain size. The results are similar 
% to the ones of an unbound case with axis limits specified. 
% The main differences are that:
%   1. If the region of interest (ROI) is only a subset of the domain, 
%   specifying the domain bounds will lead to a finer discretization. 
%   That is because if the domain bounds are not specified gridD-by-gridD 
%   points will be used in the discretization of the structured domain, 
%   while if they are specified gridD-by-gridD points will be used in the 
%   discretization of the ROI only.
%   2. If domain bounds are not set and a ROI is later specified, the 
%   results will not be scaled in relation to the ROI. That is, the maxima
%   and minima will be global, and they can be inside or outside the ROI.
%   If the bounds are specified, the results are scaled based on the ROI,
%   that is, the maxima and minima are the maxima and minima of the ROI.

% Create a dummy ROI
Xlim = [0.25*Xmax;
        0.75*Xmax];
Ylim = [0.25*Ymax;
        0.75*Ymax];

% Plot without domain size specified but with axis limits set
figure;
colormap("turbo");
surf(Xdomain,Ydomain,V,"EdgeColor","none","FaceColor","interp");
view(2);
colorbar;
title("Without domain size specified, with axis lims set");
xlim([Xlim(1) Xlim(2)]);
ylim([Ylim(1) Ylim(2)]);

% Run an bounded case with the previously defined ROI as bounds
D = [Xlim,Ylim];
[Xdomain,Ydomain,V] = spikePlot3D(S,r0,gridS,gridD,D);

% Plot with domain size specified
figure;
c = colormap("turbo");
surf(Xdomain,Ydomain,V,"EdgeColor","none","FaceColor","interp");
view(2);
colorbar;
title("With domain size specified");
xlim([min(Xdomain,[],"all") max(Xdomain,[],"all")]);
ylim([min(Ydomain,[],"all") max(Ydomain,[],"all")]);