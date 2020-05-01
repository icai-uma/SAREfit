%% Demo of ellipse fitting with outliers SAREfit
%
% Please, cite this work as:
% 
% Karl Thurnhofer-Hemsi, Ezequiel López-Rubio, Elidia Beatriz Blázquez-Parra, M. Carmen Ladrón-de-Guevara-Muñoz, Óscar David de-Cózar-MacÃas,
% Ellipse fitting by spatial averaging of random ensembles,
% Pattern Recognition, 2020, 107406, ISSN 0031-3203,
% https://doi.org/10.1016/j.patcog.2020.107406.
% (http://www.sciencedirect.com/science/article/pii/S0031320320302090)
%
% Last modification: 02/05/2020

clear all
warning off

% Prepare paths
addpath('./competitors');
addpath('./utils');
addpath('./L1mediancov/');

% Set seed
rng('default');

% Save results flag
saveV = false;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA GENERATION (SIMPLE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters of ellipse generation
NumTestSamples = 1000;      % Number of test samples on the true ellipse used for testing
NumTrainSamples=50;         % Number of training samples
NoiseLevel=0.01;            % Noise level
OutlierProbability=0.15;    % Probability of outliers
OcclusionLevel=0;           % Level of occlusion
scale = 1;                  % Scale of the points (default in [0 1]x[0 1])
[TX,X,TrueParA,TrueParG,TrueParN]=GenerateRandomTestTrainingEllipse(NumTestSamples,NumTrainSamples,NoiseLevel,OutlierProbability,OcclusionLevel,scale);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN FITTING METHODS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAREfit
SubsamplingFactor=0.15; % Proportion of the total training samples considered for training
EnsembleSize=90;        % Number of subsamplings
percentile=10;          % Percentage(%) of best fits for the final ensemble
[ParA1,ParG1,ParN1]=SAREfit(X,SubsamplingFactor,EnsembleSize,percentile);


%% Comparison with some methods
% Munoz (2014)
[ParA2,ParG2,ParN2]=EllipseFitMunoz(X);

% Fitzgibbon (1999)
A = EllipseDirectFit(X'); 
ParA3=[A(1) 0.5*A(2) A(3) 0.5*A(4) 0.5*A(5) A(6)]';
[ParG3,code3]=AtoG(ParA3);
ParN3=GtoN(ParG3);

% Taubin (1991)
A = EllipseFitByTaubin(X');
ParA4=[A(1) 0.5*A(2) A(3) 0.5*A(4) 0.5*A(5) A(6)]';
[ParG4,code4]=AtoG(ParA4);
ParN4=GtoN(ParG4);
if code4~=1
    ParG4 = NaN(5,1);
    ParA4 = NaN(6,1);
    ParN4 = NaN(5,1);
end

% Halii and Flusser (1998)
A = EllipseFitHalir(X(1,:), X(2,:)); 
ParA5=[A(1) 0.5*A(2) A(3) 0.5*A(4) 0.5*A(5) A(6)]';
ParG5=AtoG(ParA5);
ParN5=GtoN(ParG5);

% Rosin method, normalization A+C=1
[ParA6,ParG6,ParN6,code6]=EllipseFitRosin(X);
if (code6~=1)
    ParG6 = NaN(5,1);
    ParA6 = NaN(6,1);
    ParN6 = NaN(5,1);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
set(gcf,'Units','normalized');
set(gcf,'Position',[0.25,0.15,0.50,0.70]);

% Select a palette of colors
hold on
MyColors=distinguishable_colors(7);
MyColors(4,:)=MyColors(4,:)+0.5;
% Make yellow as True Ellipse
aux = MyColors(6,:);
MyColors(6,:) = MyColors(1,:);
MyColors(1,:) = aux;
% Switch color orders for the sake of clarity
aux = MyColors(6,:);
MyColors(6,:) = MyColors(5,:);
MyColors(5,:) = aux;

% Plot the ellipse, given the algebraic parameters
Handles=zeros(1,7);
if exist('TrueParG','var')
    MyHandle=PlotEllipseG(TrueParG,MyColors(1,:),4);
    Handles(1)=MyHandle(1);
end
MyHandle=PlotEllipseG(ParG2,MyColors(3,:));
Handles(3)=MyHandle(1);
MyHandle=PlotEllipseG(ParG3,MyColors(4,:));
Handles(4)=MyHandle(1);
MyHandle=PlotEllipseG(ParG4,MyColors(5,:));
Handles(5)=MyHandle(1);
MyHandle=PlotEllipseG(ParG5,MyColors(6,:));
Handles(6)=MyHandle(1);
MyHandle=PlotEllipseG(ParG6,MyColors(7,:));
Handles(7)=MyHandle(1);
% SAREfit is drawn at last
MyHandle=PlotEllipseG(ParG1,MyColors(2,:));
Handles(2)=MyHandle(1);

% Plot Training Samples
hold on
plot(X(1,:), X(2, :), '.r', 'Color', [0 0 0], 'MarkerSize',10);
xlabel('x');
ylabel('y');
axis equal
axis([-1 2.5 -1 2])
grid on

% Save plot
if saveV
    PdfFileName=sprintf('./Demo_synthetic_rng%s',num2str(rngNum));
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperOrientation','portrait');
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperSize',[8 7]);
    set(gcf,'PaperPosition',[0 0 8 7]);
    set(gca,'fontsize',10);
    saveas(gcf,PdfFileName,'pdf');
else
    legend(Handles,'True','SAREfit','Muñoz','Fitzgibbon','Taubin','Halir&Flusser','Rosin','Location','southoutside','Orientation','horizontal');
end

% Compute Natural errors
fprintf('Error for the new algorithm: %f\n',EllipseNaturalError(TrueParN,ParN1));
fprintf('Error for Munoz (2014) algorithm: %f\n',EllipseNaturalError(TrueParN,ParN2));
fprintf('Error for Fitzgibbon (1999) algorithm: %f\n',EllipseNaturalError(TrueParN,ParN3));
fprintf('Error for Taubin (1991) algorithm: %f\n',EllipseNaturalError(TrueParN,ParN4));
fprintf('Error for Halii and Flusser (1998) algorithm: %f\n',EllipseNaturalError(TrueParN,ParN5));
fprintf('Error for Rosin algorithm: %f\n',EllipseNaturalError(TrueParN,ParN6));

% Compute ParG errors
fprintf('Error for the new algorithm: Center=%f, Angle=%f, MajorA=%f, MinorA=%f, Area=%f\n',...
    EllipseParGErrors(TrueParG,ParG1));
fprintf('Error for the Munoz algorithm: Center=%f, Angle=%f, MajorA=%f, MinorA=%f, Area=%f\n',...
    EllipseParGErrors(TrueParG,ParG2));
fprintf('Error for the Fitzgibbon algorithm: Center=%f, Angle=%f, MajorA=%f, MinorA=%f, Area=%f\n',...
    EllipseParGErrors(TrueParG,ParG3));
fprintf('Error for the Taubin algorithm: Center=%f, Angle=%f, MajorA=%f, MinorA=%f, Area=%f\n',...
    EllipseParGErrors(TrueParG,ParG4));
fprintf('Error for the Halii and Flusser algorithm: Center=%f, Angle=%f, MajorA=%f, MinorA=%f, Area=%f\n',...
    EllipseParGErrors(TrueParG,ParG5));
fprintf('Error for the Rosin algorithm: Center=%f, Angle=%f, MajorA=%f, MinorA=%f, Area=%f\n',...
    EllipseParGErrors(TrueParG,ParG6));

% Compute RMS Orthogonal Error
RSMO = zeros(1,6);

[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG1);
RSMO(1) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for the new algorithm: %f\n',RSMO(1));
[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG2);
RSMO(2) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for Munoz (2014) algorithm: %f\n',RSMO(2));
[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG3);
RSMO(3) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for Fitzgibbon (1999) algorithm: %f\n',RSMO(3));
[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG4);
RSMO(4) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for Taubin (1991) algorithm: %f\n',RSMO(4));
[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG5);
RSMO(5) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for Halii and Flusser (1998) algorithm: %f\n',RSMO(5));
[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG6);
RSMO(6) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for Rosin algorithm: %f\n',RSMO(6));

