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

% Path of the images
ImagePath = './images';
%% DataEllipse2.1
Images = {'169_0015','195_0039','rueda6','rueda11'};
Suffix = {'','','',''};


%% Batch running all images
NumImages=numel(Images);
NumMethods = 6;
Labels = {'SAREfit','Muñoz','Fitzgibbon','Taubin','Halir&Flusser','Rosin'};

for NdxImage=1:NumImages
    
    ThisImage = Images{NdxImage};
    ThisSuffix = Suffix{NdxImage};
    ThisData = [ThisImage ThisSuffix];

    % Load MAT file with detected points (done with the desired algorithm)
    % It should have the same name of the image and suffix 'p'
    MatImage = strcat(sprintf('%s/p',ImagePath),ThisData);
    load(MatImage)
    X = [x;y];

    % Execute the methods
        
    % SAREfit
    SubsamplingFactor=0.15; % Proportion of the total training samples considered for training
    EnsembleSize=90;        % Number of subsamplings
    percentile=10;          % Percentage(%) of best fits for the final ensemble
    [ParA1,ParG1,ParN1]=SAREfit(X,SubsamplingFactor,EnsembleSize,percentile);

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

    
    % Plot results
    MyColors=distinguishable_colors(NumMethods+1);
    
    % Make grey the black one
    MyColors(4,:)=MyColors(4,:)+0.5;
    
    % Make yellow as True Ellipse
    aux = MyColors(6,:);
    MyColors(6,:) = MyColors(1,:);
    MyColors(1,:) = aux;
    
    % Switch color orders for the sake of clarity
    aux = MyColors(6,:);
    MyColors(6,:) = MyColors(5,:);
    MyColors(5,:) = aux;

    % Plot image with points and ellipses
    figure
    JpgImage = imread(strcat(sprintf('%s/',ImagePath),ThisImage,'.jpg'));
    imshow(JpgImage)
    PositionO = get(gca,'OuterPosition');
    PositionI = get(gca,'InnerPosition');
    alpha 0.7
    hold on
    Handles=zeros(1,NumMethods+1);
    if exist('TrueParG','var')
        MyHandle=PlotEllipseG(TrueParG,MyColors(1,:),2);
        Handles(1)=MyHandle(1);
    end
    MyHandle=PlotEllipseG(ParG2,MyColors(3,:),0.75);
    Handles(3)=MyHandle(1);
    MyHandle=PlotEllipseG(ParG3,MyColors(4,:),0.75);
    Handles(4)=MyHandle(1);
    MyHandle=PlotEllipseG(ParG4,MyColors(5,:),0.75);
    Handles(5)=MyHandle(1);
    MyHandle=PlotEllipseG(ParG5,MyColors(6,:),0.75);
    Handles(6)=MyHandle(1);
    MyHandle=PlotEllipseG(ParG6,MyColors(7,:),0.75);
    Handles(7)=MyHandle(1);
    % SAREfit is drawn at last
    MyHandle=PlotEllipseG(ParG1,MyColors(2,:),0.75);
    Handles(2)=MyHandle(1);

    % Draw samples
    plot(X(1,:), X(2, :), '+r', 'Color', MyColors(1,:),'MarkerSize',2);
    axis equal tight
    axis([0 size(JpgImage,2) 0 size(JpgImage,1)])
    set(gca,'YDir','reverse')
    set(gca,'OuterPosition',PositionO)
    set(gca,'InnerPosition',PositionI)

    PdfFileName=sprintf('%s/%s_ImageEllipses',ImagePath,ThisData);
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperOrientation','portrait');
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperSize',[8 7]);
    set(gcf,'PaperPosition',[0 0 8 7]);
    set(gca,'fontsize',10);
    saveas(gcf,PdfFileName,'pdf');
   
%     legend(Handles,'SAREfit','Muñoz','Fitzgibbon','Taubin','Halir&Flusser','Rosin','Location','southoutside','Orientation','horizontal');
    
end




