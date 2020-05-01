function [ParA,ParG,ParN]=SAREfit(X,SubsamplingFactor,EnsembleSize,percentile)
% THIS VERSION USES THE PREDICTED ELLIPSE AS GT TO REMOVE THE WORST FITTED
% ELLIPSES AND RE-COMPUTE THE MEDIAN ELLIPSE
% Fit an ellipse to a set of training points.
% Robust random ensemble of basic fits.
% Inputs:
%   X=matrix of size 2 x NumPoints with the training samples
%   SubsamplingFactor=Probability of being chosen for a training subset,
%       real scalar from 0 to 1.
%   EnsembleSize=Number of basic fits in the random ensemble, strictly positive
%       integer
%   percentile=Percentile of best fits considered in the ensemble of fitted
%   ellipses
% Outputs:
%   ParG = [Xcenter, Ycenter, a, b, AngleOfTilt]' is the vector of 
%   geometric parameters of the ellipse). a=half length of the major axis,
%   b=half length of the minor axis
%   ParA = (A,B,C,D,E,F)' is the vector of Algebraic parameters:
%           Ax^2 + 2Bxy + Cy^2 + 2Dx + 2Ey + F = 0
%   ParN= [Focus1x Focus1y Focus2x Focus2y SumDists] is the vector of
%   natural parameters of the ellipse: the two foci (Focus1 and Focus2),
%   and the sum of distances to both foci, SumDists==2a

% VERSION WITH FIXED SEED FOR PARAMETER FITTING
s = RandStream('mt19937ar','Seed',0);
% s = parallel.pool.Constant(RandStream('Threefry','Seed',0));

% Get the number of training samples
[~,NumSamples]=size(X);
NumSubsamples = floor(SubsamplingFactor*NumSamples);

% Precompute PARE algorithm initialization by Fitzgibbon (1999)
A = EllipseDirectFit(X');
ParAIni=[A(1) 0.5*A(2) A(3) 0.5*A(4) 0.5*A(5) A(6)]';
ParGIni=real(AtoG(ParAIni));

AllParN=zeros(5,EnsembleSize);
NdxNumEnsemblesUsed = 1;
% Compute the members of the ensemble
for NdxEnsemble=1:EnsembleSize
    
    % Choose a training subset at random
    RandomOrder = randperm(s,NumSamples);
    TrainingSubset=X(:,RandomOrder(1:NumSubsamples));

    % Halir&Flusser method           
    A = EllipseFitHalir(TrainingSubset(1,:), TrainingSubset(2,:));
    if ~isempty(A) && isreal(A)
        % Obtain the natural parameters of the ellipse
        ParA=[A(1) 0.5*A(2) A(3) 0.5*A(4) 0.5*A(5) A(6)]';
        [ParG,code]=AtoG(ParA);
        if (code==1)    % no error
            AllParN(:,NdxNumEnsemblesUsed)=GtoN(ParG);
            NdxNumEnsemblesUsed = NdxNumEnsemblesUsed+1;
        end
    end

end
NdxNumEnsemblesUsed=NdxNumEnsemblesUsed-1;  % Adjust number of fits

% Compute the consensus from the natural parameters of the ensemble members
[~,~,ParN]=computeConsensusEllipse(AllParN,ParGIni);

% Compute natural errors with respect to the median ellipse
ParNErrors = zeros(1,NdxNumEnsemblesUsed);
for NdxEllipse = 1:NdxNumEnsemblesUsed
    
    ParNErrors(NdxEllipse) = EllipseNaturalError(ParN,AllParN(:,NdxEllipse));
    
end

% Compute threshold of the natural errors
thresholdError = prctile(ParNErrors,percentile);

% Choose best ellipses
bestParN = AllParN(:,ParNErrors<=thresholdError);

% Compute the consensus from the best ensemble members
[ParA,ParG,ParN]=computeConsensusEllipse(bestParN,ParGIni);

end



function [ParA,ParG,ParN]=computeConsensusEllipse(AllParN,ParGIni)

% Compute the consensus from the natural parameters of the ensemble members
ParN=zeros(5,1);
ParN(1:2)=L1mediancovMEX(AllParN(1:2,:));
ParN(3:4)=L1mediancovMEX(AllParN(3:4,:));
ParN(5)=median(AllParN(5,:));
% Compute geometric parameters
Focus1=ParN(1:2);
Focus2=ParN(3:4);
Center=(Focus1+Focus2)/2;  % Ellipse center
df=norm(Focus1-Focus2); % Focal distance
m=ParN(5);

% Check whether this is an ellipse
if m>df
    % It is an ellipse, so it is OK
    a=m/2; % Length of the half major axis
    b=0.5*sqrt(m^2-df^2); % Length of the half minor axis
    phi=atan((Focus2(2)-Focus1(2))/(Focus2(1)-Focus1(1))); % Tilt angle
    ParG=[Center(1) Center(2) a b phi]';
else
    % It is not an ellipse, so we revert to the Fitzgibbon solution
    ParG=ParGIni;
    disp('ERROR: The consensus is not an ellipse. Using Fitzgibbon initialization')
end


% Convert to algebraic form
ParA=GtoA(ParG,1);

% Convert to natural form
ParN=GtoN(ParG);

end