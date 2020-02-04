function [ Results ] = roc( GT, Resp )
%=================================================
%  Receiver Operating Characteristics (ROC) curve  
%=================================================
% This function performs ROC curve, computing sensitivity and specificity.
% It additionally computes false positive rate (FPR), false negative rate
% (FNR), positive predictive value (PPV), the negative predictive value
% (NPV), the optimal point with same penalty for false positive and false
% negative (and plot d'), and the area under the ROC curve (AUC).
% All the results are stored in the output structure.
%
%   USAGE: roc( GT, Resp)
%          GT being a 1 x n or n x 1 vector corresponding to the ground
%          truth and Resp being a 1 x n or n x 1 vector corresponding to
%          the responses of the receiver. The data in the two input vectors
%          should be somehow binary.
%
%---------------------------------------------------------------------
% Renaud Marquis, with help from Sandrine Muller, July 2013
%---------------------------------------------------------------------

ModalitiesGT = unique(GT);
ModalitiesResp = unique(Resp);

if length(ModalitiesGT) > 2 || length(ModalitiesResp) > 2
    fprintf('\nError! Ground truth or Responses are not binary vectors!\n\n')
    return
end

for i = 1:length(GT)
    if GT(i) == ModalitiesGT(1) && Resp(i) == ModalitiesResp(1)
        ContingencyTable(i,1) = 1; % TP (True positive) ("hit")
    end
    if GT(i) == ModalitiesGT(1) && Resp(i) == ModalitiesResp(2)
        ContingencyTable(i,2) = 1; % FN (False negative) ("miss")
    end
    if GT(i) == ModalitiesGT(2) && Resp(i) == ModalitiesResp(1)
        ContingencyTable(i,3) = 1; % FP (False positive) ("false alarm")
    end 
    if GT(i) == ModalitiesGT(2) && Resp(i) == ModalitiesResp(2)
        ContingencyTable(i,4) = 1; % TN (True negative) ("correct reject")
    end
end

CumulativeContingencyTable = cumsum(ContingencyTable);

Sensitivity = CumulativeContingencyTable(:,1)./ (CumulativeContingencyTable(:,1)+CumulativeContingencyTable(:,2));
Specificity = CumulativeContingencyTable(:,4)./ (CumulativeContingencyTable(:,3)+CumulativeContingencyTable(:,4));
Sensitivity(isnan(Sensitivity)) = 0;
Specificity(isnan(Specificity)) = 0;
% Sensitivity(Sensitivity == 1) = 0;
% Specificity(Specificity == 1) = 0;

FalsePositiveRate = 1-Specificity; % FPR
FalseNegativeRate = 1-Sensitivity; % FNR
FalsePositiveRate(isnan(FalsePositiveRate)) = 0;
FalseNegativeRate(isnan(FalseNegativeRate)) = 0;

PositivePredictiveValue = CumulativeContingencyTable(:,1)./ (CumulativeContingencyTable(:,1)+CumulativeContingencyTable(:,3)); % PPV
NegativePredictiveValue = CumulativeContingencyTable(:,4)./ (CumulativeContingencyTable(:,2)+CumulativeContingencyTable(:,4)); % NPV
PositivePredictiveValue(isnan(Sensitivity)) = 0;
NegativePredictiveValue(isnan(NegativePredictiveValue)) = 0;

LinearInterpolation = interp1q(Sensitivity,Specificity,Sensitivity);

indicesInf = find(isinf(LinearInterpolation),length(LinearInterpolation));

LinearInterpolation(indicesInf) = nan(1,1);

PerpendicularRandomModel = 1-(0:(1/length(LinearInterpolation)):1);

Intersection = abs(LinearInterpolation-(PerpendicularRandomModel(1:end-1))');

index = find(Intersection == min(Intersection));

dPrime = [1-(LinearInterpolation(index)) LinearInterpolation(index)];

AUC = ignoreNaN(LinearInterpolation,@trapz)/length(LinearInterpolation); % Area under the ROC curve (AUC)

figure('Position',[300 300 600 600]);
scatter(Specificity,Sensitivity,'.b');
hold on;
Xscale = 0:(1/length(LinearInterpolation)):1;
plot(Xscale(1:end-1),LinearInterpolation,'m')
axis([0 1 0 1]);
xlabel('Specificity');
ylabel('Sensitivity');
RandomModel = line([1 0],[0 1]);
set(RandomModel,'Color','r');
XequalY = line([0 1],[0 1]);
set(XequalY,'Color','g');
hold on;
dPrime1 = dPrime(1);
dPrime2 = dPrime(2);
plot(dPrime1,dPrime2,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)

Results.ContingencyTable = ContingencyTable;
Results.CumulativeContingencyTable = CumulativeContingencyTable;
Results.Sensitivity = Sensitivity;
Results.Specificity = Specificity;
Results.FPR = FalsePositiveRate;
Results.FNR = FalseNegativeRate;
Results.PPV = PositivePredictiveValue;
Results.NPV = NegativePredictiveValue;
Results.LinearInterpolation = LinearInterpolation;
Results.AUC = AUC;
Results.dPrime = dPrime;

end

