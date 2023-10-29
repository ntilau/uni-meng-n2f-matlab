%%% applies the empirical interpolation to the following parameter mu 
%%% dependent function (1-x).* cos(3*pi*mu(i)*(x+1)).*exp(-(1+x).*mu(i))
%%% see Patera. testRadOperator=true enable approximation of the complex
%%% valued operator on Jsx to EtFF --> not well approximated
clc; clear all; close all;
testRadOperator=false;
if testRadOperator
    load nf2ffOp.mat; % operator and x position
end
%%
if testRadOperator
    mu=linspace(0,2*pi,size(s,1));
    x=linspace(-1,1,size(s,2));
else
    x = -1:0.01:1;
    mu = 1:(pi-1)/2000:pi;
    s=zeros(length(mu),length(x));
end
maxNbrVectors = 50;
L2=zeros(1,maxNbrVectors);
% ChosenNbrOfVectors = 10;
for ChosenNbrOfVectors=5:5:maxNbrVectors
    if ~testRadOperator
        for i=1:length(mu)
           s(i,:) = (1-x).* cos(3*pi*mu(i)*(x+1)).*exp(-(1+x).*mu(i));
        end
    end
    Size = size(s);
    % figure
    % surf(s,'EdgeColor','none')
    % figure
    % semilogy(svd(s,0),'-*');
    Indices = floor(1:length(mu)/ChosenNbrOfVectors:length(mu));
    TmpBasisA = s(Indices,:);

    [U,S,V] = svd(TmpBasisA,0);
    POD = U'*TmpBasisA;

    NbrOfVectors = length(Indices);
    TotalNbrOfVectors = length(mu);

    RowIdxA = zeros(1,NbrOfVectors);
    ColIdxA = RowIdxA;
    WeightingA = RowIdxA;

    [MaxOfCols, RowIdxOfColMax] = max(abs(TmpBasisA)); % max of each column
    [AbsoluteMax, ColIdxA(1)] = max(MaxOfCols); % column with higher val
    RowIdxA(1) = RowIdxOfColMax(ColIdxA(1)); % t= ColIdx

    for k=2:NbrOfVectors
        ValA = zeros(1, NbrOfVectors);
        FirstValA = 1;
        for i=1:TotalNbrOfVectors
            for j=1:NbrOfVectors
                if(~any(RowIdxA == j) && i~=Indices(j))
                    ValTmpA = sqrt(sum(abs(s(i,:) - TmpBasisA(j,:)).^2));
                    if(FirstValA)
                        ValA(j)= ValTmpA;
                        FirstValA = 0;
                    else
                        if(ValTmpA < ValA(j)) % Infimum
                            ValA(j)= ValTmpA;
                        end
                    end
                end
            end    
        end
        [MaxResidueVal, RowIdxA(k)] = max(ValA); % sets the next vector
    end

    OrderedIndicesA(1:NbrOfVectors) = Indices(RowIdxA(1:NbrOfVectors));

    TmpBasisA((1:NbrOfVectors),:) = TmpBasisA(RowIdxA(1:NbrOfVectors),:);
    WeightingA(1) = TmpBasisA(1,ColIdxA(1));
    BasisA = zeros(NbrOfVectors,Size(2));
    BasisA(1,:)= TmpBasisA(1,:) / WeightingA(1); % q = basis
    ResidueA = zeros(NbrOfVectors,Size(2));
    ResidueA(1,:) = TmpBasisA(1,:);

    for i=2:NbrOfVectors
    %     Sigma = zeros(i-1,1);
        TmpMatrixA = zeros(i-1,i-1);
        XsiA = zeros(i-1,1);
        for j=1:(i-1)
            for k=1:(i-1)
                TmpMatrixA(j,k) = BasisA(k,ColIdxA(j));
            end
            XsiA(j,1) = TmpBasisA(i,ColIdxA(j));
        end
        SigmaA = inv(TmpMatrixA) * XsiA;
        ApproxXsiA = zeros(size(TmpBasisA(i,:)));
        for j=1:i-1
            ApproxXsiA = ApproxXsiA + SigmaA(j)*BasisA(j,:);
        end
        ResidueA(i,:) = TmpBasisA(i,:) - ApproxXsiA;
        [AbsoluteMax, ColIdxA(i)] = max(abs(ResidueA(i,:)));
        WeightingA(i) = ResidueA(i,ColIdxA(i));
        BasisA(i,:) = ResidueA(i,:) / WeightingA(i);
    end
    % figure;
    % semilogy(svd(TmpBasisA,0),'-*');
    %% Analysis
    % for i=1:3
    %     figure
    %     plot(x,BasisA(i,:),'-*b',x,ResidueA(i,:),':.r',x,POD(i,:),':.g');
    % end
    %%
    Color = ['b','r','g','c','m','y','k'];
    % figure;
    % for i=1:6
    %     plot(x,BasisA(i,:),Color(i),'LineWidth',2);
    %     hold on;
    %     plot(x(1,ColIdxA(i)),0,['o',Color(i)],'LineWidth',2);
    % end
    %%
    QMatrix = zeros(NbrOfVectors,NbrOfVectors);
    for i=1:NbrOfVectors
        for j=1:NbrOfVectors
            QMatrix(i,j) = BasisA(j,ColIdxA(i));
        end
    end
    %% Spline
    sValues = TmpBasisA(:,ColIdxA(:));
    sValuesSplined = zeros(Size(1),NbrOfVectors);
    for i=1:NbrOfVectors
        sValuesSplined(:,i) = spline(mu(Indices), sValues(:,i), mu);
    end
    Beta = (inv(QMatrix) * sValuesSplined.');
    % figure;
    % plot(mu(Indices), sValues(:,100), mu, sValuesSplined(:,100));
    % surf(sValuesSplined,'EdgeColor','none');
    % surf(Beta.','EdgeColor','none');
    %% sApprox
    sApprox = Beta.'*BasisA;
    % figure;
    % for i=1:6
    %     plot(x,sApprox(i,:),Color(i),'LineWidth',2);
    %     hold on;
    % end
    %% 
    figure(1)
    subplot(1,2,1)
    if testRadOperator
        surf(x,mu,abs(s),'EdgeColor','none'); view(25,25);
    else
        surf(x,mu,s,'EdgeColor','none'); view(25,25)
    end
    title('Original');
    axis equal;
    subplot(1,2,2)
    if testRadOperator
        surf(x,mu,abs(sApprox),'EdgeColor','none');view(25,25);
    else
        surf(x,mu,(sApprox),'EdgeColor','none');view(25,25);
    end
    title('Approximated');
    axis equal;
    %% L2
    L2Error = sqrt(sum(sum( (abs(s-sApprox)).^2)) / sum(sum(abs(s).^2)));
    figure(2)
    if testRadOperator
         surf(x,mu,abs(s-sApprox),'EdgeColor','none');  view(25,25)
    else
         surf(x,mu,(s-sApprox),'EdgeColor','none');  view(25,25)
    end
   
    title(['L^2 Error: ',num2str(L2Error),', NbrVects: ', num2str(length(Indices)),' over ',...
        num2str(Size(1)),' x ',num2str(Size(2))]);
    axis equal;
    L2(ChosenNbrOfVectors) = L2Error;
%     close all
end
%%
svd(BasisA)
Vects = 1:ChosenNbrOfVectors;
figure(3)
semilogy(Vects, L2,'o',Vects,(svd(BasisA)).','*');
xlabel('Nbr Of Vectors');
legend('L^2_{RelError}', '\sigma_{SVD}')
title('EIM Performance')
