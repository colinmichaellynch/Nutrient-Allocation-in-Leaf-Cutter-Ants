clear

%% Set constants and free parameters 

T = 200;
N = 100; 
K = 3; %assumption: ants can only gather food
k = 2; %1.75
D = 2; %2
R = D;%D;

dietPairs = [2 10; 4 10; 6 10; 4 9; 6 9]; 
intakeTargets = [6.9 6.3 7.4 6.5 7.3]; 
%per row: cb, pb, ca, pa
foodAmounts = [.5617 .2809 .7504 .075; .6535 .1634 .7504 .075; .6883 .1152 .7504 .075; .6535 .1634 .7302 .0811; .6883 .1152 .7302 .0811];  

for i = 1:5
    
    cstar = intakeTargets(i); 
    pstar = 1;
    ca = foodAmounts(i, 3); 
    pa = foodAmounts(i, 4);
    cb = foodAmounts(i, 1);  
    pb = foodAmounts(i, 2);  
    Sa = ca/pa;
    Sb = cb/pb;
    sStar = cstar/pstar; 

    pt = zeros(T+1, 1);
    ct = pt; 
    error = pt; 
    error(1) = (sqrt(ct(1)-cstar)^2+(pt(1)-pstar)^2); 
    activeAnts = pt; 
    stateVector = zeros(N, 1)+3; %1 = A, 2 = B, 3 = R
    thresholdMatrix = unifrnd(0, D, K, N); 
    thresholdMatrix(K, :) = 0; 
    rVec = unifrnd(0, R, N, 1);
    %rVec = unifrnd(1, 1.001, N, 1);
    rVecMat(i,:) = rVec; 
    thresholdDifference(i,:) = thresholdMatrix(1,:)-thresholdMatrix(2,:);
    thresholdSum(i,:) = (thresholdMatrix(1,:))+(thresholdMatrix(2,:));

    for t = 2:(T+1)

        %assumption: individuals make decisions simutaneously, not
        %sequentially: not all are in the nest to evaluate stimulus 
        St = ct(t)/pt(t); 
        sNeg = -1*1/((Sa + Sb)/2);
        %sNeg = -1; 
        I = ct(t) - sNeg*pt(t); 
        pStar = I/(sStar-sNeg); 
        cStar = sStar*pStar; 

        if St >= sStar
            d = sqrt((ct(t)-cStar)^2+(pt(t)-pStar)^2); 
        else
            d = -sqrt((ct(t)-cStar)^2+(pt(t)-pStar)^2); 
        end

        for n = 1:N

            if stateVector(n) == 3
                thresholdVector = thresholdMatrix(:, n); 
                zVector = [-d d rVec(n)]'; %assumption: ants have perfect perception, but we can add perceptual error here later
                eVector = exp(k*(zVector-thresholdVector));  
                phiVector = eVector./sum(eVector); 

                cumPhi = cumsum(phiVector); 
                randNum = rand; 
                selectVec = randNum < cumPhi; 
                stateVector(n) = find(selectVec, 1, 'first'); 
            else
                stateVector(n) = 3; 
            end

        end

        nA = sum(stateVector == 1); 
        nB = sum(stateVector == 2); 

        ct(t+1) = ct(t) + ca*nA + cb*nB; 
        pt(t+1) = pt(t) + pa*nA + pb*nB; 

        error(t+1) = sqrt((pt(t+1)-pStar)^2+(ct(t+1)-cStar)^2); 
        stateMatrix(t,:) = stateVector;
        activeAnts(t+1) = nA+nB; 

    end
    
    foragingTrips = stateMatrix == 1 | stateMatrix == 2;
    numberForagingTrips(i, :) = sum(foragingTrips); 
    highProtein = stateMatrix == 2;
    propHighProtein(i, :) = sum(highProtein)./numberForagingTrips(i, :); 
    
end

%Stats

x = [propHighProtein(:) thresholdDifference(:) thresholdSum(:)]; 
x(any(isnan(x), 2), :) = [];
[rho,pval] = partialcorr(x); 

rho = array2table(rho, ...
    'VariableNames',{'propHighProtein','thresholdDifference','thresholdSum'},...
    'RowNames',{'propHighProtein','thresholdDifference','thresholdSum'});

pval = array2table(pval, ...
    'VariableNames',{'propHighProtein','thresholdDifference','thresholdSum'},...
    'RowNames',{'propHighProtein','thresholdDifference','thresholdSum'});

disp('Partial Correlation Coefficients for Prop. High Protein')
disp(rho(:,1))
disp('p-values')
disp(pval(:,1))

x = [numberForagingTrips(:) thresholdDifference(:) thresholdSum(:)]; 
x(any(isnan(x), 2), :) = [];
[rho,pval] = partialcorr(x); 

rho = array2table(rho, ...
    'VariableNames',{'numberForagingTrips','thresholdDifference','thresholdSum'},...
    'RowNames',{'numberForagingTrips','thresholdDifference','thresholdSum'});

pval = array2table(pval, ...
    'VariableNames',{'propHighProtein','thresholdDifference','thresholdSum'},...
    'RowNames',{'propHighProtein','thresholdDifference','thresholdSum'});

disp('Partial Correlation Coefficients for # Foraging Trips')
disp(rho(:,1))
disp('p-values')
disp(pval(:,1))

%scatter(thresholdDifference(:), numberForagingTrips(:))

%% Plot data

%Real data vs simulated data 

realData = readtable('ExperimentalDataSummary.csv');
realDataTotal = readtable('ExperimentalDataTotal.csv');

figure(1)
tiledlayout(3,2)
nexttile
h1 = histogram(realData.prop_high);
h1.NumBins = 20;
xlabel("Prop. Trips to B per Ant", 'FontSize', 15)
ylabel("Count", 'FontSize', 15) 
title("Experimental Data", 'FontSize', 18)
nexttile
h = histogram(propHighProtein(:));
h.NumBins = 20;
xlabel("Prop. Trips to B per Ant", 'FontSize', 15)
ylabel("Count", 'FontSize', 15) 
title("Simulated Data (k = 2, D = 2)", 'FontSize', 18)

nexttile
h1 = histogram(realData.freq_to_all);
h1.NumBins = 20;
xlabel("Number of Foraging Trips per Ant", 'FontSize', 15)
ylabel("Count", 'FontSize', 15) 
nexttile
h = histogram(numberForagingTrips(:));
h.NumBins = 20;
xlabel("Number of Foraging Trips per ant", 'FontSize', 15)
ylabel("Count", 'FontSize', 15) 

J = customcolormap_preset('pasteljet');  

nexttile
scatter(realData.prop_high, realData.freq_to_all,'filled')
xlabel("Prop. Trips to B per Ant", 'FontSize', 15)
ylabel("Number of Foraging Trips", 'FontSize', 15)
nexttile
scatter(propHighProtein(:), numberForagingTrips(:),'filled')
xlabel("Prop. Trips to B per Ant", 'FontSize', 15)
ylabel("Number of Foraging Trips", 'FontSize', 15)

%demonstrative plot for d 

Sa = 1; 
Sb = .5;
sStar = .75; 

proteinVector = 0:.01:100;
AVector = Sa.*proteinVector; 
BVector = Sb.*proteinVector; 
RVector = sStar.*proteinVector;
timeVector = (1:T+1); 

ptTemp = 80;
ctTemp = 70; 

St = ctTemp/ptTemp; 
sNeg = -1*1/((Sa + Sb)/2); 
I = ctTemp - sNeg*ptTemp; 
pStar = I/(sStar-sNeg); 
cStar = sStar*pStar; 

figure(2)
plot(proteinVector, AVector, 'Linewidth', 2, 'color', 'blue')
hold on
plot(proteinVector, BVector, 'Linewidth', 2, 'color', 'red')
%plot(pt, ct, 'Linewidth', 2, 'color', 'black')
plot(proteinVector, RVector, 'Linewidth', 2, 'color', 'yellow')
plot([ptTemp pStar], [ctTemp cStar], ':', 'Linewidth', 2, 'color', 'm')
scatter(ptTemp, ctTemp, 50, 'filled', 'black')
scatter(pStar, cStar, 50, 'filled', 'c')
xlabel("Protein (grams)", 'FontSize', 15)
ylabel("Carbs (grams)", 'FontSize', 15)
legend("A Rail", "B Rail", "Intake Target Rail", "d", "Current Status", "Closest Target", 'Location','northwest', 'FontSize', 12)
axis equal

figure(3)
scatter(propHighProtein(:), numberForagingTrips(:), 20*thresholdSum(:), thresholdDifference(:),'filled')
xlabel("Prop. Trips to B per Ant", 'FontSize', 15)
ylabel("Number of Foraging Trips", 'FontSize', 15)
colormap(J);  
h = colorbar;
ylabel(h, '\theta^A_{nm} - \theta^B_{nm}', 'FontSize', 20)
grid on


%% Comparison of trajectories

colonyList = string(unique(realDataTotal.colony));
dates = string(unique(realDataTotal.date));
dietPairings = (unique(realDataTotal.diet_pairing));

%best = col3, dates 4, diet2
Experiment1 = realDataTotal(realDataTotal.colony==colonyList(3) & realDataTotal.date==dates(4) & realDataTotal.diet_pairing==dietPairings(2), :);
Experiment1.protein = string(Experiment1.protein); 

dietPairs = [2 10; 4 10; 6 10; 4 9; 6 9]; 
intakeTargets = [6.9 6.3 7.4 6.5 7.3]; 
%per row: cb, pb, ca, pa
foodAmounts = [.5617 .2809 .7504 .075; .6535 .1634 .7504 .075; .6883 .1152 .7504 .075; .6535 .1634 .7302 .0811; .6883 .1152 .7302 .0811];
i = 5;
cstar = intakeTargets(i); 
pstar = 1;
ca = foodAmounts(i, 3); 
pa = foodAmounts(i, 4);
cb = foodAmounts(i, 1);  
pb = foodAmounts(i, 2);  
Sa = ca/pa;
Sb = cb/pb;
sStar = cstar/pstar; 

ctReal = 0;
ptReal = 0; 

for t = 2:(height(Experiment1))
    
    if Experiment1.protein(t) == "High"
        
        ctReal(t) = ctReal(t-1) + cb;
        ptReal(t) = ptReal(t-1) + pb;
        
    else
        
        ctReal(t) = ctReal(t-1) + ca;
        ptReal(t) = ptReal(t-1) + pa;
        
    end
    
end

proteinVector = 0:.005:max(ptReal);
AVector = Sa.*proteinVector; 
BVector = Sb.*proteinVector; 
RVector = sStar.*proteinVector;

pt = zeros(T+1, 1);
ct = pt; 
stateVector = zeros(N, 1)+3; %1 = A, 2 = B, 3 = R
thresholdMatrix = unifrnd(0, D, K, N); 
thresholdMatrix(K, :) = 0; 
rVec = unifrnd(0, R, N, 1);

for t = 2:(T+1)

    St = ct(t)/pt(t); 
    sNeg = -1*1/((Sa + Sb)/2); 
    I = ct(t) - sNeg*pt(t); 
    pStar = I/(sStar-sNeg); 
    cStar = sStar*pStar; 

    if St >= sStar
        d = sqrt((ct(t)-cStar)^2+(pt(t)-pStar)^2); 
    else
        d = -sqrt((ct(t)-cStar)^2+(pt(t)-pStar)^2); 
    end

    for n = 1:N

        if stateVector(n) == 3
            thresholdVector = thresholdMatrix(:, n); 
            zVector = [-d d rVec(n)]';
            eVector = exp(k*(zVector-thresholdVector));  
            phiVector = eVector./sum(eVector); 

            cumPhi = cumsum(phiVector); 
            randNum = rand; 
            selectVec = randNum < cumPhi; 
            stateVector(n) = find(selectVec, 1, 'first'); 
        else
            stateVector(n) = 3; 
        end

    end

    nA = sum(stateVector == 1); 
    nB = sum(stateVector == 2); 

    ct(t+1) = ct(t) + ca*nA + cb*nB; 
    pt(t+1) = pt(t) + pa*nA + pb*nB; 

end

figure(4) 
plot(proteinVector, AVector, 'Linewidth', 2, 'color', 'blue')
hold on
plot(proteinVector, BVector, 'Linewidth', 2, 'color', 'red')
plot(proteinVector, RVector, 'Linewidth', 2, 'color', 'yellow')
plot(ptReal, ctReal, 'Linewidth', 2, 'color', 'black')
plot(pt, ct, 'Linewidth', 2, 'color', 'm')
xlabel("Protein (grams)", 'FontSize', 15)
ylabel("Carbs (grams)", 'FontSize', 15)
legend("A Rail", "B Rail", "Intake Target Rail", "Real Colony Trajectory", "Simulated Colony Trajectory", 'Location','southeast', 'FontSize', 12)
xlim([0, max(ptReal)])
ylim([0, max(ctReal)])

%% Effect of k, D 

kVec = [.5 4];
dVec = [.5 4]; 

figure(5)
tiledlayout(2,2)

for j = 1:2
    for z = 1:2
        k = kVec(j);
        D = dVec(z); 
        for i = 1:5
            
            cstar = intakeTargets(i); 
            pstar = 1;
            ca = foodAmounts(i, 3); 
            pa = foodAmounts(i, 4);
            cb = foodAmounts(i, 1);  
            pb = foodAmounts(i, 2);  
            Sa = ca/pa;
            Sb = cb/pb;
            sStar = cstar/pstar; 

            pt = zeros(T+1, 1);
            ct = pt; 
            error = pt; 
            error(1) = (sqrt(ct(1)-cstar)^2+(pt(1)-pstar)^2); 
            activeAnts = pt; 
            stateVector = zeros(N, 1)+3; %1 = A, 2 = B, 3 = R
            thresholdMatrix = unifrnd(0, D, K, N); 
            thresholdMatrix(K, :) = 0; 
            rVec = unifrnd(0, R, N, 1);
            %rVec = unifrnd(1, 1.001, N, 1);
            rVecMat(i,:) = rVec; 
            thresholdDifference(i,:) = thresholdMatrix(1,:)-thresholdMatrix(2,:);
            thresholdSum(i,:) = (thresholdMatrix(1,:))+(thresholdMatrix(2,:));

            for t = 2:(T+1)

                %assumption: individuals make decisions simutaneously, not
                %sequentially: not all are in the nest to evaluate stimulus 
                St = ct(t)/pt(t); 
                sNeg = -1*1/((Sa + Sb)/2); 
                I = ct(t) - sNeg*pt(t); 
                pStar = I/(sStar-sNeg); 
                cStar = sStar*pStar; 

                if St >= sStar
                    d = sqrt((ct(t)-cStar)^2+(pt(t)-pStar)^2); 
                else
                    d = -sqrt((ct(t)-cStar)^2+(pt(t)-pStar)^2); 
                end

                for n = 1:N

                    if stateVector(n) == 3
                        thresholdVector = thresholdMatrix(:, n); 
                        zVector = [-d d rVec(n)]'; %assumption: ants have perfect perception, but we can add perceptual error here later
                        eVector = exp(k*(zVector-thresholdVector));  
                        phiVector = eVector./sum(eVector); 

                        cumPhi = cumsum(phiVector); 
                        randNum = rand; 
                        selectVec = randNum < cumPhi; 
                        stateVector(n) = find(selectVec, 1, 'first'); 
                    else
                        stateVector(n) = 3; 
                    end

                end

                nA = sum(stateVector == 1); 
                nB = sum(stateVector == 2); 

                ct(t+1) = ct(t) + ca*nA + cb*nB; 
                pt(t+1) = pt(t) + pa*nA + pb*nB; 

                error(t+1) = sqrt((pt(t+1)-pStar)^2+(ct(t+1)-cStar)^2); 
                stateMatrix(t,:) = stateVector;
                activeAnts(t+1) = nA+nB; 

            end

            foragingTrips = stateMatrix == 1 | stateMatrix == 2;
            numberForagingTrips(i, :) = sum(foragingTrips); 
            highProtein = stateMatrix == 2;
            propHighProtein(i, :) = sum(highProtein)./numberForagingTrips(i, :); 

        end
        
        str = "D = " + D + ", k = " + k;
        
        nexttile
        scatter(propHighProtein(:), numberForagingTrips(:),'filled')
        xlabel("Prop. Trips to B per Ant", 'FontSize', 15)
        ylabel("Number of Foraging Trips", 'FontSize', 15)
        title(str, 'FontSize', 15)
        
    end
    
end