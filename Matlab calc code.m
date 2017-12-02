%objective - pull in excel file, calculate optimal portfolio weights,
%single factor model weights and famma french for comparison.  Calculate
%performance too.

%importing excel file, 1-150

data=[];

optimalSharp=0
singleSharp=0   
    
    
for fileNumber=1:1:250
    %reads in excel file;
    filename=int2str(fileNumber)+"_data_practice.xls";
    historical = xlsread(filename,'returns','','basic');
    characteristics = xlsread(filename,'characteristics','','basic');
    future = xlsread(filename,'future_returns','','basic');

    equalWeights=[.1,.1,.1,.1,.1,.1,.1,.1,.1,.1];
    riskFree=future(end,end);

    %caculates optimatal future portolio weights
    futureMean=mean(future);
    futureMean=futureMean(:,2:end-1);

    %calc std dev
    futurestd=std(future,'',1);
    futurestd=futurestd(:,2:end-1);

    %calc correlation matrix
    futureCorrM=corrcoef(future);
    futureCorrM=futureCorrM(2:end-1,2:end-1);

    %calc std matrix
    futureStdM=futurestd.*futurestd';

    %calc final scalar matrix that will be multiplied against the wieghts
    %matrix to create variance matrix.
    scalarM=futureStdM.*futureCorrM;
    
    %define function which takes in a vector of size 10, and returns the sharp
    %ratio given the scalarM above

    sharpfunction = @(weightsVector, Scalar) (sum(weightsVector.*futureMean)-riskFree)/(sqrt(sum(sum((weightsVector.*weightsVector').*scalarM))));

    sharpfunction(equalWeights, scalarM);

    %optimize function f(x) that takes in vector of size 10 with bounds of -2
    %to +2 and which sum to 1.

    %need to create negative as matlab can only minimize
    negativeSharpFunction = @(weightsVector, Scalar) -(sum(weightsVector.*futureMean)-riskFree)/(sqrt(sum(sum((weightsVector.*weightsVector').*scalarM))));
    negativeSharpFunction(equalWeights, scalarM);

    %defining characteristics
    con = @sumWeightsEqualOne; %need to define as =0
    lowerBound=[-2,-2,-2,-2,-2,-2,-2,-2,-2,-2];
    upperBound=-lowerBound;


    %final optimization function
    idealWeights = fmincon(negativeSharpFunction,equalWeights,[],[],[],[],lowerBound,upperBound,con);
    optimalSharp=optimalSharp+sharpfunction(idealWeights, scalarM);

    data(fileNumber,:)=idealWeights;
    
        %CALCULATE SINGLE FACTOR MODEL
    
    historicalReturnSF=historical(:,5:14);
    for j=5:1:14
        historicalReturnSF(:,j-4)=historicalReturnSF(:,j-4)-historical(:,4);
    end
    MeanSF=mean(historicalReturnSF);
    stdSF=std(historicalReturnSF,'',1);
    
    MarketReturnSF=historical(:,1);
    MarketReturnAverage=mean(MarketReturnSF);
    MarketReturnSTD=std(MarketReturnSF);
    
    riskFreeSF=historical(60,4);
    
    factors=zeros(0,10);
    X = [ones(length(MarketReturnSF),1) MarketReturnSF];
    for i=1:1:10
            b = X\historicalReturnSF(:,i);
            factors=[factors,b];
    end
    
    predicted=X*factors;
    epsilon=historicalReturnSF-predicted;
    epsilonMean=mean(epsilon);
    epsilonSTD=std(epsilon,'',1);
    
    %future expected    
    X = [ones(length(MarketReturnAverage),1) MarketReturnAverage];
    returnNoRF=X*factors;
    expectedReturn=returnNoRF+riskFreeSF;
    
    alpha=factors(1,:);
    beta=factors(2,:);
    
    
    SFFinalVar=(beta.^2)*MarketReturnSTD^2+epsilonSTD.^2;
    SFFinalSTD=SFFinalVar.^.5;
    
    %matrices
    varianceMatrix=beta.*beta'*MarketReturnSTD^2;
    for ii=1:1:10
       varianceMatrix(ii,ii)=SFFinalVar(ii);
    end
    stdMatrix=SFFinalSTD.*SFFinalSTD';  
    corrMatrix=varianceMatrix./stdMatrix;
    
    %final factors for optimization single index
    %corrMatrix, expectedReturn, SFFinalSTD

    %calc final scalar matrix that will be multiplied against the wieghts
    %matrix to create variance matrix.
    scalarM=corrMatrix.*stdMatrix;
    
    %define function which takes in a vector of size 10, and returns the sharp
    %ratio given the scalarM above

    sharpfunction2 = @(weightsVector, Scalar) (sum(weightsVector.*expectedReturn)-riskFreeSF)/(sqrt(sum(sum((weightsVector.*weightsVector').*scalarM))));

    sharpfunction2(equalWeights, scalarM);

    %optimize function f(x) that takes in vector of size 10 with bounds of -2
    %to +2 and which sum to 1.

    %need to create negative as matlab can only minimize
    negativeSharpFunction2 = @(weightsVector, Scalar) -(sum(weightsVector.*expectedReturn)-riskFreeSF)/(sqrt(sum(sum((weightsVector.*weightsVector').*scalarM))));
    negativeSharpFunction2(equalWeights, scalarM);

    %defining characteristics
    con = @sumWeightsEqualOne; %need to define as =0
    lowerBound=[-2,-2,-2,-2,-2,-2,-2,-2,-2,-2];
    upperBound=-lowerBound;

    %final optimization function
    idealWeightsSingle = fmincon(negativeSharpFunction,equalWeights,[],[],[],[],lowerBound,upperBound,con);
    
    singleSharp=singleSharp+sharpfunction2(idealWeightsSingle, scalarM);

    
    dataSingleFactor(fileNumber,:)=idealWeightsSingle;
end

singleSharp/250
optimalSharp/250
    