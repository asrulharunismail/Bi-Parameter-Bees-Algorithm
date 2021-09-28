clc;
clear;
close all;
%% Problem Definition
[typeOfFunction] = 'EilA101';
Instance=TsplibVRP(typeOfFunction);
Dims=Instance.dim;
ObjFunction=@(x) Instance.evaluation( x );             
VarSize=[1 Dims];                                      
%% Bees Algorithm Parameters
n= 10; nep = 40; 
MaxEval = 1000000;
recruitment = round(linspace(nep,1,n));
assigntment = round(linspace(1,Dims,n)); 
ColonySize=sum(recruitment)        
MaxIt=round(MaxEval/ColonySize);

%% Initialization
Empty_Patch.Position=[]; Empty_Patch.Cost=[]; Empty_Patch.Sol=[];
Empty_Patch.Counter=[];
Patch=repmat(Empty_Patch,n,1);
counter=0;
for i=1:n
    [Patch(i).Position,Patch(i).Cost, Patch(i).Sol]=BiBA_Clustering_using_Coverage(Instance);
	counter = counter + 5000;
    Patch(i).Counter = counter;
end
% Sort
[~, SortOrder]=sort([Patch.Cost]);
Patch=Patch(SortOrder);
BestSol=Patch(1);
% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);
Counter=zeros(MaxIt,1);
OptSol.Cost=inf;
%% Bees Algorithm Main Loop
P=1;
for it=1:MaxIt
    if counter >= MaxEval
        break;
    end
    % Patch Sites
    for i=1:n
        BestForager.Cost=inf;
        for j=1:recruitment(i)
            ForagerBees.Position= Foraging_VRP(Patch(i).Position,assigntment(i));
            [ForagerBees.Cost,ForagerBees.Sol]=ObjFunction(ForagerBees.Position);
            ForagerBees.Position= BRO_0_VRP(ForagerBees,Instance,1,1);
            [ForagerBees.Cost,ForagerBees.Sol]=ObjFunction(ForagerBees.Position);
            counter = counter + 1;
            ForagerBees.Counter = counter;
            if ForagerBees.Cost<BestForager.Cost && sum (ForagerBees.Sol.UC)== sum (Instance.r) %%%
                BestForager=ForagerBees;
            end
            
        end
        if BestForager.Cost<Patch(i).Cost
            Patch(i)=BestForager;
        end
        
    end    
    % Sort
    [~, SortOrder]=sort([Patch.Cost]);
    Patch=Patch(SortOrder);
    BestSol=Patch(1);
    if BestSol.Cost < OptSol.Cost
        OptSol=BestSol;
    end
    % Store Best Cost Ever Found
    BestCost(it)=OptSol.Cost;
    Counter(it)=OptSol.Counter;
    OPTSol(it)=OptSol;
    %% Display Iteration Information
    if BestSol.Sol.IsFeasible && sum (BestSol.Sol.UC)== sum (Instance.r) 
        FLAG=' *';
    else
        FLAG='';
    end
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it)) FLAG '; Fittness Evaluations = ' num2str(Counter(it))]);
    figure(1);
    PlotSolution(BestSol.Sol,Instance);
    pause(0.01);
end
%% Results
figure;
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');

