function [ z, L,sol] = BiBA_Clustering_using_Coverage( model )
I=model.I;
J=model.J;
C=model.c;
X=[];
X(1,:)= model.x;
X(2,:)= model.y;
X=transpose(X);
[idx,~] = BiBA_Clust_VRP( X ,J ,model );
q=[];

for i=1:J
    q{i,1}=transpose(find(idx==i));
end
for j=1: J-1
    q{j}=[q{j} I+j];
end

z=q{1};
for j=2: J
    z =[z q{j}];
end

[L,sol] = MyCost_VRP(z,model);

end

function [ id,pos ] = BiBA_Clust_VRP( X,k,model )
    CostFunction=@(m) ClusteringCost_VRP(m, X, model);     % Cost Function
    VarMin= repmat(min(X),k,1);      % Lower Bound of Variables
    VarMax= repmat(max(X),k,1);      % Upper Bound of Variables
    range=VarMax(1)-VarMin(1);
    %% SBA
    MaxEval = 5000; %%%sebelumnya 5000
    n=10;
    nep =40;
    recruitment = round(linspace(nep,1,n));
    ColonySize=sum(recruitment);  
    MaxIt=round(MaxEval/ColonySize);
    %% Initialization
    Empty_Bees.Position=[];
    Empty_Bees.Cost=[];
    Empty_Bees.Out=[];
    Empty_Bees.Size=[];
    Empty_Bees.counter=[];
    Bees=repmat(Empty_Bees,n,1);
    counter=0;
    rad=rand()*abs((max(model.x)-model.x0));
    
    % Generate Initial Solutions
    for i=1:n
        Bees(i).Position=Coverage_Init(model.x0,model.y0,rad,k);
        [Bees(i).Cost,Bees(i).Out]=CostFunction(Bees(i).Position);
        Bees(i).Size = range;
        counter=counter+1;
        Bees(i).counter= counter;
    end
    sz= linspace(0,1,n);
    %% Sites Selection 
    [~, RankOrder]=sort([Bees.Cost]);
    Bees=Bees(RankOrder);
    P=1;
    BestCost=zeros(MaxIt,1);
    %% Bees Algorithm Local and Global Search
    for it=1:MaxIt
        if counter >= MaxEval
            break;
        end
        % All Sites (Exploitation and Exploration)
        for i=1:n
            bestnewbee.Cost=inf;
            assigntment=D_Tri_real_array(0,sz(i),1,1,recruitment(i));
            for j=1:recruitment(i)
                if P==1
                    newbee.Position= Integrated_Foraging_stlim_unif(Bees(i).Position,assigntment(j),VarMax(1),VarMin(1),Bees(i).Size);
                else
                    newbee.Position= Integrated_Foraging_stlim(Bees(i).Position,assigntment(j),VarMax(1),VarMin(1),Bees(i).Size);
                end
                [newbee.Cost,newbee.Out]=CostFunction(newbee.Position);
                newbee.Size= Bees(i).Size;
                counter=counter+1;
                newbee.counter= counter;
                if newbee.Cost<bestnewbee.Cost
                    bestnewbee=newbee;
                end
            end
            if bestnewbee.Cost<Bees(i).Cost
                Bees(i)=bestnewbee;
            end
        end
        % SORTING
        [~, RankOrder]=sort([Bees.Cost]);
        Bees=Bees(RankOrder);
        % Update Best Solution Ever Found
        OptSol=Bees(1);
        BestCost(it)=OptSol.Cost;
    end
    figure(1);
    PlotSolution_C(X, OptSol);
    pause(0.01);
    id=OptSol.Out.ind;
    pos=OptSol.Position;
end

function [z, out] = ClusteringCost_VRP(m, X, model)
    J=model.J;
    C=model.r;
    maxC=model.c(1);
    
    d = pdist2(X, m);
    
    [dmin, ind] = min(d, [], 2);
    
    WCD = sum(dmin);
    
        [~, SortOrder]=sort([ind]);
    dmin=dmin(SortOrder);
    tes=zeros(1,J);
    VC=zeros(1,J);
    for i=1:J
        tes(i)=sum(C(find(ind==i)));
    end
    for j=1:J
        if tes(j) <= maxC
            VC(j)=0;
        else
            VC(j)=(tes(j) - maxC)*10;
        end
    end
    WVC=sum(VC);
    
    z=WCD+WVC*10^2;

    out.d=d;
    out.dmin=dmin;
    out.ind=ind;
    out.WCD=WCD; 
end


function z = Coverage_Init(x,y,r,n1)
hold on
n=n1/2;
RA=rand();
th = 0+RA:pi/n:2*pi+RA;

for i=1:n1
    for j=1:2
        if j==1
            z(i,j) = r * cos(th(i)) + x;
        else
            z(i,j) = r * sin(th(i)) + y;
        end
    end    
end

h = plot(z(:,1), z(:,2), 'o');
hold off
end

function y=Integrated_Foraging_stlim_unif(x,ass,Vmx,Vmn,size)
    r=ass*size;
    
    y=x; 
    
    y = y + (random('unif',-r,r)); %.*pert);
    y(y>Vmx)=Vmx;
    y(y<Vmn)=Vmn;
end


function y=Integrated_Foraging_stlim(x,ass,Vmx,Vmn,size)
    r=ass*size;
    
    nVar=numel(x);
    
    k=randi([1 nVar]);
  
    y=x;
    y(k)=y(k)+ r*((-1)^randi(2));
    y(y>Vmx)=Vmx;
    y(y<Vmn)=Vmn; 
end

function [ M ] = D_Tri_real_array(k,t,b,baris,kolom)
    M=zeros(baris,kolom);
    for i=1:baris
        for j=1:kolom
            M(i,j)=D_Tri_real(k,t,b);
        end
    end
     
end

function [ angka ] = D_Tri_real(k,t,b)
m=randi([1 10]);
    a=(t-k)/10;
    b=(b-t)/10;
    switch m
        case 1
            angka=lapis1(t,a,b);
        case 2
            angka=lapis2(t,a,b);
        case 3
            angka=lapis3(t,a,b);
        case 4
            angka=lapis4(t,a,b);
        case 5
            angka=lapis5(t,a,b);
        case 6
            angka=lapis6(t,a,b);
        case 7
            angka=lapis7(t,a,b);
        case 8
            angka=lapis8(t,a,b);
        case 9
            angka=lapis9(t,a,b);
        case 10
            angka=lapis10(t,a,b);
    end
end

function angka=lapis1(t,a,b)
    angka=unifrnd((t-a),(t+b),1);
end

function angka=lapis2(t,a,b)
    angka=unifrnd((t-2*a),(t+2*b),1);
end

function angka=lapis3(t,a,b)
    angka=unifrnd((t-3*a),(t+3*b),1);
end

function angka=lapis4(t,a,b)
    angka=unifrnd((t-4*a),(t+4*b),1);
end

function angka=lapis5(t,a,b)
    angka=unifrnd((t-5*a),(t+5*b),1);
end

function angka=lapis6(t,a,b)
    angka=unifrnd((t-6*a),(t+6*b),1);
end

function angka=lapis7(t,a,b)
    angka=unifrnd((t-7*a),(t+7*b),1);
end

function angka=lapis8(t,a,b)
    angka=unifrnd((t-8*a),(t+8*b),1);
end

function angka=lapis9(t,a,b)
    angka=unifrnd((t-9*a),(t+9*b),1);
end

function angka=lapis10(t,a,b)
    angka=unifrnd((t-10*a),(t+10*b),1);
end

function PlotSolution_C(X, sol)
    m = sol.Position;
    k = size(m,1);
    ind = sol.Out.ind;
    Colors = hsv(k);
    for j=1:k
        Xj = X(ind==j,:);
        plot(Xj(:,1),Xj(:,2),'x','LineWidth',1,'Color',Colors(j,:));
        hold on;
    end
    plot(m(:,1),m(:,2),'ok','LineWidth',2,'MarkerSize',12);
    hold off;
    grid on;
end