function q = BRO_0_VRP( p, model, BatasB, BatasA)
I=model.I;
J=model.J;
l=randi([1 J]);
x0=model.x0;
y0=model.y0;

PTemp=p.Sol.L{l};
x=[model.x(PTemp) x0];
y=[model.y(PTemp) y0];
PTemp=[PTemp I+l];
t(:,1)=x;t(:,2)=y;
D=round(pdist2(t,t));
n=numel(PTemp);

if n<=1
    q=p.Position;
    
else
    u=1:n;
    
    idrem = randperm(n,[randi([BatasB BatasA])]);

    rem = u(idrem);
    u(idrem)=[];

    u = TwoOpt_0(u,D);

    [u,~] = Insert_Forgotten( rem,u,D,x,y);
    
    u(u==max(u))=[];

    sequence=PTemp(u);
    p.Sol.L{l}=sequence;
    for i=1:J-1
        p.Sol.L{i}=[p.Sol.L{i} I+i];
    end
    Temp=transpose(p.Sol.L);
    q=cell2mat(Temp);

end
end

function p = TwoOpt_0(p,D)
n = numel(p);
if n==0
    return
end
zmin = -1;
% Iterate until the tour is 2-optimal
while zmin < 0
    zmin = 0;
    i = 0;
    b = p(n);
    % Loop over all edge pairs (ab,cd)
    while i < n-2
        a = b;
        i = i+1;
        b = p(i);
        Dab = D(a,b);
        j = i+1;
        d = p(j);
        while j < n
            c = d;
            j = j+1;
            d = p(j);
            z = (D(a,c) - D(c,d)) + D(b,d) - Dab;
            % Keep best exchange
            if z < zmin
                zmin = z;
                imin = i;
                jmin = j;
            end
        end
    end
    % Apply exchange
    if zmin < 0
        p(imin:jmin-1) = p(jmin-1:-1:imin);
    end
end
end

function [p,L] = Insert_Forgotten( rem,init,D,x,y )
p=init;
nrem = numel (rem);
for i = 1 : nrem
    np = numel (p);
    
    Dis = zeros(np,1);
    center.x(1) = (x(p(1))+x(p(end)))*0.5;
    center.y(1) = (y(p(1))+y(p(end)))*0.5;
    Dis(1) = (sqrt((x(rem(i))-center.x(1))^2+(y(rem(i))-center.y(1))^2)); 
    for j=2:np
        center.x(j) = (x(p(j))+x(p(j-1)))*0.5;
        center.y(j) = (y(p(j))+y(p(j-1)))*0.5;
        Dis(j) = (sqrt((x(rem(i))-center.x(j))^2+(y(rem(i))-center.y(j))^2));
    end
   [a,b] = min (Dis);
    s = p(b);
    idx=find(p==s);
    if idx==1
        p = [rem(i) p];
    else
        p = [p(1:idx-1) rem(i) p(idx:end)];    
    end
end    
L=TourLength_0(p,D);
end

function L=TourLength_0(tour,D)

    n=numel(tour);

    tour=[tour tour(1)];
    
    L=0;
    for i=1:n
        L=L+D(tour(i),tour(i+1));
    end

end