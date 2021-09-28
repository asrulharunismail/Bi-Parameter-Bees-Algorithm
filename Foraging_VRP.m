function qnew=Foraging_VRP(q, as)
%as=assignment
    m=randi([1 3]);
    
    switch m
        case 1
            % Do Swap
            qnew=Swap(q,as);
            
        case 3
            % Do Reversion
            qnew=Reversion(q,as);
            
        case 2
            % Do Insertion
            qnew=Insertion(q,as);
    end

end

function qnew=Swap(q,as)

    n=numel(q);
    ass = D_Tri_dis(1,as,n);
    i=randi(n);
    i1=i;
    i2=i(1)+ass;
    if i2>n
        i2=n;
    end
        
    qnew=q;
    qnew([i1 i2])=q([i2 i1]);
    
end

function qnew=Reversion(q,as)

    n=numel(q);
    
    n=numel(q);
    ass = D_Tri_dis(1,as,n);
    i=randi(n);
    i1=i;
    i2=i(1)+ass;
    if i2>n
        i2=n;
    end
    
    qnew=q;
    qnew(i1:i2)=q(i2:-1:i1);

end

function qnew=Insertion(q,as)

    n=numel(q);
    ass = D_Tri_dis(1,as,n);
    i=randi(n);
    i1=i;
    i2=i(1)+ass;
    if i2>n
        i2=n;
    end
    
    if i1<i2
        qnew=[q(1:i1-1) q(i1+1:i2) q(i1) q(i2+1:end)];
    else
        qnew=[q(1:i2) q(i1) q(i2+1:i1-1) q(i1+1:end)];
    end

end

function [ angka ] = D_Tri_dis(k,t,b)
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
    angka=randi([round(t-a),round(t+b)]);
end

function angka=lapis2(t,a,b)
    angka=randi([round(t-2*a),round(t+2*b)]);
end

function angka=lapis3(t,a,b)
    angka=randi([round(t-3*a),round(t+3*b)]);
end

function angka=lapis4(t,a,b)
    angka=randi([round(t-4*a),round(t+4*b)]);
end

function angka=lapis5(t,a,b)
    angka=randi([round(t-5*a),round(t+5*b)]);
end

function angka=lapis6(t,a,b)
    angka=randi([round(t-6*a),round(t+6*b)]);
end

function angka=lapis7(t,a,b)
    angka=randi([round(t-7*a),round(t+7*b)]);
end

function angka=lapis8(t,a,b)
    angka=randi([round(t-8*a),round(t+8*b)]);
end

function angka=lapis9(t,a,b)
    angka=randi([round(t-9*a),round(t+9*b)]);
end

function angka=lapis10(t,a,b)
    angka=randi([round(t-10*a),round(t+10*b)]);
end