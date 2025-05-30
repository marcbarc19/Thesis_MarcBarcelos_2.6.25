function [alphas,betas,st,sm]=compute_bGLS3_version27OCT_multiperson(r,es,ITER,Kratio)
ITER=20; % PK added
Kratio=1; % PK added

if isempty(ITER)    
    ITER=20;
end
if isempty(Kratio)
    Kratio=1;
end
TRESH=1e-3;

N=size(r,1)-2;
P=size(es,2);
assert(size(r,1)==size(es,1));


for p=1:P,
    pos=~isnan(es(:,p));
    es(:,p)=es(:,p)-mean(es(pos,p));
    es(~pos,p)=0;
end

%es=es-repmat(mean(es),[size(es,1),1]);


    
    b3=r(3:end)-r(2:end-1);
    A3=[es(2:end-1,:),es(1:end-2,:)];
    posr=(~isnan(r(3:end))) & (~isnan(r(2:end-1)));
    assert(sum(isnan(b3(posr)))==0);
    
    K11=2;
    K12=-1;
    K13=0;

    zold=zeros(2*P,1)-9999;
    for iter=1:ITER,
        CC=diag(K11*ones(1,N),0)+ diag(K12*ones(1,N-1),1) + diag(K12*ones(1,N-1),-1) +diag(K13*ones(1,N-2),2) + diag(K13*ones(1,N-2),-2);
        iC=inv(CC);
        z=inv((A3(posr,:)')*iC(posr,posr)*A3(posr,:))*((A3(posr,:)')*iC(posr,posr)*b3(posr));
        d=A3*z-b3;
        d(~posr)=mean(d(posr));
        
        K=cov(d(1:(end-1)),d(2:end));    
        K11=(K(1,1)+K(2,2))/2;
        K12=K(1,2);

        if K12<(-(4+Kratio)/(2*Kratio+6))*K11
            K12=(-(4+Kratio)/(2*Kratio+6))*K11;
        end
        if K12>(-1/2)*K11
            K12=(-1/2)*K11;
        end

        K13=(K11+2*K12)/(-2);
        
        
        if (sum(abs(z-zold))<TRESH)
            break;
        end
        zold=z;
        
    end
    alphas=z((P+1):(P+P));
    betas=-(z(1:P)+z((P+1):(P+P)));
    sm=sqrt((K11+2*K12)/(-2));
    st=sqrt((2*K11+3*K12));
    %totstd=sqrt(K11);

    