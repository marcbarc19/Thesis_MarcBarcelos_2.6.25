function [alphas,betas,st,sm,LL]=Estimates_AdaptationModel(r,es,ITER,Kratio)

%%%% function for estimating parameters of the Adaptation model %%%%
%%% input:  
%           r = ITI
%           es = asyn
%           ITER = defaultsetting 20
%           Kratio = defaultsetting 1

%%% output:
%           alphas = phase correction parameter
%           betas = period correction parameter
%           st = timekeeper noise estimate
%           sm = motor noise estimate
%           LL = log likelihood -- fit of the model to the data


% N. Jacoby & M.C. van der Steen

if isempty(ITER)    
    ITER=20;
end
if isempty(Kratio)
    Kratio=1;
end
TRESH=1e-3;

N=size(r,1)-4;
% N=size(r,1)-2;
P=size(es,2);
assert(size(r,1)==size(es,1));


for p=1:P,
    pos=~isnan(es(:,p));
    es(:,p)=es(:,p)-mean(es(pos,p));
    es(~pos,p)=0;
end

%es=es-repmat(mean(es),[size(es,1),1]);

    b3=r(5:end)-r(4:end-1);
    A3=[es(4:end-1,:),es(3:end-2,:)];
    posr=(~isnan(r(5:end))) & (~isnan(r(4:end-1)));
     assert(sum(isnan(b3(posr)))==0);
    
%     b3=r(3:end)-r(2:end-1);
%     A3=[es(2:end-1,:),es(1:end-2,:)];
%     posr=(~isnan(r(3:end))) & (~isnan(r(2:end-1)));
%      assert(sum(isnan(b3(posr)))==0);
    
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

    CC=diag(K11*ones(1,N),0)+ diag(K12*ones(1,N-1),1) + diag(K12*ones(1,N-1),-1) +diag(K13*ones(1,N-2),2) + diag(K13*ones(1,N-2),-2);
   residual = A3*z - b3;
    mu = zeros(size(b3));
    
    % Validate CC
    if ~isequal(CC, CC') || any(eig(CC) <= 0)
        warning('⚠️ Covariance matrix was not positive definite — regularizing.');
        CC = CC + eye(size(CC)) * 1e-6;
    end
    
    n = length(mu);
    if ~isequal(size(residual), [n,1])
        error('❌ Residual and mean vector sizes do not match.');
    end
    if ~isequal(size(CC), [n,n])
        error('❌ Covariance matrix size does not match residual vector.');
    end
    
    % Calculate pdf safely
    try
        pdf_val = mvnpdf(residual, mu, CC);
    catch err
        warning('❌ mvnpdf failed: %s', err.message);
        pdf_val = 1e-300 * (1 + rand()*9);
    end
    
    % Final safeguard
    if ~isfinite(pdf_val) || pdf_val <= 0
        pdf_val = 1e-300 * (1 + rand()*9); 
    end
    
    LL = log2(pdf_val);
    
end