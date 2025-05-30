function model = SMS_fitModel(modelParameters ,ITI, Asyn ,IOI, mod_method,TestRange)

alpha = NaN;
beta = NaN;
delta = NaN;
phi = NaN;
sT = NaN;
sM = NaN;
LLE = NaN;


    
switch mod_method

    case 'nori'    %Noris function

        
        ITER = modelParameters.ITER;
        Kratio = modelParameters.Kratio;
        
        
       [alpha, beta, sT, sM] = compute_bGLS3_version27OCT_multiperson( ...
        ITI(TestRange),Asyn(TestRange),ITER,Kratio);
    
    case 'adapt'   %Mariekesadaptationmaodel

        ITER = modelParameters.ITER;
        Kratio = modelParameters.Kratio;
        
        
        model(1).parameters = 'adapt';
       [alpha, beta, sT, sM, LLE ] = Estimates_AdaptationModel( ...
        ITI(TestRange),Asyn(TestRange),ITER,Kratio);       

    
    case 'hybrid' 


       [delta, alpha, sT, sM, LLE] = Estimates_HybridModel( ...
        IOI(TestRange),ITI(TestRange),Asyn(TestRange));                   

    case 'jointAlpha'    % MAriekes Joint Model Alpha   %fixed spelling

        L = modelParameters.LowBound;
        H = modelParameters.UpBound;
        

       [phi,delta,alpha,sT,sM,LLE] = Estimates_JointModelAlpha( ...
        IOI(TestRange),ITI(TestRange),Asyn(TestRange),L,H);   

    case 'jointBeta'     %Mariekes Joint Model Beta 

        L = modelParameters.LowBound;
        H = modelParameters.UpBound;

        [phi,delta,beta,sT,sM,LLE] = Estimates_JointModelBeta( ...
        IOI(TestRange),ITI(TestRange),Asyn(TestRange),L,H);

end

model = table( alpha, beta, delta, phi, sT, sM, LLE);


