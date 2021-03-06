function [mte,varargout] = lesMTE(r,K,sigma,dt,numtrials,x0)



if ~exist('x0', 'var') || isempty(x0)
    x0 = 10;
end



rng('shuffle');

Tmax = 1000;
bigstep = round(Tmax/dt);
%pop0 = 10;%sqrt(K);

te_vec = zeros(1,numtrials);
popthresh = 1;

for n =1:numtrials
    %if mod(n,10) ==0
        %disp(num2str(n))
    %end
    
    lextinct = 0;
    thispop = x0;
    s = 0;
    random_numbers = randn(bigstep,1);

    while lextinct == 0
        if thispop <= popthresh
            te_vec(n) = s*dt;
            lextinct = 1;    
        end
        
        s = s+1;
        
        if s > bigstep
            %te_vec(n) = bigstep*dt;
            %disp('Long run')
            %break
            this_rn = randn;
        else
            this_rn = random_numbers(s);
        end
        
        thispop =  thispop + dt.*r.*thispop.*(1-thispop./K) + sigma.*thispop.*this_rn.*sqrt(dt) + .5.*thispop.*sigma.^2.*dt.*(this_rn.^2-1);
        
        

    end

end

mte = mean(te_vec);
varargout(1) = {te_vec};

end

