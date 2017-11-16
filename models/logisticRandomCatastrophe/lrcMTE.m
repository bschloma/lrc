function [mte,varargout] = lrcMTE(r,K,lambda,f,dt,numtrials,lplot,method,x0)

if ~exist('method', 'var') || isempty(method)
    method = 'pdmp';
end

if ~exist('x0', 'var') || isempty(x0)
    x0 = 10;
end
rng('shuffle');

switch method
    case 'pdmp'
        eta = 1.;
    case 'euler'
        eta = 0.;
    otherwise
        disp('Error in lrcMoments:  Invalid method type');
        return
end

bigstep = round(10000/dt);
%random_numbers = rand(bigstep,numtrials);
%pop0 = 10;%sqrt(K);

te_vec = zeros(1,numtrials);

for n =1:numtrials
    %if mod(n,10)==0
    %    disp(num2str(n))
    %end
    lextinct = 0;
    thispop = x0;
    s = 0;
    random_numbers = rand(bigstep,1);

    while lextinct == 0
        if thispop <= 1
            te_vec(n) = s*dt;
            lextinct = 1;    
        end
        
        s = s+1;
        
        if s > bigstep
            %te_vec(s) = dt*bigstep;
            %disp('Long wait')
            %break
            this_rn = rand;
        else
            this_rn = random_numbers(s);
        end
        
        dNt = this_rn < lambda*dt;
        
        thispop = (1-dNt).*thispop + (1-eta.*dNt).*dt*r.*thispop.*(1-thispop./K) + dNt.*f.*thispop;

        %thispop = rcUpdate(thispop,dt,r,K,lambda,f,this_rn);
        
        

    end

end

mte = mean(te_vec);
varargout(1) = {te_vec};

if lplot
    figure;
    hold on;
    nbins =100;
    hist(te_vec,nbins);
    set(gca,'fontsize',22)
    xlabel('Extinction Time','fontsize',22)
    ylabel('Probability','fontsize',22)
    
end

end

function [x] = rcUpdate(x,dt,r,K,lambda,f,random_number)

if random_number < lambda*dt
    x = f*x;
else
    x = x+ dt*x*r*(1-x/K);
end

end