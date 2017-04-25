function [mte] = lrcMTE(r,K,lambda,f,dt,numtrials,lplot)

rng('shuffle');

bigstep = round(10000/dt);
%random_numbers = rand(bigstep,numtrials);
pop0 = 10;%sqrt(K);

te_vec = zeros(1,numtrials);

for n =1:numtrials
    lextinct = 0;
    thispop = pop0;
    s = 0;
    random_numbers = rand(bigstep,1);

    while lextinct == 0
        s = s+1;
        
        if s > bigstep
            %te_vec(s) = dt*bigstep;
            %disp('Long wait')
            %break
            this_rn = rand;
        else
            this_rn = random_numbers(s);
        end
        
        thispop = rcUpdate(thispop,dt,r,K,lambda,f,this_rn);
        
        if thispop <= 1
            te_vec(n) = s*dt;
            lextinct = 1;    
        end

    end

end

mte = mean(te_vec);

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