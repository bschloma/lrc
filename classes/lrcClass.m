classdef lrcClass
    
    properties
        
        savename = '';
        savedir = '';
        params;
        popmat;
        moments;
        momentSweep = struct();
        dis = struct();
        
    end
    
    methods
        
        function obj = lrcClass(varargin)
            
            
    
            if nargin == 0
                c= clock;
                runname = '';
                for cc = 1:5   %leave off seconds
                    runname = [runname num2str(c(cc)) '_'];
                end
                
                runname = runname(1:(end-1));    %remove trailing _
                
                obj.savename= ['lrc_' runname];
                
                obj.savedir = pwd;
            elseif nargin == 1
                obj.savename = varargin{1};
                obj.savedir = pwd;
            elseif nargin == 2
                obj.savename = varargin{1};
                obj.savedir = varargin{2};
            else
                disp('lrcClass takes at most 2 inputs')
            end
            
            
            obj.params = lrcParamsClass(obj.savename,obj.savedir);
            
        end
        
        function obj = computeSamplePaths(obj)

            %params = obj.params;

            %% Unpack params locally
            mu = obj.params.mu;
            %meanK = params.meanK;
            %sigK = params.sigK;
            Kparams = obj.params.Kparams;
            Tmax = obj.params.Tmax;
            dt = obj.params.dt;
            numtrials = obj.params.numtrials;
            lplot = obj.params.lplot;
            poiscolor =obj.params.poiscolor;
            bcolor = obj.params.bcolor ;
            f = obj.params.f;
            lambda = obj.params.lambda;
            lextinct = obj.params.lextinct;
            M = obj.params.M;
            
            
            %% Compute
            
            obj.popmat = lrcPaths(mu,Kparams,f, lambda, numtrials, Tmax, dt,lextinct);
            tvec = 0:dt:Tmax;
            
            %% Plot
            if lplot
                
                figure; hold on;
                semilogy(mu.*tvec,obj.popmat(:,1:10),'linewidth',3);
                set(gca,'fontsize',24,'linewidth',4,'yscale','log')
                xlabel('rt','fontsize',24);
                ylabel('X_t','fontsize',24);
                axis([0 Tmax 1e-1 10^(log10(Kparams(1))+1)]);
                
            end
            
            
        end
        
        function obj = computeMoments(obj)
            
           
            
            %% Unpack params locally
            mu = obj.params.mu;
            %meanK = params.meanK;
            %sigK = params.sigK;
            Kparams = obj.params.Kparams;
            Tmax = obj.params.Tmax;
            dt = obj.params.dt;
            numtrials = obj.params.numtrials;
            lplot = obj.params.lplot;
            poiscolor = obj.params.poiscolor;
            bcolor = obj.params.bcolor ;
            f = obj.params.f;
            lambda = obj.params.lambda;
            lextinct = obj.params.lextinct;
            M = obj.params.M;
            
            
            %% Compute
            
            obj.moments = lrcMoments(mu,Kparams,f, lambda, numtrials, Tmax, dt,lextinct,M);
            
            %% Plot
            if lplot
                
                tvec = 0:dt:Tmax;

                %mvec = 1:M;
                scaled_moments = scaleByPowersOfK(obj.moments,Kparams(1));
                
                figure; hold on;
                plot(mu.*tvec,scaled_moments,'linewidth',3);
                set(gca,'fontsize',24,'linewidth',4)
                xlabel('rt','fontsize',24);
                ylabel('E[X^n_t/K^n]','fontsize',24);
                %axis([0 Tmax 0 .5);
                
                legendcell = cell(1,M);
                for mm =1:M
                    legendcell{mm} = ['n = ' num2str(mm)];
                end
                
                legend(legendcell,'fontsize',18,'location','se');
                
            end

        end
    function obj = computeMomentsSweepParams(obj)

       
        
        %% Unpack params locally
        mu = obj.params.mu;
        %meanK = params.meanK;
        %sigK = params.sigK;
        Kparams = obj.params.Kparams;
        Tmax = obj.params.Tmax;
        dt = obj.params.dt;
        numtrials = obj.params.numtrials;
        lplot = obj.params.lplot;
        poiscolor =obj.params.poiscolor;
        bcolor = obj.params.bcolor ;
        f = obj.params.f;
        lambda = obj.params.lambda;
        lextinct = obj.params.lextinct;
        M = obj.params.M;
        
        obj.params = obj.params.initZSweep();
        
        zmin = obj.params.zsweep.zmin;
        zmax = obj.params.zsweep.zmax;
        numzs = obj.params.zsweep.numzs;
        
        zmin_an = obj.params.zsweep.zmin_an;
        zmax_an = obj.params.zsweep.zmax_an;
        numzs_an = obj.params.zsweep.numzs_an;
        numsteps = length(0:obj.params.dt:obj.params.Tmax);
        
        %% Arrays
        zarray = linspace(zmin,zmax,numzs);
        larray = -zarray./log(f);
        farray = exp(-zarray./lambda);
        zarray_an = linspace(zmin_an,zmax_an,numzs_an);
        larray_an = -zarray_an./log(f);
        farray_an = exp(-zarray_an./lambda);
%         momentsL = zeros(numsteps,M,numel(larray));
%         momentsF = zeros(numsteps,M,numel(larray));
        statMomL = zeros(numzs_an,M);
        statMomF = zeros(numzs_an,M);
        
        %% Compute
        % sweep lambda
        for j = 1:numzs
            disp(['j = ' num2str(j)])
            
            obj.momentSweep.momentsL(:,:,j) = lrcMoments(mu,Kparams,f,larray(j),numtrials,Tmax,dt,lextinct,M);
            
        end
        
        % sweep f
        for k = 1:numzs
            disp(['k = ' num2str(k)])
            
            obj.momentSweep.momentsF(:,:,k) = lrcMoments(mu,Kparams,farray(k),lambda,numtrials,Tmax,dt,lextinct,M);
            
        end
        
        
        %% Plot
        if lplot
            
            meanEndL = reshape(obj.momentSweep.momentsL(end,1,:),1,numzs);
            varEndL = reshape(obj.momentSweep.momentsL(end,2,:) - obj.momentSweep.momentsL(end,1,:).^2,1,numzs);
            
            meanEndF = reshape(obj.momentSweep.momentsF(end,1,:),1,numzs);
            varEndF = reshape(obj.momentSweep.momentsF(end,2,:) - obj.momentSweep.momentsF(end,1,:).^2,1,numzs);
            
            for jj = 1:numzs_an
                statMomL(jj,:) = lrcExactStationaryMoments(mu,Kparams(1),larray_an(jj),f,M);
            end
            
            for kk = 1:numzs_an
                statMomF(kk,:) = lrcExactStationaryMoments(mu,Kparams(1),lambda,farray_an(kk),M);
            end
            
            meanEndL_an = statMomL(:,1);
            meanEndF_an = statMomF(:,1);
            
            varEndL_an = statMomL(:,2) - meanEndL_an.^2;
            varEndF_an = statMomF(:,2) - meanEndF_an.^2;
            
            
            figure; hold on;
            plot(meanEndL./Kparams(1),varEndL./Kparams(1)./Kparams(1),'kd','markersize',24,'markerfacecolor',poiscolor);
            plot(meanEndF./Kparams(1),varEndF./Kparams(1)./Kparams(1),'ko','markersize',24,'markerfacecolor',poiscolor);
            plot(meanEndL_an./Kparams(1),varEndL_an./Kparams(1)./Kparams(1),'k-','linewidth',4)
            plot(meanEndF_an./Kparams(1),varEndF_an./Kparams(1)./Kparams(1),'k-','linewidth',4)
            
            set(gca,'fontsize',24,'linewidth',4)
            xlabel('E[X]/\kappa','fontsize',24);
            ylabel('Var[X]/\kappa^2','fontsize',24);
            axis([0 1 0 .3])
            
        end
    end
    
    function obj = disTrafo(obj)
        
       
        
        %% Unpack params locally
        mu = obj.params.mu;
        Kparams = obj.params.Kparams;
        K = Kparams(1);
        Tmax = obj.params.Tmax;
        dt = obj.params.dt;
        numtrials = obj.params.numtrials;
        f = obj.params.f;
        lambda = obj.params.lambda;
        lextinct = obj.params.lextinct;
        lplot = obj.params.lplot;
        poiscolor = obj.params.poiscolor;
        bcolor = obj.params.bcolor;
        sigma = obj.params.sigma;
        M = obj.params.M;

        if numel(fieldnames(obj.params.distrafo))==0
            obj.params = obj.params.initDisTrafo();
        end
        
        numbins = obj.params.distrafo.numbins;
        numalphs = obj.params.distrafo.numalphs;
        lmultmax = obj.params.distrafo.lmultmax;
        ymax = obj.params.distrafo.ymax;
        
        
        
        %% Arrays
        %scalevec = linspace(1,lmultmax,numalphs);
        %scalevec = linspace(1,sqrt(lmultmax),numalphs); 
        %scalevec = scalevec.^2;
        
        scalevec = logspace(0,log10(lmultmax),numalphs);
        obj.dis.scalevec = scalevec;

        lvec = lambda.*scalevec; 
        fvec = exp(-sigma./sqrt(lvec)); 
        rvec = mu -lvec.*log(fvec) - lvec.*log(fvec).*log(fvec)./2; 
        Kvec = K.*(1-lvec.*log(fvec)./mu - lvec.*log(fvec).*log(fvec)./2./mu); 
        dtvec = dt./sqrt(scalevec);    %adaptive timestep
        
        obj.dis.lvec = lvec;
        obj.dis.fvec = fvec;
        obj.dis.rvec = rvec;
        obj.dis.Kvec = Kvec;
        obj.dis.dtvec = dtvec;
        mendvec = zeros(1,numalphs);
        varvec = zeros(1,numalphs);
        
        obj.dis.lrcHist = zeros(numalphs,numbins);
        obj.dis.lesHist = zeros(1,numbins);
        obj.dis.lrcMeans = zeros(numalphs,1);
        obj.dis.lrcVars = zeros(numalphs,1);
        obj.dis.dkl = zeros(numalphs,1);
        
        rB = mu;
        KB = K;
        dtB = obj.params.dt;
        
        bins = linspace(0, 5, numbins);
        [obj.dis.lesHist,obj.dis.bins,obj.dis.lesMeans,obj.dis.lesVars] = computeLESendHist(rB,KB,sigma,dtB,Tmax,numtrials,lextinct,bins,false);
        lesExactMoms = lesExactStationaryMoments(rB,KB,sigma,2);
        obj.dis.lesExactMoms = lesExactMoms;
        obj.dis.lrcExactMoms = zeros(numalphs,M);
        
        
        
        %% Loop through alphas
        for s = 1:numalphs

            [pophist,~,mendvec(s),varvec(s)] = computeLRCendHist(rvec(s),Kvec(s),lvec(s),fvec(s),dtvec(s),Tmax,numtrials,lextinct,bins,false);
                       
            obj.dis.lrcHist(s,:) = pophist;
            obj.dis.lrcMeans(s) = mendvec(s);
            obj.dis.lrcVars(s) = varvec(s);
            
            obj.dis.dkl(s) = DKL(pophist,obj.dis.lesHist);
            obj.dis.lrcExactMoms(s,:) = lrcExactStationaryMoments(rvec(s),Kvec(s),lvec(s),fvec(s),M);
                        
            
        end
        
        if lplot
            obj.plotDisTrafo();
        end
        
        
        
    end
    
    function obj = plotDisTrafo(obj)

        if numel(fieldnames(obj.params.distrafo))==0
            disp('Error in plotDisTrafo:  Need to compute disTrafo first')
            return
        end
        
        % Probability density evolution
        figure; hold on;
        bins = obj.dis.bins;

        
        for a = 1:obj.params.distrafo.numalphs
            ax_s = subplot(2,3,a);
            bar(ax_s,bins,obj.dis.lrcHist(a,:)./sum(obj.dis.lrcHist(a,:))./diff(obj.dis.bins(2:3)),'facecolor',obj.params.poiscolor,'edgecolor',obj.params.poiscolor)
               
            set(gca,'fontsize',22,'linewidth',3,'ticklength',[.02;.02])
            
            if a > 2
                axis([0 5 0 obj.params.distrafo.ymax])
            else
                axis([0 5 0 3])
            end
            
        end
        
        a = a+1;
        ax_s = subplot(2,3,a);

        bar(ax_s,bins,obj.dis.lesHist./sum(obj.dis.lesHist)./diff(obj.dis.bins(2:3)),'facecolor',obj.params.bcolor,'edgecolor',obj.params.bcolor)            
        set(gca,'fontsize',22,'linewidth',3,'ticklength',[.02;.02])
                
        axis([0 5 0 obj.params.distrafo.ymax])

        % Mean


        xmin = .8; xmax_mult = 7; numxs = 10;
        xvec = linspace(xmin,xmax_mult*obj.dis.scalevec(end),numxs);
        K = obj.params.Kparams(1);
        
        
        figure; hold on;
        plot(obj.dis.scalevec,obj.dis.lrcMeans./K,'ko','markerfacecolor',obj.params.poiscolor,'markersize',22,'linewidth',2)
        plot(obj.dis.scalevec,obj.dis.lrcExactMoms(:,1)./K,'color',obj.params.poiscolor,'linewidth',4)
        plot(xmax_mult*obj.dis.scalevec(end),obj.dis.lesMeans./K,'ks','markerfacecolor',obj.params.bcolor,'linewidth',2,'markersize',22)
        plot(xvec,obj.dis.lesExactMoms(1).*ones(size(xvec))./K,'color',obj.params.bcolor,'linewidth',4,'linestyle','--')

        set(gca,'fontsize',22,'xscale','log','linewidth',4)
        xlabel('\alpha','fontsize',22)
        ylabel('E[X]/K','fontsize',22)
        axis([xmin xmax_mult*obj.dis.scalevec(end) 0 1]);
        
        % Variance
        %vend_theory = Kvec.^2.*(-log(fvec)-(1-fvec)).*lvec./rvec.*(1+lvec.*log(fvec)./rvec);
        %vend_theory_lim = Kvec.^2.*(log(fvec).^2).*lvec./2./rvec.*(1+lvec.*log(fvec)./rvec);
        lrc_var_theory = obj.dis.lrcExactMoms(:,2) - obj.dis.lrcExactMoms(:,1).^2;
        les_var_theory =obj.dis.lesExactMoms(2) - obj.dis.lesExactMoms(1)^2;
        figure; hold on;
        plot(obj.dis.scalevec,obj.dis.lrcVars./K^2,'ko','markerfacecolor',obj.params.poiscolor,'linewidth',2,'markersize',22)
        plot(7*obj.dis.scalevec(end),obj.dis.lesVars./K^2,'ks','markerfacecolor',obj.params.bcolor,'linewidth',2,'markersize',22)
        plot(obj.dis.scalevec,lrc_var_theory/K^2,'color',obj.params.poiscolor,'linewidth',4)

        plot(xvec,les_var_theory./K^2.*ones(size(xvec)),'color',obj.params.bcolor,'linewidth',4,'linestyle','--')
       
        set(gca,'fontsize',22,'linewidth',4)
        xlabel('\alpha','fontsize',22)
        ylabel('Var[X]/K^2','fontsize',22)
        set(gca,'xscale','log')
        axis([.8 7*obj.dis.scalevec(end) 0 .3])
        
        figure; hold on;
        plot(obj.dis.scalevec,obj.dis.dkl,'ko','markerfacecolor',[.2 .2 .2],'markersize',22,'linewidth',2)

        %plot(xmax_mult*scalevec(end),obj.dis.lesMeans./K,'ks','markerfacecolor',bcolor,'linewidth',2,'markersize',22)
        %plot(xvec,lesExactMoms(1).*ones(size(xvec)),'kd','color',bcolor,'markerfacecolor',.5*bcolor,'linewidth',1,'markersize',18)
        %plot(xvec,lesExactMoms(1).*ones(size(xvec))./K,'color',bcolor,'linewidth',4,'linestyle','--')

        set(gca,'fontsize',22,'xscale','log','linewidth',4)
        xlabel('\alpha','fontsize',22)
        ylabel('D_{KL}','fontsize',22)
        axis([xmin xmax_mult*obj.dis.scalevec(end) 0 1.2*max(obj.dis.dkl)]);

    
    end

    end

end
  