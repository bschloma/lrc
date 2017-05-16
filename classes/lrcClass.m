% Program:  lrcClass.m
%
% Usage:    lrc = lrcClass();
%
% Summary:  
%   
%   This is the primary class for performing simulations of the
% Logistic Random Catastrophe (LRC) and Logistic Environmental
% Stochasticity Model.  Contains wrappers around standalone subroutines.
%
% Author:   Brandon Schlomann
%
% Date:     April 2017 - First written.
%



classdef lrcClass
    
    properties
        
        % base properites
        savename = '';                  % name of class instance
        savedir = '';                   % full path to save directory
        params;                         % lrcParams class, for storing params
        popmat;                         % 2D array of sample paths
        moments;                        % 2D array of moments
        cumulants;                      % 2D array of cumulants
        
        % other properties, for calculations related to paper
        momentSweep = struct();         % struct for moments, sweeping parameters
        dis = struct();                 % 'distributions', struct for distribution transformation 
        mte = struct();                 % 'mean time to extinction', struct for MTE calculations
        
        
    end
    
    methods
        
        function obj = lrcClass(varargin)
            % initialize lrcClass.  Options for 0, 1 (savename), or 2 (savename & savedir) inputs.
            
                    
            if nargin == 0
                
                % create savename based on date and time
                c= clock;
                runname = '';
                for cc = 1:5   %leave off seconds
                    runname = [runname num2str(c(cc)) '_'];
                end
                
                runname = runname(1:(end-1));    %remove trailing _
                
                obj.savename= ['lrc_' runname];
                
                % let present working directory be save directory
                obj.savedir = pwd;
            elseif nargin == 1
                
                % input only savename
                obj.savename = varargin{1};
                
                % savedir = present working directory
                obj.savedir = pwd;
            elseif nargin == 2
                
                %input both savename and savedir
                obj.savename = varargin{1};
                obj.savedir = varargin{2};
            else
                disp('lrcClass takes at most 2 inputs')
            end
            
            % initialize params class
            obj.params = lrcParamsClass(obj.savename,obj.savedir);
            
        end
        
        function obj = computeSamplePathsLRC(obj)
            
            % create sample paths from LRC model.  Acts as a wrapper around
            % lrcPaths.m
            %
            % Stored in obj.popmat, 2D array with rows = time, collumns =
            % trial number.
            
            
            % Compute           
            obj.popmat = lrcPaths(obj.params.mu,obj.params.Kparams,obj.params.f, ...
                obj.params.lambda, obj.params.numtrials, obj.params.Tmax, ...
                obj.params.dt,obj.params.lextinct,obj.params.method);
                        
            % Plot
            if obj.params.lplot              
                obj = obj.plotPathsLRC;
            end           
            
        end
        
        function obj = plotPathsLRC(obj,varargin)
            
            % plot sample LRC paths.  Takes an optional argument for the
            % number of paths to plot (popmat(:,1:numpaths).  Default is
            % the smaller of 10 and all paths.
            
            % check to ensure paths were computed first.
            if isempty(obj.popmat)
                disp('Error in plotPathsLRC:  Need to computeSamplePathsLRC first')
                return
            end
            
            if nargin==1
                % no input (only obj, internally)
                numpaths = min(size(obj.popmat,1),10);
            elseif nargin == 2
                % user input number of paths to plot
                numpaths = varargin{1};
            else
                disp('Error in plotPathsLRC:  Takes at most 1 argument');
                return
            end
            
            % compute time array for plotting
            tvec = 0:obj.params.dt:obj.params.Tmax;
            
            % plot
            figure; hold on;
            semilogy(obj.params.mu.*tvec,obj.popmat(:,1:numpaths),'linewidth',3);
            set(gca,'fontsize',24,'linewidth',4,'yscale','log')
            xlabel('rt','fontsize',24);
            ylabel('X_t','fontsize',24);
            axis([0 obj.params.Tmax 1e-1 10^(log10(obj.params.Kparams(1))+1)]);
            
        end
        
        function obj = computeMoments(obj)
          
            % Compute moments over time of LRC model.  Wrapper around
            % LRCmoments.m.  Computes first obj.params.M moments.
                        
            % Compute
            
            obj.moments = lrcMoments(obj.params.mu,obj.params.Kparams, ...
                obj.params.f, obj.params.lambda, obj.params.numtrials, obj.params.Tmax, ... 
                obj.params.dt,obj.params.lextinct,obj.params.method,obj.params.M);
            
            obj.cumulants = zeros(size(obj.moments));
            for s = 1:size(obj.moments,1)
                for n = 1:obj.params.M
                    obj.cumulants(s,n) = obj.m2c(obj.moments(s,1:n));
                end
            end
            
            % Plot
            if obj.params.lplot
                
                obj.plotMomentsLRC();
                obj.plotCumulantsLRC();
                
                
            end
            
        end
        
        function obj = plotMomentsLRC(obj)
            
            % Plot LRC moments over time, scaled by powers of K.
            
            % check to ensure paths were computed first.
            if isempty(obj.moments)
                disp('Error in plotMomentsLRC:  Need to computeMomentsLRC first')
                return
            end
            
             tvec = 0:obj.params.dt:obj.params.Tmax;
                
             %mvec = 1:M;
             scaled_moments = scaleByPowersOfK(obj.moments,obj.params.Kparams(1));
             
             figure; hold on;
             plot(obj.params.mu.*tvec,scaled_moments,'linewidth',3);
             set(gca,'fontsize',24,'linewidth',4)
             xlabel('rt','fontsize',24);
             ylabel('E[X^n_t/K^n]','fontsize',24);
             %axis([0 Tmax 0 .5);
             
             legendcell = cell(1,obj.params.M);
             for mm =1:obj.params.M
                 legendcell{mm} = ['n = ' num2str(mm)];
             end
             
             legend(legendcell,'fontsize',18,'location','se');
        end
        
        function obj = plotCumulantsLRC(obj)
            
            % Plot LRC cumulants over time, scaled by powers of K.
            
            % check to ensure paths were computed first.
            if isempty(obj.cumulants)
                disp('Error in plotCumulantsLRC:  Need to computeMomentsLRC first')
                return
            end
            
             tvec = 0:obj.params.dt:obj.params.Tmax;
                
             %mvec = 1:M;
             scaled_cumulants = scaleByPowersOfK(obj.cumulants,obj.params.Kparams(1));
             
             figure; hold on;
             plot(obj.params.mu.*tvec,scaled_cumulants,'linewidth',3);
             set(gca,'fontsize',24,'linewidth',4)
             xlabel('rt','fontsize',24);
             ylabel('C^n/K^n','fontsize',24);
             %axis([0 Tmax 0 .5);
             
             legendcell = cell(1,obj.params.M);
             for mm =1:obj.params.M
                 legendcell{mm} = ['n = ' num2str(mm)];
             end
             
             legend(legendcell,'fontsize',18,'location','se');
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
            method = obj.params.method;
            
            obj.params = obj.params.initZSweep();
            
            zmin = obj.params.zsweep.zmin;
            zmax = obj.params.zsweep.zmax;
            numzs = obj.params.zsweep.numzs;
            
            %zmin_an = obj.params.zsweep.zmin_an;
            %zmax_an = obj.params.zsweep.zmax_an;
            %numzs_an = obj.params.zsweep.numzs_an;
            numsteps = length(0:obj.params.dt:obj.params.Tmax);
            
            %% Arrays
            zarray = linspace(zmin,zmax,numzs);
            larray = -zarray./log(f);
            farray = exp(-zarray./lambda);
            %zarray_an = linspace(zmin_an,zmax_an,numzs_an);
            %larray_an = -zarray_an./log(f);
            %farray_an = exp(-zarray_an./lambda);
            %         momentsL = zeros(numsteps,M,numel(larray));
            %         momentsF = zeros(numsteps,M,numel(larray));
            %statMomL = zeros(numzs_an,M);
            %statMomF = zeros(numzs_an,M);
            
            %% Compute
            % sweep lambda
            for j = 1:numzs
                disp(['Beginning lambda # ' num2str(j) ' of ' num2str(obj.params.zsweep.numzs)])
                
                obj.momentSweep.momentsL(:,:,j) = lrcMoments(mu,Kparams,f,larray(j),numtrials,Tmax,dt,lextinct,method,M);
                
            end
            
            % sweep f
            for k = 1:numzs
                disp(['Beginning f # ' num2str(k) ' of ' num2str(obj.params.zsweep.numzs)])
                
                obj.momentSweep.momentsF(:,:,k) = lrcMoments(mu,Kparams,farray(k),lambda,numtrials,Tmax,dt,lextinct,method,M);
                
            end
            
            
            %% Plot
            if lplot
                
                obj.plotMomentSweepLRC();
                
            end
        end
        
        function obj = plotMomentSweepLRC(obj)
            statMomL = zeros(obj.params.zsweep.numzs_an,obj.params.M);
            statMomF = zeros(obj.params.zsweep.numzs_an,obj.params.M);
            zarray_an = linspace(obj.params.zsweep.zmin_an,obj.params.zsweep.zmax_an,...
                obj.params.zsweep.numzs_an);
            larray_an = -zarray_an./log(obj.params.f);
            farray_an = exp(-zarray_an./obj.params.lambda);
            
            meanEndL = reshape(obj.momentSweep.momentsL(end,1,:),1,obj.params.zsweep.numzs);
            varEndL = reshape(obj.momentSweep.momentsL(end,2,:) - obj.momentSweep.momentsL(end,1,:).^2,1,obj.params.zsweep.numzs);
            
            meanEndF = reshape(obj.momentSweep.momentsF(end,1,:),1,obj.params.zsweep.numzs);
            varEndF = reshape(obj.momentSweep.momentsF(end,2,:) - obj.momentSweep.momentsF(end,1,:).^2,1,obj.params.zsweep.numzs);
            
            for jj = 1:obj.params.zsweep.numzs_an
                statMomL(jj,:) = lrcExactStationaryMoments(obj.params.mu,obj.params.Kparams(1),...
                    larray_an(jj),obj.params.f,obj.params.M);
            end
            
            for kk = 1:obj.params.zsweep.numzs_an
                statMomF(kk,:) = lrcExactStationaryMoments(obj.params.mu,obj.params.Kparams(1),...
                    obj.params.lambda,farray_an(kk),obj.params.M);
            end
            
            meanEndL_an = statMomL(:,1);
            meanEndF_an = statMomF(:,1);
            
            varEndL_an = statMomL(:,2) - meanEndL_an.^2;
            varEndF_an = statMomF(:,2) - meanEndF_an.^2;
            
            
            figure; hold on;
            plot(meanEndL./obj.params.Kparams(1),varEndL./obj.params.Kparams(1)./obj.params.Kparams(1),'kd','markersize',24,'markerfacecolor',obj.params.poiscolor);
            plot(meanEndF./obj.params.Kparams(1),varEndF./obj.params.Kparams(1)./obj.params.Kparams(1),'ko','markersize',24,'markerfacecolor',obj.params.poiscolor);
            plot(meanEndL_an./obj.params.Kparams(1),varEndL_an./obj.params.Kparams(1)./obj.params.Kparams(1),'k-','linewidth',4)
            plot(meanEndF_an./obj.params.Kparams(1),varEndF_an./obj.params.Kparams(1)./obj.params.Kparams(1),'k-','linewidth',4)
            
            set(gca,'fontsize',24,'linewidth',4)
            xlabel('E[X]/\kappa','fontsize',24);
            ylabel('Var[X]/\kappa^2','fontsize',24);
            axis([0 1 0 .3])
            
        end
        
         
       function obj = computeMomentsSweepParamsLES(obj)
            
            
            
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
            method = obj.params.method;
            
            if isempty(obj.params.zsweep)
                obj.params = obj.params.initZSweep();
            end
            
            zmin = obj.params.zsweep.zmin;
            zmax = obj.params.zsweep.zmax;
            numzs = obj.params.zsweep.numzs;
            
            
            sigmin = obj.params.zsweep.sigmin;
            sigmax = obj.params.zsweep.sigmax;
            numsigs = obj.params.zsweep.numsigs;
            
            %zmin_an = obj.params.zsweep.zmin_an;
            %zmax_an = obj.params.zsweep.zmax_an;
            %numzs_an = obj.params.zsweep.numzs_an;
            numsteps = length(0:obj.params.dt:obj.params.Tmax);
            
            %% Arrays
            zarray = linspace(zmin,zmax,numzs);
            %larray = -zarray./log(f);
            %farray = exp(-zarray./lambda);
            %zarray_an = linspace(zmin_an,zmax_an,numzs_an);
            %larray_an = -zarray_an./log(f);
            %farray_an = exp(-zarray_an./lambda);
            %         momentsL = zeros(numsteps,M,numel(larray));
            %         momentsF = zeros(numsteps,M,numel(larray));
            %statMomL = zeros(numzs_an,M);
            %statMomF = zeros(numzs_an,M);
            
            sigarray = sqrt(2*zarray);
            
            obj.momentSweep.momentsS = zeros(numsteps,M,numsigs);
            %% Compute
            % sweep sigma
            for j = 1:numsigs
                disp(['Beginning sigma # ' num2str(j) ' of ' num2str(obj.params.zsweep.numsigs)])
                
                obj.momentSweep.momentsS(:,:,j) = lesMoments(mu,Kparams,sigarray(j),numtrials,Tmax,dt,lextinct,M);
                
            end
            
            %% Plot
            if lplot
                
                obj.plotMomentSweepLES();
                
            end
            
       end
        
       function obj = plotMomentSweepLES(obj)
            statMomS = zeros(obj.params.zsweep.numsigs_an,obj.params.M);
            sigarray_an = linspace(obj.params.zsweep.sigmin_an,obj.params.zsweep.sigmax_an,...
                obj.params.zsweep.numsigs_an);
            
            
            meanEndS = reshape(obj.momentSweep.momentsS(end,1,:),1,obj.params.zsweep.numzs);
            varEndS = reshape(obj.momentSweep.momentsS(end,2,:) - obj.momentSweep.momentsS(end,1,:).^2,1,obj.params.zsweep.numsigs);
            
           
            for jj = 1:obj.params.zsweep.numsigs_an
                statMomS(jj,:) = lesExactStationaryMoments(obj.params.mu,obj.params.Kparams(1),...
                    sigarray_an(jj),obj.params.M);
            end
            
            
            meanEndS_an = statMomS(:,1);
            
            varEndS_an = statMomS(:,2) - meanEndS_an.^2;
            
            
            figure; hold on;
            plot(meanEndS./obj.params.Kparams(1),varEndS./obj.params.Kparams(1)./obj.params.Kparams(1),'ks','markersize',24,'markerfacecolor',obj.params.bcolor);
            plot(meanEndS_an./obj.params.Kparams(1),varEndS_an./obj.params.Kparams(1)./obj.params.Kparams(1),'-','linewidth',4,'color',[.3 .3 .3])
            
            set(gca,'fontsize',24,'linewidth',4)
            xlabel('E[X]/K','fontsize',24);
            ylabel('Var[X]/K^2','fontsize',24);
            axis([0 1 0 .3])
            
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
            method = obj.params.distrafo.method;
            
            
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
                
                [pophist,~,mendvec(s),varvec(s)] = computeLRCendHist(rvec(s),Kvec(s),lvec(s),fvec(s),dtvec(s),Tmax,numtrials,lextinct,method,bins,false);
                
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
        
        function obj = computeMTElrc(obj)
            
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
            
            obj.mte.lrcMTE = lrcMTE(mu,K,lambda,f,dt,numtrials,lplot);
            
            
            
            
        end
        
        function obj = compareMTElrcles_sigma(obj)
            
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
            
            if numel(fieldnames(obj.params.mte))==0
                obj.params = obj.params.initMTE();
            end
            
            lambdavec = linspace(lambda,obj.params.mte.lmultmax*lambda,obj.params.mte.numls);
            
            
            [obj.mte.lrcMTE_sig,obj.mte.lesMTE_sig] = compareMTE_RC_ES(mu,K,dt,lambdavec,f,numtrials);
            
        end

        function obj = compareMTElrcles_K(obj)
            
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
            
            if numel(fieldnames(obj.params.mte))==0
                obj.params = obj.params.initMTE();
            end
            
            Kvec = logspace(log10(K),obj.params.mte.Kmultmax*log10(K),obj.params.mte.numKs);
            
            
            [obj.mte.lrcMTE_K,obj.mte.lesMTE_K] = compareMTE_RC_ES_K(mu,Kvec,dt,lambda,f,numtrials);
            
        end
        
        function save(obj)
            thisobj = obj;
          
            if(~isdir(obj.savedir))
                mkdir(obj.savedir);
            end
            
            % in case savedir was copied with a '\' at the end, must be
            % consistent
            if obj.savedir(end) == filesep
                obj.savedir = obj.savedir(1:end-1);
            end
            
            save(strcat(obj.savedir,filesep,obj.savename),'thisobj');
            
         end            

    end
    
    methods(Static)
        
        function [c,mat] = m2c(m)
            nmoms = numel(m);
            mat = zeros(nmoms);
            mat(:,1) = m;
            
            % make c mat
            for n = 2:nmoms
                mat(:,n) = circshift(m,n-1);
                mat(n-1,n) = 1;
                
                if n > 2
                    
                    mat(1:n-2,n) = 0;
                    for p = n:nmoms
                        mat(p,n) = nchoosek(p-1,n-2)*mat(p,n);
                    end
                end
                
            end
            
            c = (-1)^(nmoms-1)*det(mat);
            
                        
        end
        
    end
    
end
  