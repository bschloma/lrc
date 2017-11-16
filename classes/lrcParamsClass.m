% Program:  lrcParamsClass.m
%
% Summary:  Class for keeping track of lrc and les model parameters.
%           Minimal methods.  
%
%           Properties:
%                 %%main params, used with lrc.m and les.m%%
%                 savename - string, filename;
%                 savedir - string, fullpath to save directory;
%                 mu - growth rate;
%                 Kparams - 1x2 array, mean and stddev of carrying cap;
%                 Tmax - Max simulation time;
%                 dt - timestep;
%                 numtrials - number of trials in ensemble;
%                 f - collapse fraction;
%                 l - average collaspe rate;
%                 lextinct - logical for extinction;
%                 lplot - logical for plotting;
%                 sigma - les model volatility;
%                 poiscolor - 1x3 array for lrc color;
%                 bcolor - 1x3 array for les color;
% 
%                 %% structs for additional computation params%%
%                 zsweep = struct();
%                 distrafo = struct();
%                 mte = struct();
% 
% Usage:    params = lrcParamsClass(savename,savedir)
%           savename and savedir optional.
%
% Author:   Brandon Schlomann
%
% Date:     4/14/17 - First written.
% 



classdef lrcParamsClass
    
    properties
        
        savename = '';
        savedir = '';
        mu = 0.;
        Kparams = [0. 0.];
        Tmax = 0.;
        dt = 0.;
        numtrials = 0;
        f = 0.;
        lambda = 0.;
        lextinct = 0;
        lplot = 0.;
        sigma = 0;
        poiscolor = [0. 0. 0.];
        bcolor = [0. 0. 0.];
        M = 0;
        method = '';
        x0 = 0.;
      

        zsweep = struct();
        distrafo = struct();
        mte = struct();
    end
    
    methods
        
        function obj = lrcParamsClass(varargin)
            
            if nargin == 0
                c= clock;
                runname = '';
                for cc = 1:5   %leave off seconds
                    runname = [runname num2str(c(cc)) '_'];
                end
                
                runname = runname(1:(end-1));    %remove trailing _
                
                obj.savename= ['lrcParams_' runname];
                
                obj.savedir = pwd;
            elseif nargin == 1
                obj.savename = varargin{1};
                obj.savedir = pwd;
            elseif nargin == 2
                obj.savename = varargin{1};
                obj.savedir = varargin{2};
            else
                disp('lrcParamsClass takes at most 2 inputs')
            end
            
            
          
            obj.mu = 1;
            obj.Kparams = [1e4 0];
            obj.f = 1e-2;
            obj.lambda = 1e-1;
            obj.Tmax = 50;
            obj.dt = .01;
            obj.numtrials = 1e3;
            obj.lextinct = false;
            obj.lplot = true;
            obj.sigma = .8;
            obj.M = 8;
            obj.method = 'pdmp';
            obj.x0 = 10;
            
            %obj.zsweep = struct();
            %obj.distrafo = struct();
            %obj.mte = struct();
           obj.poiscolor =[.4,.1,.8];
            obj.bcolor = [.1,.7,.2];
%             obj.zmin = 1e-3;
%             obj.zmax = .8;
%             obj.numzs = 6;
%            
%             obj.M = 8;
%             obj.zmin_an = 1e-3;
%             obj.zmax_an = 1.1;
%             obj.numzs_an = 300;
%             obj.poiscolor =[.4,.1,.8];
%             obj.bcolor = [.1,.7,.2];
%             
%             obj.sigmin = .0047;
%             obj.sigmax = 1.2649;
%             obj.numsigs = 6;
%             obj.sigmin_an = .0047;
%             obj.sigmax_an = 1.2649;
%             obj.numsigs_an = 300;
%             
%             obj.numbins = 20;
%             
%             obj.numalphs = 5;
%             obj.lmultmax = 75;
%             
%             obj.ymax = .3;
%             
        end
        
        function obj = initZSweep(obj)
            obj.zsweep.zmin = 1e-3;
            obj.zsweep.zmax = .8;
            obj.zsweep.numzs = 6;
           
            obj.zsweep.zmin_an = 1e-3;
            obj.zsweep.zmax_an = 1.1;
            obj.zsweep.numzs_an = 300;
            
            obj.zsweep.sigmin = .0047;
            obj.zsweep.sigmax = 1.5;%:1.2649;
            obj.zsweep.numsigs = 6;
            obj.zsweep.sigmin_an = .0047;
            obj.zsweep.sigmax_an = 1.2649;
            obj.zsweep.numsigs_an = 300;
        end
        
        function obj = initDisTrafo(obj)
            
            obj.distrafo.numbins = 20;
            
            obj.distrafo.numalphs = 5;
            obj.distrafo.lmultmax = 75;
            
            obj.distrafo.ymax = 3;
            obj.distrafo.method = 'euler';
        end
        
         function obj = initMTE(obj)
                        
            obj.mte.numls = 5;
            obj.mte.lmultmax = 75;
            obj.mte.numKs = 5;
            obj.mte.Kmultmax = 2.5;
            obj.mte.numx0s = 5;
            obj.mte.x0multmax = 100;
            
            obj.mte.maptype = 1;
            
         end
        
         function obj = initMTEtrafo(obj)
             obj.mte.mtetrafo = struct();
             
             obj.mte.mtetrafo.numalphs = 5;
             obj.mte.mtetrafo.lmultmax = 75;
             obj.mte.mtetrafo.scalevec = [];
                          
         end
             
%         

%         function obj = mapSigma(obj)
%             
%             obj.distrafo.sigmin = sqrt(2*params.zmin);
%             obj.distrafo.sigmax = sqrt(2*params.zmax);
%             obj.distrafo.numsigs = params.numzs;
%             
%             obj.distrafo.sigmin_an = sqrt(2*params.zmin_an);
%             obj.distrafo.sigmax_an = sqrt(2*params.zmax_an);
%             obj.distrafo.numsigs_an = params.numzs_an;

            
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
    
end