function varargout=bsp_pcm_simulate_data(what,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseDir         ='/srv/diedrichsen/data/Cerebellum/Pontine7T';
regDir          ='/RegionOfInterest';
simDir          = '/simulations';
resDir          = '/results';
subj_name       = {'S98','S97','S96','S95','S01','S03','S04','S07'};
numDummys       = 3;
numTRs          = 328;
numRuns         = 16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(what)

case 'simulate_data'
    num_subj = [1:length(subj_name)];
    numSim = 1; %number of voxels to simulate
    signal = 1; %integer signal level
    noise = 0; %integer noise level
    suffix = 'allsignal'; %suffix for rawts_signal filename
    vararginoptions(varargin,{'num_subj','numSim','signal','noise','suffix'});
    
    %Usage: bsp_pcm_simulate_data('simulate_data','num_subj',1,'numSim',50,'signal',1,'noise',0,'suffix','highSNR');
    

    load(fullfile(fullfile(baseDir,resDir,'/test_GLM_physio_task_instruc_model_tikhonov.mat')));
%     condition = [1 1 1 1 1 1 1 1 1 2]; %condition vector and model G parameters need to match
%     condition = [1];
    condition = ones(((numTRs-numDummys)*numRuns),1);

%     thetaSubj = theta(1:numRuns:end,1:2);  %only grab thetas for task and instruction, not constant
    thetaSubj = [10];
    design = X(1:numRuns:end,:);
    
    signalVector = signal*ones(numSim,1);
    noiseVector = noise*ones(numSim,1);

    for s = 1:length(num_subj);
        featureTemp = reshape(X(s,:),((numTRs-numDummys)*numRuns),[]);
%         feature = sum(feature,2);
        feature = num2cell(featureTemp,[1,2]);
        [M,Z] = pcm_buildModelFromFeatures(feature,'name','pontine','style','rsa_style','type','component');
%         M.Gd = [1 1];   %manually add Gd, Gc, and numGparams, as pcm_buildModelFromFeatures doesn't
%         M.Gc = [1 0; 0 1];
%         M.numGparams = 2;
%         M.Gd = [1];
%         M.Gc = [1];
        M.numGparams = 1;
        
        numVox = size(feature{1,1},1);
    
        Ysim = pcm_makeDataset(M,thetaSubj','design',condition','numVox',numVox,'numSim',numSim,'signal',signalVector,'noise',noiseVector);

        Y = cell2mat(Ysim)';
        Y = sum(Y,2); %sum along condition dimension to produce single time series per simulated voxel
        Y = reshape(Y,[],numSim);
        resMS = ones(1,size(Y,1));

        filename = fullfile(fullfile(baseDir,regDir,'data',subj_name{s},sprintf('rawts_simulate_%s.mat',suffix)));
        save(filename,'Y','resMS','-v7.3');
    end
end
