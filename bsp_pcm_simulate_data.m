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
    numSim = 100; %number of subjects to simulate
    numVox = 50; %number of voxels to simulate
    signal = 1; %integer signal level (default = 0.1)
    noise = 2.0; %integer noise level (default = 1)
    suffix = 'regress2.0c'; %suffix for rawts_signal filename
    vararginoptions(varargin,{'num_subj','numSim','signal','noise','suffix'});
    
    %Usage: bsp_pcm_simulate_data('simulate_data','num_subj',1,'numSim',100,'numVox',50,'signal',1,'noise',0,'suffix','highSNR');
    

%     load(fullfile(fullfile(baseDir,simDir,'test_GLM_physio_task_instruc_model_tikhonov.mat')));
    load(fullfile(fullfile(baseDir,simDir,'test_GLM_physio_all_tikhonov_cerebellum.mat')));

    % exp = 2.0; theta = 0.7
    % exp = 1.75; theta = 0.56
    % exp = 1.5; theta = 0.41
    % exp = 1.25; theta = 0.22
    % exp = 1.0; theta = 1e-6
    % exp = 0.75; theta = -0.29
    % exp = 0.5; theta = -0.7
    % exp = 0.25; theta = -1.4
    
    thetaSubj = [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 -1e-6 -1e-6 -1e-6 -1e-6 -1e-6 -1e-6 -1e-6 -1e-6 -1e-6 -1e-6];
%     thetaSubj = [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7];
    X = design(1:numRuns:end,:);
    
    signalVector = signal*ones(numSim,1);
    noiseVector = noise*ones(numSim,1);

    for s = 1:length(num_subj);
        featureTemp = reshape(X(s,:),((numTRs-numDummys)*numRuns),[]);
        feature = num2cell(featureTemp,[1,2]);
        [M,Z] = pcm_buildModelFromFeatures(feature,'name','pontine');
        M.numGparams = 20; %manually add numGparams, as pcm_buildModelFromFeatures doesn't
        M.Gc = M.Gc .* thetaSubj;
        M.Gd = M.Gd .* thetaSubj;
    
        Ysim = pcm_makeDataset(M,thetaSubj','design',Z,'numVox',numVox,'numSim',numSim,'signal',signalVector,'noise',noiseVector);

        for k = 1:length(Ysim)
            Y = Ysim{1,k};
            resMS = ones(1,size(Y,2));

            filename = fullfile(fullfile(baseDir,simDir,'data','S98',sprintf('rawts_simulate_%s_%04d.mat',suffix,k)));
            save(filename,'Y','resMS','-v7.3');
        end
    end
end
