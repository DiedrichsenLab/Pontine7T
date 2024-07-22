function XX=bsp_get_feature(what,sn,SPM,INFO,varargin)

% Function that gets the regressors / features for a specific GLM 
if isdir('/Volumes/diedrichsen_data$/data')
    workdir='/Volumes/diedrichsen_data$/data';
elseif isdir('/srv/diedrichsen/data')
    workdir='/srv/diedrichsen/data';
else
    fprintf('Workdir not found. Mount or connect to server and try again.');
end

baseDir=(sprintf('%s/Cerebellum/Pontine7T',workdir));
pinfo = dload(fullfile(baseDir,'participants_new_format.tsv')); 
subj_name = pinfo.participant_id;

separate = 0;  % Return seperate regressor for each run (or concatinated)
sscale = 0;   % Seperate z-normalization for each run 
zscale = 1;   % Common z-normalization across all runs 

vararginoptions(varargin,{'separate','sscale','zscale'}); 

nRuns = length(SPM.nscan);
switch (what)
    case 'Tasks'
        for rn = 1:nRuns
            indx = SPM.Sess(rn).col;
            ii=INFO.inst(indx)==0;
            X{rn} = SPM.xX.X(SPM.Sess(rn).row,indx(ii));
        end
    case 'Instruct'
        for rn = 1:nRuns
            indx = SPM.Sess(rn).col;
            ii=INFO.inst(indx)==1;
            X{rn} = SPM.xX.X(SPM.Sess(rn).row,indx(ii));
        end
    case 'InstructC'
        for rn = 1:nRuns
            indx = SPM.Sess(rn).col;
            ii=INFO.inst(indx)==1;
            X{rn} = sum(SPM.xX.X(SPM.Sess(rn).row,indx(ii)),2);
        end
    case 'Null'         % Null model - row/column permuted Tasks+InstructC
        numTRs = 325;
        for rn = 1:nRuns
            A = load(fullfile(baseDir,'simulations','design.mat'));
            Atemp = A.X(1:numTRs,:);
            Atemp = Atemp(randperm(size(Atemp,1)),:); % randomly permute rows
            Atemp = Atemp(:,randperm(size(Atemp,2))); % randomly permute columns
            X{rn} = Atemp;
        end
    case 'Hpass'        % High pass filter
        for rn = 1:nRuns
            k.HParam = 128;
            k.RT     = SPM.xY.RT;
            k.row  = SPM.Sess(rn).row;
            HPF = spm_filter(k);
            X{rn} = HPF.X0;
        end
    case 'Mov'          % Realignment parameters
        % load the motion files of every run
        for rn = 1:nRuns
            X{rn} = load(fullfile(baseDir,'imaging_data',subj_name{sn},sprintf('rp_run_%02d.txt',rn)));
        end
    case 'MovPCA'       % 2 PCs of realigment paramenters
        % load the motion files of every run
        MOV  = [];
        for rn = 1:nRuns
            mov = load(fullfile(baseDir,'imaging_data',subj_name{sn},sprintf('rp_run_%02d.txt',rn)));
            MOV = [MOV;mov];
        end
        % Get the principal components
        [~,score] = pca(MOV);
        for rn = 1:nRuns
            X{rn} = score(SPM.Sess(rn).row,1:2);
        end
    case 'Retro_HR'    % Retroicor of cardio and resp 18 comp
        % Load the physio regressors from TAPAS
        for rn = 1:nRuns
            A = load(fullfile(baseDir,'physio',subj_name{sn},sprintf('run%02d',rn),'reg_retro_hr.txt'));
            X{rn} = A(:,1:4); % Two fundamentals
        end
    case 'Retro_RESP'    % Retroicor of cardio and resp 18 comp
        % Load the physio regressors from TAPAS
        for rn = 1:nRuns
            A = load(fullfile(baseDir,'physio',subj_name{sn},sprintf('run%02d',rn),'reg_retro_resp.txt'));
            X{rn} = A(:,1:4); % Two fundamentals
        end
    case 'HR'    % Raw HEart rate 
        for rn = 1:nRuns
            A = load(fullfile(baseDir,'physio',subj_name{sn},sprintf('run%02d',rn),'reg_hr.txt'));
            X{rn} = A(:,2);
        end
    case 'HRconv'    % HEart rate convolved and z-normalized
        % Load the physio regressors from TAPAS
        for rn = 1:nRuns
            A = load(fullfile(baseDir,'physio',subj_name{sn},sprintf('run%02d',rn),'reg_hr.txt'));
            X{rn} = A(:,1);
        end
    case 'RV'    % Retroicor of cardio and resp 18 comp
        % Load the physio regressors from TAPAS
        for rn = 1:nRuns
            A = load(fullfile(baseDir,'physio',subj_name{sn},sprintf('run%02d',rn),'reg_rvt.txt'));
            X{rn} = A(:,1);
        end    
    case 'aCompCor'    % anatomical compcor estimated from fmriprep
        % Load the acompcor regressors from fmriprep confounds file
        for rn = 1:nRuns
            A = load(fullfile(baseDir,'confounds',subj_name{sn},sprintf('run%02d',rn),'acompcor.txt'));
            X{rn} = A(4:328,1:2); % grab first two components
        end
    case 'CSF'          % Mean signal of the CSF around brainstem
        % Get the CSF data
        csf = load(fullfile(baseDir,'RegionOfInterest','data',subj_name{sn},'rawts_csf.mat'));
        % mean CSF signal
        for rn = 1:nRuns
            mcsf = mean(csf.Y(SPM.Sess(rn).row,:),2);
            X{rn} = bsxfun(@rdivide,mcsf,sqrt(sum(mcsf.^2)));
        end
    case 'CSFPCAindiv'       % 2 Pcs of CSF
        % Get the CSF data
        csf = load(fullfile(baseDir,'RegionOfInterest','data',subj_name{sn},'rawts_csf.mat'));
        % Compute the PCs
        for rn = 1:nRuns
            runcsf = csf.Y(SPM.Sess(rn).row,:);
            % Get the principal components
            [~,score] = pca(runcsf);
            X{rn} = score(:,1:2);
        end
    case 'CSFPCAall'   % 2 PCs of CSF computed over 4 runs
        % Get the CSF mask
        csf = load(fullfile(baseDir,'RegionOfInterest','data',subj_name{sn},'rawts_csf.mat'));
        [~,score] = pca(csf.Y);
        % Include the first 2 PCs of CSF in design
        for rn = 1:nRuns
            X{rn} = score(SPM.Sess(rn).row,1:2);
        end
        
    case 'CSFpons'          % Mean signal of the CSF around brainstem
        % Get the CSF data
        csf = load(fullfile(baseDir,'RegionOfInterest','data',subj_name{sn},'rawts_pons.mat'));
        % mean CSF signal
        for rn = 1:nRuns
            mcsf = mean(csf.Y(SPM.Sess(rn).row,:),2);
            X{rn} = bsxfun(@rdivide,mcsf,sqrt(sum(mcsf.^2)));
        end
    case 'CSFponsPCAindiv'       % 2 Pcs of CSF
        % Get the CSF data
        csf = load(fullfile(baseDir,'RegionOfInterest','data',subj_name{sn},'rawts_pons.mat'));
        % Compute the PCs
        for rn = 1:nRuns
            runcsf = csf.Y(SPM.Sess(rn).row,:);
            % Get the principal components
            [~,score] = pca(runcsf);
            X{rn} = score(:,1:2);
        end
    case 'CSFponsPCAall'   % 2 PCs of CSF computed over 4 runs
        % Get the CSF mask
        csf = load(fullfile(baseDir,'RegionOfInterest','data',subj_name{sn},'rawts_pons.mat'));
        [~,score] = pca(csf.Y);
        % Include the first 2 PCs of CSF in design
        for rn = 1:nRuns
            X{rn} = score(SPM.Sess(rn).row,1:2);
        end
    case 'CSFmedulla'          % Mean signal of the CSF around brainstem
        % Get the CSF data
        csf = load(fullfile(baseDir,'RegionOfInterest','data',subj_name{sn},'rawts_medulla.mat'));
        % mean CSF signal
        for rn = 1:nRuns
            mcsf = mean(csf.Y(SPM.Sess(rn).row,:),2);
            X{rn} = bsxfun(@rdivide,mcsf,sqrt(sum(mcsf.^2)));
        end
    case 'CSFmedullaPCAindiv'       % 2 Pcs of CSF
        % Get the CSF data
        csf = load(fullfile(baseDir,'RegionOfInterest','data',subj_name{sn},'rawts_medulla.mat'));
        % Compute the PCs
        for rn = 1:nRuns
            runcsf = csf.Y(SPM.Sess(rn).row,:);
            % Get the principal components
            [~,score] = pca(runcsf);
            X{rn} = score(:,1:2);
        end
    case 'CSFmedullaPCAall'   % 2 PCs of CSF computed over 4 runs
        % Get the CSF mask
        csf = load(fullfile(baseDir,'RegionOfInterest','data',subj_name{sn},'rawts_medulla.mat'));
        [~,score] = pca(csf.Y);
        % Include the first 2 PCs of CSF in design
        for rn = 1:nRuns
            X{rn} = score(SPM.Sess(rn).row,1:2);
        end
    case 'CSFpostdrain'          % Mean signal of the CSF around brainstem
        % Get the CSF data
        csf = load(fullfile(baseDir,'RegionOfInterest','data',subj_name{sn},'rawts_postdrain.mat'));
        % mean CSF signal
        for rn = 1:nRuns
            mcsf = mean(csf.Y(SPM.Sess(rn).row,:),2);
            X{rn} = bsxfun(@rdivide,mcsf,sqrt(sum(mcsf.^2)));
        end
    case 'CSFpostdrainPCAindiv'       % 2 Pcs of CSF
        % Get the CSF data
        csf = load(fullfile(baseDir,'RegionOfInterest','data',subj_name{sn},'rawts_postdrain.mat'));
        % Compute the PCs
        for rn = 1:nRuns
            runcsf = csf.Y(SPM.Sess(rn).row,:);
            % Get the principal components
            [~,score] = pca(runcsf);
            X{rn} = score(:,1:2);
        end
    case 'CSFpostdrainPCAall'   % 2 PCs of CSF computed over 4 runs
        % Get the CSF mask
        csf = load(fullfile(baseDir,'RegionOfInterest','data',subj_name{sn},'rawts_postdrain.mat'));
        [~,score] = pca(csf.Y);
        % Include the first 2 PCs of CSF in design
        for rn = 1:nRuns
            X{rn} = score(SPM.Sess(rn).row,1:2);
        end
end
if (zscale) % subtract mean from each regressor - seperate for each regressor
    for rn = 1:nRuns
        X{rn}=bsxfun(@minus,X{rn},mean(X{rn}));
    end
end
if (sscale)
    for rn = 1:nRuns
        VAR(rn,:) = sum(X{rn}.*X{rn})/size(X{rn},1);
    end
    for rn = 1:nRuns
        X{rn}=bsxfun(@rdivide,X{rn},sqrt(mean(VAR,1)));
    end
end
N = size(SPM.xX.X,1);
K = size(X{1},2);
if (separate)
    XX = zeros(N,K*nRuns);
    for rn = 1:nRuns
        XX(SPM.Sess(rn).row,[1:K]+(rn-1)*K) = X{rn};
    end
else
    XX = zeros(N,K);
    for rn = 1:nRuns
        num_rows = length(SPM.Sess(rn).row);
        XX(SPM.Sess(rn).row,:) = X{rn}(1:num_rows, :);
    end
end
end