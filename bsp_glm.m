function varargout=bsp_glm(what,varargin)
% GLM analysis for the pontine7T project, starting from preprocessed data
% (see bsp_imana).

%

numDummys = 3;                                                              % per run
numTRs    = 328;                                                            % per run (includes dummies)
%========================================================================================================================
% PATH DEFINITIONS
global baseDir;
global subj_name;

baseDir         ='/srv/diedrichsen/data/Cerebellum/Pontine7T';
if ~exist(baseDir,'dir')
    baseDir         ='/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T';
end
imagingDir      ='/imaging_data';
imagingDirRaw   ='/imaging_data_raw';
anatomicalDir   ='/anatomicals';
suitDir         ='/suit';
regDir          ='/RegionOfInterest';
%========================================================================================================================
% PRE-PROCESSING
subj_name = {'S99','S98','S97','S96'};
%========================================================================================================================
% GLM INFO
funcRunNum  = [1,16];  % first and last behavioural run numbers
run         = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
runB        = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];  % Behavioural labelling of runs
sess        = [1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2];                  % session number

%========================================================================================================================
switch(what)
    case 'GLM:do'
        bsp_glm('GLM:glm1',[4],[1:16]);
        bsp_glm('GLM:glm2',[4],[1:16]);
        bsp_glm('GLM:estimate',[4],1);
        bsp_glm('GLM:estimate',[4],2);
        bsp_glm('GLM:Fcontrast','sn', [4], 'glm', 1, 'type', 'task');
        bsp_glm('GLM:Fcontrast','sn', [4], 'glm', 2, 'type', 'task');
    
    case 'GLM:makeMask' % Changed to include CSF
        sn=varargin{1}; % subjNum
        tissues = [1:3];
        
        P{1} = fullfile(fullfile(baseDir,imagingDir,subj_name{sn},'rmeanrun_01.nii'));
        for i=1:length(tissues)
            P{i+1}=fullfile(baseDir,anatomicalDir,subj_name{sn},sprintf('c%danatomical.nii',tissues(i)));
        end
        out =  fullfile(fullfile(baseDir,imagingDir,subj_name{sn},'brain_mask.nii'));
        spm_imcalc(char(P),out,'i1>800 & (i2+i3+i4)>0.7');
        
        out =  fullfile(fullfile(baseDir,imagingDir,subj_name{sn},'gray_mask.nii'));
        spm_imcalc(char(P),out,'i1>800 & i2>0.4');
    
    case 'GLM:glm1'                   % FAST glm w/out hpf one regressor per task and per instruction
        % GLM with FAST and no high pass filtering
        % 'spm_get_defaults' code modified to allow for -v7.3 switch (to save>2MB FAST GLM struct)
        % EXAMPLE: bsp_imana('GLM:glm1',[1:XX],[1:XX])
        sn=varargin{1};
        runs=varargin{2}; % [1:16]
        
        announceTime=5;  % Instruction duration
        glm=1;
        subjs=length(sn);
        
        % load in task information
        C=dload(fullfile(baseDir,'pontine_taskConds_GLM.tsv'));
        Cc=getrow(C,C.StudyNum==1 & C.condNum==1); % Only get the first condition for each tasks
        nTask      = max(Cc.taskNum);                                 % how many tasks there are?
        
        for s=1:subjs,
            T=[];
            A = dload(fullfile(baseDir,'data','fmri',subj_name{sn(s)},['fmri_',subj_name{sn(s)},'.tsv'])); % get scanning timing and order
            A = getrow(A,A.run_num>=funcRunNum(1) & A.run_num<=funcRunNum(2)); % get only the runs we need (remove any test or Behav training)
            
            glmSubjDir =[baseDir, sprintf('/GLM_firstlevel_%d/',glm),subj_name{sn(s)}];dircheck(glmSubjDir); % Create GLM folder
            
            % Fill up struct for glm
            J.dir = {glmSubjDir};
            J.timing.units = 'secs';
            J.timing.RT = 1.0;
            J.timing.fmri_t = 16;
            J.timing.fmri_t0 = 1;
            
            for r=1:numel(runs) % loop through runs
                P=getrow(A,A.run_num==runB(r));
                for i=1:(numTRs-numDummys) % get the filenames of the nifti volumes for each run
                    N{i} = [fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('run_%2.2d.nii,%d',runs(r),i))];
                end;
                J.sess(r).scans= N; % number of scans in run
                
                % loop through tasks
                itt = 1;
                for it = 1:nTask
                    % The order of tasks are different for each run, to
                    % have a common order for the tasks, I will be reading
                    % from the Cc file for all the runs and subjects
                    ST = find(strcmp(lower(P.task_name),lower(Cc.taskNames{it})));
                    if (isempty(ST) || length(ST)>1)
                        keyboard;
                    end;
                    for taskType = 1:2 % there are two taskTypes: instruction; not instructions
                        if taskType == 1 % instructions
                            % get the isntruction onset
                            instruct_onset = P.real_start_time(ST)- J.timing.RT*numDummys; %% get the instruction start time for the first task
                            
                            % filling in the fields for SPM_info.mat
                            S.task      = it;              % task Number of the task after
                            S.taskName   = {'Instruct'};   % task name (TN)
                            S.inst      = 1;              % is it instruction (1) or not (0)?
                            S.taskOrder = ST;             % instOrder Number of task
                            S.time      = instruct_onset; % instruction onset time
                            S.sess      = sess(r);
                            S.run       = runB(r);
                            % Determine taskName_after and taskName_before
                            % this instruction
                            S.taskName_after  = P.task_name(ST);
                            if ST > 1
                                S.taskName_before = P.task_name(ST-1);
                            elseif ST == 1
                                S.taskName_before = {'NaN'};
                            end
                            T  = addstruct(T, S);
                            
                            % filling in the fields for SPM.mat
                            J.sess(r).cond(itt).name     = 'Instruct';
                            J.sess(r).cond(itt).onset    = instruct_onset; % correct start time for numDummys and announcetime included (not for instruct)
                            J.sess(r).cond(itt).duration = 5;              % instructions last for 5 sec
                            J.sess(r).cond(itt).tmod     = 0;
                            J.sess(r).cond(itt).orth     = 0;
                            J.sess(r).cond(itt).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                            
                            itt = itt + 1;
                        elseif taskType == 2 % not instructions
                            % get the task onset (instruction onset + announceTime)
                            onset = P.real_start_time(ST) - J.timing.RT*numDummys +announceTime;
                            
                            % filling in the fields for SPM_info.mat
                            S.task      = it;
                            S.taskName  = {Cc.taskNames{it}};
                            S.inst      = 0;
                            S.taskOrder = ST;
                            S.time      = onset;
                            S.taskName_after  = {'none'}; % taskName before and after are only defined for instructions
                            S.taskName_before = {'none'};
                            S.sess      = sess(r);
                            S.run       = runB(r);
                            
                            T  = addstruct(T, S);
                            
                            % filling in the fields for SPM.mat
                            J.sess(r).cond(itt).name     = Cc.taskNames{it};
                            J.sess(r).cond(itt).onset    = onset;
                            J.sess(r).cond(itt).duration = 30;             % each task lasts for 30 sec
                            J.sess(r).cond(itt).tmod     = 0;
                            J.sess(r).cond(itt).orth     = 0;
                            J.sess(r).cond(itt).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                            
                            itt = itt + 1;
                        end % if it's instructions or not?
                    end % taskType (looping over instruction and non-instructions)
                end % it (tasks)
                J.sess(r).multi = {''};
                J.sess(r).regress = struct('name', {}, 'val', {});
                J.sess(r).multi_reg = {''};
                J.sess(r).hpf = inf;                                        % set to 'inf' if using J.cvi = 'FAST'. SPM HPF not applied
            end
            J.fact = struct('name', {}, 'levels', {});
            J.bases.hrf.derivs = [0 0];
            J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
            J.volt = 1;
            J.global = 'None';
            J.mask = {fullfile(baseDir,imagingDir,subj_name{sn(s)},'brain_mask.nii')};
            J.mthresh = 0.01;
            J.cvi_mask = {fullfile(baseDir,imagingDir,subj_name{sn(s)},'gray_mask.nii')};
            J.cvi =  'fast';
            
            spm_rwls_run_fmri_spec(J);
            
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            fprintf('glm_%d has been saved for %s \n',glm, subj_name{sn(s)});
        end
        
    case 'GLM:glm2'                   % FAST glm w/out hpf one regressor per task including instruction
        % GLM with FAST and no high pass filtering
        % 'spm_get_defaults' code modified to allow for -v7.3 switch (to save>2MB FAST GLM struct)
        % EXAMPLE: bsp_imana('GLM:glm1',[1:XX],[1:XX])
        sn=varargin{1};
        runs=varargin{2}; % [1:16]
        
        announceTime=5;  % Instruction duration
        glm=2;
        subjs=length(sn);
        
        % load in task information
        C=dload(fullfile(baseDir,'pontine_taskConds_GLM.tsv'));
        Cc=getrow(C,C.StudyNum==1 & C.condNum==1); % Only get the first condition for each tasks
        nTask      = max(Cc.taskNum);                                 % how many tasks there are?
        
        for s=1:subjs,
            T=[];
            A = dload(fullfile(baseDir,'data','fmri',subj_name{sn(s)},['fmri_',subj_name{sn(s)},'.tsv'])); % get scanning timing and order
            A = getrow(A,A.run_num>=funcRunNum(1) & A.run_num<=funcRunNum(2)); % get only the runs we need (remove any test or Behav training)
            
            glmSubjDir =[baseDir, sprintf('/GLM_firstlevel_%d/',glm),subj_name{sn(s)}];dircheck(glmSubjDir); % Create GLM folder
            
            % Fill up struct for glm
            J.dir = {glmSubjDir};
            J.timing.units = 'secs';
            J.timing.RT = 1.0;
            J.timing.fmri_t = 16;
            J.timing.fmri_t0 = 1;
            
            for r=1:numel(runs) % loop through runs
                P=getrow(A,A.run_num==runB(r));
                for i=1:(numTRs-numDummys) % get the filenames of the nifti volumes for each run
                    N{i} = [fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('run_%2.2d.nii,%d',runs(r),i))];
                end;
                J.sess(r).scans= N; % number of scans in run
                
                instruct_onset = P.real_start_time - J.timing.RT*numDummys;
                % Add overall Instruction regressor
                S.task      = 0;              % task Number of the task after
                S.taskName   = {'Instruct'};   % task name (TN)
                S.inst      = 1;              % is it instruction (1) or not (0)?
                S.taskOrder = 0;             % instOrder Number of task
                S.time      = instruct_onset(1); % instruction onset time
                S.sess      = sess(r);
                S.run       = runB(r);
                T  = addstruct(T, S);
                
                % filling in the fields for SPM.mat
                J.sess(r).cond(1).name     = 'Instruct';
                J.sess(r).cond(1).onset    = instruct_onset; % correct start time for numDummys and announcetime included (not for instruct)
                J.sess(r).cond(1).duration = 5;              % instructions last for 5 sec
                J.sess(r).cond(1).tmod     = 0;
                J.sess(r).cond(1).orth     = 0;
                J.sess(r).cond(1).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                
                % loop through tasks
                for it = 1:nTask
                    % The order of tasks are different for each run, to
                    % have a common order for the tasks, I will be reading
                    % from the Cc file for all the runs and subjects
                    ST = find(strcmp(lower(P.task_name),lower(Cc.taskNames{it})));
                    
                    % get the task onset (instruction onset + announceTime)
                    onset = P.real_start_time(ST) - J.timing.RT*numDummys +announceTime;
                    
                    % filling in the fields for SPM_info.mat
                    S.task      = it;
                    S.taskName  = {Cc.taskNames{it}};
                    S.inst      = 0;
                    S.instOrder = 0;
                    S.time      = onset;
                    T  = addstruct(T, S);
                    
                    % filling in the fields for SPM.mat
                    J.sess(r).cond(it+1).name     = Cc.taskNames{it};
                    J.sess(r).cond(it+1).onset    = onset;
                    J.sess(r).cond(it+1).duration = 30;             % each task lasts for 30 sec
                    J.sess(r).cond(it+1).tmod     = 0;
                    J.sess(r).cond(it+1).orth     = 0;
                    J.sess(r).cond(it+1).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                    
                end % it (tasks)
                J.sess(r).multi = {''};
                J.sess(r).regress = struct('name', {}, 'val', {});
                J.sess(r).multi_reg = {''};
                J.sess(r).hpf = inf;                                        % set to 'inf' if using J.cvi = 'FAST'. SPM HPF not applied
            end
            J.fact = struct('name', {}, 'levels', {});
            J.bases.hrf.derivs = [0 0];
            J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
            J.volt = 1;
            J.global = 'None';
            J.mask = {fullfile(baseDir,imagingDir,subj_name{sn(s)},'brain_mask.nii')};
            J.mthresh = 0.01;
            J.cvi_mask = {fullfile(baseDir,imagingDir,subj_name{sn(s)},'gray_mask.nii')};
            J.cvi =  'fast';
            
            spm_rwls_run_fmri_spec(J);
            
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            fprintf('glm_%d has been saved for %s \n',glm, subj_name{sn(s)});
        end
        
    case 'GLM:estimate'               % Estimate GLM depending on subjNum & glmNum
        % example: bsp_imana('GLM:estimate_glm',1,1)
        sn=varargin{1};
        glm=varargin{2};
        
        subjs=length(sn);
        
        for s=1:subjs,
            glmDir =[baseDir,sprintf('/GLM_firstlevel_%d/',glm),subj_name{sn(s)}];
            load(fullfile(glmDir,'SPM.mat'));
            SPM.swd=glmDir;
            spm_rwls_spm(SPM);
        end
    case 'GLM:contrast'               % Create Contrast images for each task - and instruction 
        %%% Calculating contrast images.
        % 'SPM_light' is created in this step (xVi is removed as it slows
        % down code for FAST GLM).
        % Example1: bsp_imana('GLM:contrast', 'sn', 3, 'glm', 1, 'type', 'task')
        
        sn             = 3;             %% list of subjects
        glm            = 1;             %% The glm number
        
        vararginoptions(varargin, {'sn', 'glm'})
        
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, sprintf('GLM_firstlevel_%d', glm));
        
        for s = sn
            fprintf('******************** calculating contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'))
            
            SPM  = rmfield(SPM,'xCon');
            cd(fullfile(glmDir, subj_name{s}))
            T = load('SPM_info.mat');
            
            % t contrast for tasks
            ucondition = unique(T.task);
            num_tasks = length(ucondition); 
            idx = 1;
            for tt = 1:num_tasks % 0 is "instruct" regressor
                con  = zeros(1,size(SPM.xX.X,2));
                indx = find(T.task == tt & T.inst==0); 
                con(:,indx) = 1/16;
                con(:,T.task ~= tt & T.inst==0) = -1./(16*(num_tasks-1));                
                name = T.taskName{indx(1)}; 
                SPM.xCon(idx) = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);
                idx=idx+1;
            end 
            % Add Instruction regressor - against mean of tasks 
            con  = zeros(1,size(SPM.xX.X,2));
            indx = find(T.inst>0); 
            con(:,indx) = 1/length(indx);
            con(:,T.task >0 & T.inst==0) = -1./(16*num_tasks);                
            name = 'Instruct'; 
            SPM.xCon(idx) = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);
            
            SPM = spm_contrasts(SPM,1:length(SPM.xCon));
            save('SPM.mat', 'SPM','-v7.3');
            
            % rename contrast images and spmT images
            conName = {'con','spmT'};
            for i = 1:length(SPM.xCon)
                for n = 1:numel(conName)
                    oldName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%2.4d.nii',conName{n},i));
                    newName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                    movefile(oldName{i}, newName{i});
                end % conditions (n, conName: con and spmT)
            end % i (contrasts)
        end % sn
    
    case 'GLM:Fcontrast'               % Create Contrast images
        %%% Calculating contrast images.
        % 'SPM_light' is created in this step (xVi is removed as it slows
        % down code for FAST GLM).
        % Example1: bsp_glm('GLM:Fcontrast', 'sn', 3, 'glm', 1, 'type', 'task')
        % Example2: bsp_glm('GLM:Fcontrast', 'sn', 3, 'glm', 1, 'type', 'cond')
        
        sn             = 3;             %% list of subjects
        glm            = 1;             %% The glm number
        type           = 'task';        %% it can be set to task, ....
        
        vararginoptions(varargin, {'sn', 'glm', 'type'})
        
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, sprintf('GLM_firstlevel_%d', glm));
        
        for s = sn
            fprintf('******************** calculating contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'))
            
            if isfield(SPM,'xCon')
                SPM  = rmfield(SPM,'xCon');
            end
            cd(fullfile(glmDir, subj_name{s}))
            T = load('SPM_info.mat');
            
            % F contrast
            name = sprintf('%s', type);
            switch(type)
                case 'task'
                    numTasks = max(T.task);
                    con = zeros(numTasks,size(SPM.xX.X,2));
                    for i=1:numTasks
                        con(i,T.task==i & T.inst==0)=1-1/numTasks;
                        con(i,T.task~=i & T.inst==0)=-1/numTasks;
                    end
            end
            
            SPM.xCon(1) = spm_FcUtil('Set',name, 'F', 'c',con',SPM.xX.xKXs);
            SPM = spm_contrasts(SPM,1:length(SPM.xCon));
            save('SPM.mat', 'SPM','-v7.3');
            SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
            save(fullfile(glmDir, subj_name{s}, 'SPM_light.mat'), 'SPM');
            
        end % sn
    
    case 'GLM:contrast_F_summary'
        % Calculating contrast images for overall F-contrast between
        % tasks / conditions
        % EXAMPLE: bsp_glm('GLM:contrast_F_summary',[1:XX],1)
        sn=varargin{1};
        glm = varargin{2};           %% The glm number
        
        D=[];
        subjs=length(sn);
        
        for s = 1:subjs
            glmSubjDir =[baseDir, sprintf('/GLM_firstlevel_%d/',glm), subj_name{sn(s)}];
            regDir = [baseDir,'/suit/','anatomicals/',subj_name{sn(s)}];
            
            fprintf('******************** doing calc %s ********************\n', subj_name{sn(s)});
            R{1}=region('image',fullfile(regDir,'rc7anatomical.nii'),0.5);  % 'cortical_mask_grey_corr.nii'
            R{2}=region('image',fullfile(regDir,'rc_anatomical_seg1.nii'),0.5);  % 'regions_cerebellum_suit.nii'
            R = region_calcregions(R);
            
            V= spm_vol(fullfile(glmSubjDir,'spmF_0001.nii'));
            data = region_getdata(V,R);
            T.sn = subj_name{sn(s)};
            T.numzero = [sum(data{1}==0), sum(data{2}==0)];
            data{1}=data{1}(data{1}>0);
            data{2}=data{2}(data{2}>0);
            T.avrgF = [mean(data{1}),mean(data{2})];
            T.prcF = [prctile(data{1},95),prctile(data{2},95)];
            D=addstruct(D,T);
        end % sn
        subplot(1,2,1);
        myboxplot([],D.avrgF);
        title('Average F-value');
        set(gca,'Ylim',[0 max(D.avrgF(:))+0.2]);
        drawline(1,'dir','horz')
        
        subplot(1,2,2);
        myboxplot([],D.prcF)
        set(gca,'Ylim',[0 max(D.prcF(:))+1])
        title('95% F-value');
        drawline(finv(0.95,16,1000),'dir','horz')
        varargout={D};
    
    case 'GLM:restats'
        sn=varargin{1};
        glm = varargin{2};           %% The glm number
        glmSubjDir =fullfile(baseDir, sprintf('GLM_firstlevel_%d',glm), subj_name{sn});
        dataSubjDir =fullfile(baseDir, 'imaging_data',subj_name{sn});
        load(fullfile(glmSubjDir,'SPM.mat'));
        start=[0 cumsum(SPM.nscan)]+1;
        
        for i=1:length(SPM.nscan)
            [dir,filename]=spm_fileparts(SPM.xY.VY(start(i)).fname);
            movparam_name{i}=fullfile(dataSubjDir,['rp_' filename(1:end) '.txt']);
        end;
        spm_rwls_resstats(SPM,[],movparam_name);
        keyboard;
        
    case 'TS:getRawTs'                % Get raw timeseries and save them
        % bsp_imana('TS:getRawTs',1,1,'dentate');
        sn=varargin{1}; % subjNum
        glm=varargin{2}; % glmNum
        
        glmDir =fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm));
        
        for s=sn,
            glmDirSubj=fullfile(glmDir, subj_name{s});
            load(fullfile(glmDirSubj,'SPM.mat'));
            
            % load data
            load(fullfile(baseDir,regDir,'data',subj_name{s},sprintf('regions.mat')));
            % SPM=spmj_move_rawdata(SPM,fullfile(baseDir,imagingDir,subj_name{s}));
            
            % Get the raw data files
            V=SPM.xY.VY;
            VresMS = spm_vol(fullfile(glmDirSubj,'ResMS.nii'));
            
            % Get time series data
            for r = 1:length(R)
                Y = region_getdata(V,R{r});  % Data is N x P
                resMS = region_getdata(VresMS,R{r});
                filename=(fullfile(baseDir,regDir,'data',subj_name{s},sprintf('rawts_%s.mat',R{r}.name)));
                save(filename,'Y','resMS','-v7.3');
                fprintf('Raw ts saved for %s for %s \n',subj_name{s},R{r}.name);
            end
        end
    case 'test_GLM' 
        % Get crossval R2 and R from GLM for different designs
        % As an evaluation, it uses only the differences between the
        % task-related regressors across all voxel
        % Input arguments :
        %    sn: Subject number
        %    roi: Cell array of ROI names 
        %    inK: List of terms in the pre-filtering matrix (cell array)
        %    inX: List of terms in the design matrix
        %    reg: Regression method: OLS, GLS, Ridge_pcm, Tikhonov_pcm..
        % Output arguments: 
        %   D: Data Frame holding the evaluation results with fields
        %   D.sn: Subject number
        %   D.roi: Cell array of ROI names 
        %   D.run: Run that served as test set
        %   D.R: Predictive correlation (across voxels and tasks)
        %   D.R2: Predictive correlation (across voxel and tasks) 
        sn = 2;
        glm = 1;
        roi = {'dentate'};%{'cerebellum_gray','dentate','brainstem','pontine'};
        inK = {};                   % Terms in filtering matrix - except for intercept
        inX = {{'Tasks','Instruct'}}; % Terms in the design matrix
        reg = {'OLS'};  % Regression methods
        runs = [1:16];
        D = []; % Result structure
        % Get the posible options to test
        vararginoptions(varargin,{'sn','roi','inK','inX','reg','evalX','runs'});
        
        % Load SPM file
        for s = sn
            glmDir = fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm));
            glmDirSubj = fullfile(glmDir, subj_name{s});
            load(fullfile(glmDirSubj,'SPM.mat'));
            INFO = load(fullfile(glmDirSubj,'SPM_info.mat'));
            nRuns = length(SPM.nscan);
            
            % Get the Data (Y)
            for r=1:length(roi) % Loop over ROIs
                fprintf('SN: %d, ROI: %s\n',s,roi{r});
                % Load raw time series
                load(fullfile(baseDir,regDir,'data',subj_name{s},sprintf('rawts_%s.mat',roi{r})));
                % Voxel-wise prewhiten the data
                Y = bsxfun(@rdivide,Y,sqrt(resMS));
                checksum = sum(abs(Y),1);
                
                badindx = isnan(checksum) | checksum==0;
                if sum(badindx)>0
                    warning('%d Nans or 0 in ts file',sum(badindx));
                    Y = Y(:,~badindx);
                end
                
                for model = 1:length(inX)
                    if (inX{model}{1}~='Tasks')
                        error('First group of regressors in design matrix needs to be Tasks'); 
                    end
                    
                    % Add regressors of no interest to X0
                    X0 = [];
                    if (~isempty(inK) && ~isempty(inK{model}))
                        for t=1:length(inK{model})
                            X0 = [X0 get_feature(inK{model}{t},s,SPM,INFO,1,1,1)];
                        end
                    end
                    
                    % Add intercept to X0
                    X0 = [X0 SPM.xX.X(:,SPM.xX.iB)];
                    
                    % Get the design matrix (only task related and intercepts)
                    X = [];
                    group = []; % regressor group
                    row = [];
                    i = 0 ;
                    
                    % Add regressor of no interest to X
                    for t=1:length(inX{model})
                        x = get_feature(inX{model}{t},s,SPM,INFO,0,1,1);
                        k = size(x,2);
                        indx = [i+1:i+k];
                        group(indx)=t;
                        X(:,indx) = x;
                        i = i + k;
                    end
                    N = size(X,1);
                    
                    % Get weight/whitening matrix
                    W = SPM.xX.W;
                    
                    % Make run indicator
                    row=zeros(N,1);
                    for rn=1:nRuns
                        row(SPM.Sess(rn).row,1)=rn;
                    end;
                    row(~ismember(row,runs))=0;
                    
                    for method=1:length(reg) % Loop over different methods
                                                
                        % Filtering design matrix and data
                        switch(reg{method})
                            case {'OLS'}
                                R         = eye(N)-X0*pinv(X0);
                                Xr        = R*X;    % KWX
                                Yr        = R*Y;
                            case {'GLS','ridge_fixed','tikhonov_pcm'}
                                R         = eye(N)-X0*pinv(X0);
                                Yr        = R*W*Y;
                                Xr        = R*W*X;
                        end
                        fprintf('%s:',reg{method});
                        % Crossvalidated approach
                        for rn=runs
                            
                            trainI = find(row~=rn & row>0);
                            testI  = find(row==rn);
                            
                            tic;
                            % Now estimate with the favorite regression approach
                            switch(reg{method})
                                case {'OLS','GLS','WLS','GLS_WLS'}
                                    Btrain = pinv(Xr(trainI,:))*Yr(trainI,:);
                                    Btest  = pinv(Xr(testI, :))*Yr(testI,:);
                                    theta = [];
                                case 'ridge_fixed'
                                    alpha = 1;
                                    Xtrain = Xr(trainI,:);
                                    Xtest  = Xr(testI, :);
                                    Btrain = (Xtrain' * Xtrain + eye(size(Xtrain,2))* alpha) \ (Xtrain' * Yr(trainI,:));
                                    Btest  = (Xtest'  * Xtest  + eye(size(Xtest,2)) * alpha) \ (Xtest'  * Yr(testI,:));
                                    theta = [];
                                case 'ridge_pcm'
                                    group0 = ones(1,size(Xr,2));
                                    if (rn==1)
                                        [theta,fitINFO]=pcm_fitModelRegression(Xr,Yr,group0,X0);
                                    end
                                    Btrain = pcm_estimateRegression(Xr(trainI,:),Yr(trainI,:),group0,X0(trainI,:),theta);
                                    Btest  = pcm_estimateRegression(Xr(testI,:), Yr(testI,:), group0,X0(testI,:), theta);
                                case 'tikhonov_pcm'
                                    if (rn==1)
                                        [theta,fitINFO]=pcm_fitModelRegression(Xr,Yr,group,X0);
                                    end
                                    Btrain = pcm_estimateRegression(Xr(trainI,:),Yr(trainI,:),group,X0(trainI,:),theta);
                                    Btest  = pcm_estimateRegression(Xr(testI,:), Yr(testI,:), group,X0(testI,:), theta);
                            end
                            fprintf('.');
                            % Performance valuation using only task-related regressors
                            time = toc;
                            evindx = find(group==1);
                            if ~isempty(evindx)
                                q = length(evindx); 
                                C = eye(q)-ones(q)/q; 
                                Ypred = C*Btrain(evindx,:);
                                Ytest = C*Btest(evindx,:);
                                % Record performance
                                T.roi = roi(r);
                                T.run  = rn;
                                T.sn   = s;
                                T.method  = reg(method);
                                T.methodN = method;
                                T.model = model;
                                T.theta = nan(1,5);
                                T.theta(1:length(theta)) = theta;
                                T.time = time;
                                T.R     = sum(sum(Ypred.*Ytest))./sqrt(sum(sum(Ypred.*Ypred)).*sum(sum(Ytest.*Ytest))); % R with differences between task-related regressors
                                T.R2    = 1-sum(sum((Ypred-Ytest).^2))./sum(sum(Ytest.^2));  % R2-on the differences of task related regressors
                                D = addstruct(D,T);
                            end % evalX
                        end % runs
                        fprintf('\n');
                    end % regression methods
                end % Model terms
            end % ROI
        end % Subjects
        varargout = {D};

    case 'test_GLM_timeseries' 
        % Get crossval R2 and R from GLM for different designs
        % example: bsp_imana('test_GLM','inK',{'Hpass','CSF'...},'inX',{'Mov',...},'ridge',0.1);
        % This version is more complex, as it uses the predicted time
        % series: Use only if you really understand what's going on. 
        % Input arguments :
        %    sn: Subject number
        %    roi: Cell array of ROI names 
        %    inK: List of terms in the pre-filtering matrix (cell array)
        %    inX: List of terms in the design matrix
        %    reg: Regression method: OLS, GLS, Ridge_pcm, Tikhonov_pcm..
        %    eval: Which regressors to use from the design matrix ot
        %    evaluate
        % Output arguments: 
        %   
        sn = 2;
        glm = 1;
        roi = {'dentate'};%{'cerebellum_gray','dentate','brainstem','pontine'};
        inK = {};                   % Terms in filtering matrix - except for intercept
        inX = {{'Tasks','Instruct'}}; % Terms in the design matrix
        reg = {'OLS'};  % Regression methods
        evalX = {1}; % Evaluation on what term in inX
        runs = [1:16];
        D = []; % Result structure
        % Get the posible options to test
        vararginoptions(varargin,{'sn','roi','inK','inX','reg','evalX','runs'});
        
        % Load SPM file
        for s = sn
            glmDir = fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm));
            glmDirSubj = fullfile(glmDir, subj_name{s});
            load(fullfile(glmDirSubj,'SPM.mat'));
            INFO = load(fullfile(glmDirSubj,'SPM_info.mat'));
            nRuns = length(SPM.nscan);
            
            % Get the Data (Y)
            for r=1:length(roi) % Loop over ROIs
                fprintf('SN: %d, ROI: %s\n',s,roi{r});
                % Load raw time series
                load(fullfile(baseDir,regDir,'data',subj_name{s},sprintf('rawts_%s.mat',roi{r})));
                % Voxel-wise prewhiten the data
                Y = bsxfun(@rdivide,Y,sqrt(resMS));
                checksum = sum(abs(Y),1);
                
                badindx = isnan(checksum) | checksum==0;
                if sum(badindx)>0
                    warning('%d Nans or 0 in ts file',sum(badindx));
                    Y = Y(:,~badindx);
                end
                
                for model = 1:length(inX)
                    
                    % Add regressors of no interest to X0
                    X0 = [];
                    if (~isempty(inK) && ~isempty(inK{model}))
                        for t=1:length(inK{model})
                            X0 = [X0 get_feature(inK{model}{t},s,SPM,INFO,1,1,1)];
                        end
                    end
                    
                    % Add intercept to X0
                    X0 = [X0 SPM.xX.X(:,SPM.xX.iB)];
                    
                    % Get the design matrix (only task related and intercepts)
                    X = [];
                    group = []; % regressor group
                    row = [];
                    i = 0 ;
                    
                    % Add regressor of no interest to X
                    for t=1:length(inX{model})
                        x = get_feature(inX{model}{t},s,SPM,INFO,0,1,1);
                        k = size(x,2);
                        indx = [i+1:i+k];
                        group(indx)=t;
                        X(:,indx) = x;
                        i = i + k;
                    end
                    N = size(X,1);
                    
                    % Get weight/whitening matrix
                    W = SPM.xX.W;
                    
                    % Make run indicator
                    row=zeros(N,1);
                    for rn=1:nRuns
                        row(SPM.Sess(rn).row,1)=rn;
                    end;
                    row(~ismember(row,runs))=0;
                    
                    for method=1:length(reg) % Loop over different methods
                                                
                        % Filtering design matrix and data
                        switch(reg{method})
                            case {'OLS'}
                                R         = eye(N)-X0*pinv(X0);
                                Xr        = R*X;    % KWX
                                Yr        = R*Y;
                            case {'GLS','ridge_fixed','tikhonov_pcm'}
                                R         = eye(N)-X0*pinv(X0);
                                Yr        = R*W*Y;
                                Xr        = R*W*X;
                        end
                        fprintf('%s:',reg{method});
                        % Crossvalidated approach
                        for rn=runs
                            
                            trainI = find(row~=rn & row>0);
                            testI  = find(row==rn);
                            
                            tic;
                            % Now estimate with the favorite regression approach
                            switch(reg{method})
                                case {'OLS','GLS','WLS','GLS_WLS'}
                                    Btrain = pinv(Xr(trainI,:))*Yr(trainI,:);
                                    Btest  = pinv(Xr(testI, :))*Yr(testI,:);
                                    theta = [];
                                case 'ridge_fixed'
                                    alpha = 1;
                                    Xtrain = Xr(trainI,:);
                                    Xtest  = Xr(testI, :);
                                    Btrain = (Xtrain' * Xtrain + eye(size(Xtrain,2))* alpha) \ (Xtrain' * Yr(trainI,:));
                                    Btest  = (Xtest'  * Xtest  + eye(size(Xtest,2)) * alpha) \ (Xtest'  * Yr(testI,:));
                                    theta = [];
                                case 'ridge_pcm'
                                    group0 = ones(1,size(Xr,2));
                                    if (rn==1)
                                        [theta,fitINFO]=pcm_fitModelRegression(Xr,Yr,group0,X0);
                                    end
                                    Btrain = pcm_estimateRegression(Xr(trainI,:),Yr(trainI,:),group0,X0(trainI,:),theta);
                                    Btest  = pcm_estimateRegression(Xr(testI,:), Yr(testI,:), group0,X0(testI,:), theta);
                                case 'tikhonov_pcm'
                                    if (rn==1)
                                        [theta,fitINFO]=pcm_fitModelRegression(Xr,Yr,group,X0);
                                    end
                                    Btrain = pcm_estimateRegression(Xr(trainI,:),Yr(trainI,:),group,X0(trainI,:),theta);
                                    Btest  = pcm_estimateRegression(Xr(testI,:), Yr(testI,:), group,X0(testI,:), theta);
                            end
                            fprintf('.');
                            % Performance valuation using only task-related regressors
                            time = toc;
                            for ev=1:length(evalX)
                                evindx = find(ismember(group,evalX{ev}));
                                if ~isempty(evindx)
                                    q = length(evindx); 
                                    C = eye(q)-ones(q)/q; 
                                    Ypred = Xr(testI,evindx)*Btrain(evindx,:);
                                    Ypredc = Xr(testI,evindx)*C*Btrain(evindx,:);
                                    Ytestp = Xr(testI,evindx)*Btest(evindx,:);
                                    Ytest  = Yr(testI,:);
                                    % Record performance
                                    T.roi = roi(r);
                                    T.run  = rn;
                                    T.sn   = s;
                                    T.method  = reg(method);
                                    T.methodN = method;
                                    T.evalX = ev;
                                    T.model = model;
                                    T.theta = nan(1,5);
                                    T.theta(1:length(theta)) = theta;
                                    T.time = time;
                                    T.R     = sum(sum(Ypred.*Ytest))./sqrt(sum(sum(Ypred.*Ypred)).*sum(sum(Ytest.*Ytest))); % R with timeseries 
                                    T.Rc     = sum(sum(Ypredc.*Ytest))./sqrt(sum(sum(Ypredc.*Ypredc)).*sum(sum(Ytest.*Ytest))); % R with timeseries - only contrast
                                    T.Rp     = sum(sum(Ypred.*Ytestp))./sqrt(sum(sum(Ypred.*Ypred)).*sum(sum(Ytestp.*Ytestp))); % R of predicted time series 
                                    T.R2    = 1-sum(sum((Ypred-Ytest).^2))./sum(sum(Ytest.^2));
                                    D = addstruct(D,T);
                                end
                            end % evalX
                        end % runs
                        fprintf('\n');
                    end % regression methods
                end % Model terms
            end % ROI
        end % Subjects
        varargout = {D};
    
    case 'F-test'                     % F-test to compare between p01 and p02
        % example: bsp_imana('F-test',1);
        sn = varargin{1};
        glm = 1;
        glmDir = fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm));
        
        % Load SPM file
        glmDirSubj = fullfile(glmDir, subj_name{sn});
        load(fullfile(glmDirSubj,'SPM.mat'));
        
        % F contrast for all regressor
        con = eye(length(SPM.Sess)*length(SPM.Sess(1).U));  % include each regressors
        con = [con;zeros(length(SPM.Sess),length(con))];    % add zeros for the intercepts
        SPM.xCon=spm_FcUtil('Set','SignalEffects', 'F', 'c',con,SPM.xX.xKXs);
        spm_contrasts(SPM);        
    
    case 'test_GLM_lowfreq'
        % Compare different methods to deal with auto-correlated noise
        % And sporardic artifacts
        sn = [2];
        model = {{'Tasks','Instruct'},...
            {'Tasks','Instruct'}};
        inK   = {{'Hpass'},...
            {}};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'OLS','GLS'};
        
        D=bsp_glm('test_GLM','roi',roi,'reg',method,'inX',model,'inK',inK,...
            'runs',[1:16],'sn',sn);
        save(fullfile(baseDir,'results','test_GLM_lowfreq.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_lowfreq'
        D=load('test_GLM_lowfreq.mat');
        sn = [2 3]; 
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1]}; 
        style={':',':','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.R,'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'HpassOLS','OLS','HpassGLS','GLS'},'facecolor',color); 
            title('Different ROIs');
            ylabel(sprintf('SN %d',sn(s)));
            subplot(num_subj,2,(s-1)*2+2); 
            lineplot(D.run,D.R,'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'HpassOLS','OLS','HpassGLS','GLS'},'linecolor',color,...
                    'linestyle',style,'linewidth',2); % {'HpassOLS','HpassGLS','OLS','GLS'} 
            title('Different Runs across ROIs');
        end; 
    
    case 'test_GLM_Physio'
        sn = [2 3];
        model = {{'Tasks','InstructC','Retro_HR'},...
            {'Tasks','InstructC','Retro_RESP'},...
            {'Tasks','InstructC','HR'},...
            {'Tasks','InstructC','RV'},...
            {'Tasks','InstructC'}};
        inK   = {{},{},{},{},{}};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','roi',roi,'reg',method,'inX',model,'inK',inK,...
            'runs',[1:16],'sn',sn);
        save(fullfile(baseDir,'results','test_GLM_physio.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_Physio'
        D=load('test_GLM_physio.mat');
        
        sn = [2 3]; 
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        % style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.Rc,'split',[D.model],'subset',D.sn==sn(s) & D.evalX==1,'leg',{'Retro HR','Retro Resp','HR','RV','none'},'facecolor',color); 
            title('performance of task differences');
            ylabel(sprintf('SN %d',sn(s)));
            subplot(num_subj,2,(s-1)*2+2); 
            barplot(D.roi,D.R,'split',[D.methodN D.model],'subset',D.sn==sn(s) & D.evalX==2,'leg',{'Retro HR','Retro Resp','HR','RV'},'facecolor',color); 
                    % {'HpassOLS','HpassGLS','OLS','GLS'} 
            title('R of Physio-regressor');
        end; 
    
    case 'test_GLM_script'
        model = {{'Tasks','Instruct'},...
            {'Tasks','Instruct','Retro_HR'},...
            {'Tasks','Instruct','Retro_RESP'},...
            {'Tasks','Instruct','HR'},...
            {'Tasks','Instruct','RV'}};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'OLS','GLS','ridge_pcm','tikhonov_pcm'};
        
        D=bsp_glm('test_GLM','roi',roi,'reg',method,'inX',model,'evalX',{[1 2]},'runs',[1:10]);
        save(fullfile(baseDir,'results','test_GLM_5.mat'),'-struct','D');
        varargout={D};
        
    case 'test_GLM_Physio_Filter'
        sn = [2 3];
        model = {{'Tasks','InstructC'},...
            {'Tasks','InstructC'},...
            {'Tasks','InstructC'},...
            {'Tasks','InstructC'},...
            {'Tasks','InstructC'}};
        inK   = {{'Retro_HR'},{'Retro_RESP'},{'HR'},{'RV'},{}};
        evalX = {[1]};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','roi',roi,'reg',method,'inX',model,'inK',inK,'evalX',evalX,...
            'runs',[1:16],'sn',sn);
        save(fullfile(baseDir,'results','test_GLM_physio_filter.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_Physio_Filter'
        D=load('test_GLM_physio_filter.mat');
        
        sn = [2 3]; 
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        % style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.Rc,'split',[D.model],'subset',D.sn==sn(s),'leg',{'Retro HR','Retro Resp','HR','RV','none'},'facecolor',color); 
            title('performance of task differences');
            ylabel(sprintf('SN %d',sn(s)));
            subplot(num_subj,2,(s-1)*2+2); 
            barplot(D.roi,D.R,'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'Retro HR','Retro Resp','HR','RV'},'facecolor',color); 
                    % {'HpassOLS','HpassGLS','OLS','GLS'} 
            title('R of Physio-filter');
        end; 
        
    case 'test_GLM_CSF'
        sn = [2 3];
        model = {{'Tasks','InstructC','CSF'},...
            {'Tasks','InstructC','CSFPCAindiv'},...
            {'Tasks','InstructC','CSFPCAall'},...
            {'Tasks','InstructC'}};
        inK   = {{},{},{},{}};
        evalX = {[1]};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','roi',roi,'reg',method,'inX',model,'inK',inK,'evalX',evalX,...
            'runs',[1:16],'sn',sn);
        save(fullfile(baseDir,'results','test_GLM_csf.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_CSF'
        D=load('test_GLM_csf.mat');
        
        sn = [2 3]; 
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        % style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.Rc,'split',[D.model],'subset',D.sn==sn(s),'leg',{'CSF','CSF PCAindiv','CSF PCAall','none'},'facecolor',color); 
            title('performance of task differences');
            ylabel(sprintf('SN %d',sn(s)));
            subplot(num_subj,2,(s-1)*2+2); 
            barplot(D.roi,D.R,'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'CSF','CSF PCAindiv','CSF PCAall','none'},'facecolor',color); 
                    % {'HpassOLS','HpassGLS','OLS','GLS'} 
            title('R of CSF');
        end; 
    
    case 'test_GLM_CSF_Filter'
        sn = [2 3];
        model = {{'Tasks','InstructC'},...
            {'Tasks','InstructC'},...
            {'Tasks','InstructC'},...
            {'Tasks','InstructC'}};
        inK   = {{'CSF'},{'CSFPCAindiv'},{'CSFPCAall'},{}};
        evalX = {[1]};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','roi',roi,'reg',method,'inX',model,'inK',inK,'evalX',evalX,...
            'runs',[1:16],'sn',sn);
        save(fullfile(baseDir,'results','test_GLM_csf_filter.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_CSF_Filter'
        D=load('test_GLM_csf_filter.mat');
        
        sn = [2 3]; 
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        % style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.Rc,'split',[D.model],'subset',D.sn==sn(s),'leg',{'CSF','CSF PCAindiv','CSF PCAall','none'},'facecolor',color); 
            title('performance of task differences');
            ylabel(sprintf('SN %d',sn(s)));
            subplot(num_subj,2,(s-1)*2+2); 
            barplot(D.roi,D.R,'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'CSF','CSF PCAindiv','CSF PCAall','none'},'facecolor',color); 
                    % {'HpassOLS','HpassGLS','OLS','GLS'} 
            title('R of CSF-filter');
        end; 
        
    case 'physio_reg' % Examines the covariance of physiological regressors with main regressors 
        sn = 2; 
        glm = 1; 
        vararginoptions(varargin,{'sn'})
        reg = {'Tasks','InstructC','Retro_HR','Retro_RESP','HR','RV'};

        for s = sn
            glmDir = fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm));
            glmDirSubj = fullfile(glmDir, subj_name{s});
            load(fullfile(glmDirSubj,'SPM.mat'));
            INFO = load(fullfile(glmDirSubj,'SPM_info.mat'));
            nRuns = length(SPM.nscan);
            
            for i=1:nRuns 
                tr(SPM.Sess(i).row,1)=[1:SPM.nscan(i)]; 
                rn(SPM.Sess(i).row,1)=i; 
            end; 
    
            for r = 1:length(reg)
                X{r}=get_feature(reg{r},s,SPM,INFO,0,1,1); 
            end 
            
            onsets = vertcat(SPM.Sess(1).U([1:2:18]).ons); 
            % Average run-averaged 
            to_plot=[1,2,5,6]; 
            for i = 1:4 
                subplot(4,2,(i-1)*2+1);
                Y=reshape(sum(X{to_plot(i)},2),SPM.nscan(1),16); 
                traceplot([1:SPM.nscan(1)],Y','errorfcn','stderr'); 
                drawlines(onsets,'k'); 
                set(gca,'XLim',[-3 SPM.nscan(1)+3]);
            end; 
            XI = [X{1} X{2}]; 
            to_evaluate=[3,4,5,6]; 
            for i = 1:4
                for r=1:16 
                    indx = find(rn==r); 
                    b(:,r)=pinv(XI(indx,:))*X{to_evaluate(i)}(indx,1);
                end;
                subplot(4,2,(i-1)*2+2);
                myboxplot([],b');
                drawline(0,'dir','horz'); 
                set(gca,'XTickLabel',{'VS','AO','FE','FS','TOM','NB','SP','R','RM','In'});
                ylabel(reg{to_evaluate(i)});
            end; 
            
            keyboard; 
        end
end
end

function XX=get_feature(what,sn,SPM,INFO,separate,sscale,zscale)
% Function that gets the regressors
global baseDir
global subj_name
glm = 1;
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
    case 'HR'    % Retroicor of cardio and resp 18 comp
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
        XX(SPM.Sess(rn).row,:) = X{rn};
    end
end
end
