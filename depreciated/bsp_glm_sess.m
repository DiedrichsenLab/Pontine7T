function varargout=bsp_glm(what,varargin)
% Second level analysis for the pontine7T project, starting from the level
% of the GLM.

%

numDummys = 3;                                                              % per run
numTRs    = 328;                                                            % per run (includes dummies)
%========================================================================================================================
% PATH DEFINITIONS
global baseDir; 
global subj_name; 

baseDir         ='/srv/diedrichsen/data/Pontine7T/op-coreg-1';
% baseDir         ='/Volumes/diedrichsen_data$/data/Pontine7T';
imagingDir      ='/imaging_data';
imagingDirRaw   ='/imaging_data_raw';
anatomicalDir   ='/anatomicals';
suitDir         ='/suit';
regDir          ='/RegionOfInterest';
%========================================================================================================================
% PRE-PROCESSING
subj_name = {'S99'};
%========================================================================================================================
% GLM INFO
sess        = {'sess1','sess2'};
funcRunNum  = [1,8];  % first and last behavioural run numbers
%funcRunNum  = [9,16];
run         = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
runB        = [1,2,3,4,5,6,7,8];  % Behavioural labelling of runs
%runB        = [9,10,11,12,13,14,15,16];

%========================================================================================================================
switch(what)
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
        % EXAMPLE: bsp_imana('GLM:glm1',[1:XX],1,[1:XX])
        sn=varargin{1};
        sessNum=varargin{2};
        runs=varargin{3}; % [1:8]
        
        announceTime=5;  % Instruction duration
        glm=1;
        subjs=length(sn);
        
        % load in task information
        C=dload(fullfile(baseDir,'pontine_taskConds_GLM.tsv'));
        Cc=getrow(C,C.StudyNum==1 & C.condNum==1); % Only get the first condition for each tasks
        nTask      = max(Cc.taskNum);                                 % how many tasks there are?
        
        for s=1:subjs,
            T=[];
            A = dload(fullfile(baseDir,'data','fmri',subj_name{sn(s)},['fmri_',subj_name{sn(s)},'_op','.tsv'])); % get scanning timing and order
            A = getrow(A,A.run_num>=funcRunNum(1) & A.run_num<=funcRunNum(2)); % get only the runs we need (remove any test or Behav training)
            
            glmSubjDir =[baseDir, sprintf('/GLM_firstlevel_%d/',glm), sprintf('%s/',sess{sessNum}), subj_name{sn(s)}];dircheck(glmSubjDir); % Create GLM folder
            
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
        C=dload(fullfile(baseDir,'bsp_taskConds_GLM.txt'));
        Cc=getrow(C,C.StudyNum==1);
        Tasks = unique(Cc.taskNames,'rows','stable');                       % get the task names
        nTask      = unique(length(Tasks));                                 % how many tasks there are?
        
        for s=1:subjs,
            T=[];
            A = dload(fullfile(baseDir,'data', subj_name{sn(s)},['sc1_',subj_name{sn(s)},'.dat'])); % get scanning timing and order
            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2)); % get only the runs we need (remove any test or Behav training)
            
            glmSubjDir =[baseDir, sprintf('/GLM_firstlevel_%d/',glm),subj_name{sn(s)}];dircheck(glmSubjDir); % Create GLM folder
            
            % Fill up struct for glm
            J.dir = {glmSubjDir};
            J.timing.units = 'secs';
            J.timing.RT = 1.0;
            J.timing.fmri_t = 16;
            J.timing.fmri_t0 = 1;
            
            for r=1:numel(runs) % loop through runs
                P=getrow(A,A.runNum==runB(r));
                for i=1:(numTRs-numDummys) % get the filenames of the nifti volumes for each run
                    N{i} = [fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('run_%2.2d.nii,%d',runs(r),i))];
                end;
                J.sess(r).scans= N; % number of scans in run
                
                % loop through tasks
                for it = 1:nTask
                    % The order of tasks are different for each run, to
                    % have a common order for the tasks, I will be reading
                    % from the Cc file for all the runs and subjects
                    ST = find(strcmp(P.taskName,Tasks{it}));
                    
                    % get the task onset (instruction onset + announceTime)
                    onset = P.realStartTime(ST) - J.timing.RT*numDummys +announceTime;
                    
                    % filling in the fields for SPM_info.mat
                    S.task      = it;
                    S.TN        = {Tasks{it}};
                    S.inst      = 0;
                    S.instOrder = 0;
                    S.time      = onset;
                    T  = addstruct(T, S);
                    
                    % filling in the fields for SPM.mat
                    J.sess(r).cond(it).name     = Tasks{it};
                    J.sess(r).cond(it).onset    = onset;
                    J.sess(r).cond(it).duration = 30;             % each task lasts for 30 sec
                    J.sess(r).cond(it).tmod     = 0;
                    J.sess(r).cond(it).orth     = 0;
                    J.sess(r).cond(it).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                    
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
            J.mask = {fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'c_anatomical_pcereb.nii')};
            J.mthresh = 0.05;
            J.cvi_mask = {fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'c_anatomical_pcereb.nii')};
            J.cvi =  'fast';
            
            spm_rwls_run_fmri_spec(J);
            
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            fprintf('glm_%d has been saved for %s \n',glm, subj_name{sn(s)});
        end
    case 'GLM:estimate'               % Estimate GLM depending on subjNum & glmNum
        % example: bsp_imana('GLM:estimate_glm',1,1,1)
        sn=varargin{1};
        glm=varargin{2};
        sessNum=varargin{3};
        
        subjs=length(sn);
        
        for s=1:subjs,
            glmDir =[baseDir,sprintf('/GLM_firstlevel_%d/',glm),sprintf('%s/',sess{sessNum}),subj_name{sn(s)}];
            load(fullfile(glmDir,'SPM.mat'));
            SPM.swd=glmDir;
            spm_rwls_spm(SPM);
        end
    case 'GLM:contrast'               % Create Contrast images
        %%% Calculating contrast images.
        % 'SPM_light' is created in this step (xVi is removed as it slows
        % down code for FAST GLM).
        % Example1: bsp_imana('GLM:contrast', 'sn', 3, 'glm', 1, 'type', 'task')
        % Example2: bsp_imana('GLM:contrast', 'sn', 3, 'glm', 1, 'type', 'cond')
        
        sn             = 3;             %% list of subjects
        glm            = 1;             %% The glm number
        con_vs         = 'rest_task';   %'average_4';   %% set it to 'rest' or 'average' (depending on the contrast you want)
        type           = 'task';        %% it can be set to either cond or task.
        
        vararginoptions(varargin, {'sn', 'glm', 'con_vs', 'type'})
        
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, sprintf('GLM_firstlevel_%d', glm));
        
        for s = sn
            fprintf('******************** calculating contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'))
            
            SPM  = rmfield(SPM,'xCon');
            cd(fullfile(glmDir, subj_name{s}))
            T = load('SPM_info.mat');
            
            % t contrast for tasks
            ucondition = unique(T.(type));
            idx = 1;
            for tt = 1:length(ucondition) % 0 is "instruct" regressor
                switch con_vs
                    case 'rest_task' % contrast against rest
                        con = zeros(1,size(SPM.xX.X,2));
                        con(:,logical(T.(type) == ucondition(tt))) = 1;
                        tmp = zeros(1, size(SPM.xX.X, 2));
                        tmp(:, strcmp(T.taskName, 'rest')) = 1;
                        sum_rest = sum(tmp);
                        tmp      = tmp./sum_rest;
                        con      = con/abs(sum(con));
                        con      = con - tmp;
                    case 'average_1' % contrast against the average of the other tasks including the instructions
                        % New: Eva, Oct 2nd
                        con                     = zeros(1,size(SPM.xX.X, 2));
                        con(1,logical(T.(type) == ucondition(tt))) = 1./sum(logical(T.(type) == ucondition(tt)));
                        con(1,logical(T.(type) ~= ucondition(tt))) = -1./sum(logical(T.(type) ~= ucondition(tt)));
                    case 'average_2' % contrast against the average of the other tasks not including the instructions
                        con        = zeros(1,size(SPM.xX.X, 2));
                        % TO TRY below - no instructions as contrasts
                        con(1,logical(T.(type) == ucondition(tt))) = 1./sum(logical(T.(type) == ucondition(tt)));
                        con(1,logical((T.(type) ~= ucondition(tt)) .* (T.inst == 0))) = -1./sum(logical((T.(type) ~= ucondition(tt)) .* (T.inst == 0)));
                    case 'average_4' % contrast against the average of all the tasks (not including the instructions)
                        con        = zeros(1,size(SPM.xX.X, 2));
                        % TO TRY below - no instructions as contrasts
                        con(1,logical(T.(type)      == ucondition(tt))) = 1./sum(logical(T.(type) == ucondition(tt)));
                        con(1, logical(T.inst == 0)) = con(1, logical(T.inst == 0)) - 1./sum(T.inst == 0);
                    case 'rest'
                        con                     = zeros(1,size(SPM.xX.X,2));
                        con(:,logical(T.(type) == ucondition(tt))) = 1;
                        con                     = con/abs(sum(con));
                end
                
                name = sprintf('%s-%s',char(unique(T.taskName(T.(type) == ucondition(tt)))), con_vs);
                
                SPM.xCon(idx) = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);
                idx=idx+1;
            end % tt (conditions)
            SPM = spm_contrasts(SPM,1:length(SPM.xCon));
            save('SPM.mat', 'SPM','-v7.3');
            SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
            save(fullfile(glmDir, subj_name{s}, 'SPM_light.mat'), 'SPM');
            
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
        % Example1: bsp_glm('GLM:Fcontrast', 'sn', 3, 'glm', 1, 'sesn', 1, 'type', 'task')
        % Example2: bsp_glm('GLM:Fcontrast', 'sn', 3, 'glm', 1, 'sesn', 1, 'type', 'cond')
        
        sn             = 3;             %% list of subjects
        glm            = 1;             %% The glm number
        sesn           = 1;             %% session label
        type           = 'task';        %% it can be set to task, ....
        
        vararginoptions(varargin, {'sn', 'glm', 'sesn', 'type'})
        
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, sprintf('GLM_firstlevel_%d', glm),sprintf('%s/',sess{sesn}));
        
        for s = sn
            fprintf('******************** calculating contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'))
            
            SPM  = rmfield(SPM,'xCon');
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
        % EXAMPLE: bsp_glm('GLM:contrast_F_summary',[1:XX],1,1)
        sn=varargin{1};
        sessNum=varargin{2};           
        glm = varargin{3};           %% The glm number
        
        D=[]; 
        subjs=length(sn);
        
        for s = 1:subjs
            glmSubjDir =[baseDir, sprintf('/GLM_firstlevel_%d/',glm), sprintf('%s/',sess{sessNum}), subj_name{sn(s)}];
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
    
    case 'PHYS:extract'               % Extract puls and resp files from dcm
        sn=varargin{1};
        
        cd (fullfile(baseDir,'physio',subj_name{sn}));
        logfiles = dir('*.dcm');
        
        for lFile = 1:length(logfiles)
            fname = strsplit(logfiles(lFile).name,'.');
            mkdir(fname{1});
            movefile(logfiles(lFile).name,fname{1})
            cd (fname{1})
            extractCMRRPhysio(logfiles(lFile).name);
            cd ..
        end
    case 'PHYS:createRegressor'        % Create Retroicor regressors using TAPAS (18 components)
        sn=1; 
        run = [1:16]; 
        stop = true; 
        vararginoptions(varargin,{'sn','run','stop'}); 
        
        for nrun= run
            
            logDir = fullfile(baseDir,'physio',subj_name{sn},sprintf('run%2.2d',nrun));
            cd (logDir);
            PULS = dir('*PULS.log');
            RSP  = dir('*RESP.log');
            log  = dir('*Info.log');
            
            % 2. Define physio model structure
            physio = tapas_physio_new();
            
            % 3. Define input files
            physio.save_dir = {logDir};                         % enter directory to save physIO output
            physio.log_files.cardiac = {PULS.name};             % .puls file
            physio.log_files.respiration = {RSP.name};          % .resp file
            physio.log_files.scan_timing = {log.name};          % Log file
            physio.log_files.vendor = 'Siemens_Tics';           % Vendor
            physio.log_files.relative_start_acquisition = 0;    % only used for Philips
            physio.log_files.align_scan = 'last';               % sync log and EPI files from the last scan
            
            % 4. Define scan timing parameters
            physio.scan_timing.sqpar.Nslices = 60;              % number of slices in EPI volume
            physio.scan_timing.sqpar.TR = 1;                    % TR in secs
            physio.scan_timing.sqpar.Ndummies = numDummys;      % Real dummy are not recorded
            physio.scan_timing.sqpar.Nscans = numTRs-numDummys; % number of time points/scans/volumes
            physio.scan_timing.sqpar.onset_slice = 30;          % slice of interest (choose middle slice if unsure?)
            physio.scan_timing.sync.method = 'scan_timing_log'; % using info file
            
            % 5. Define cardiac data parameters
            physio.preproc.cardiac.modality = 'PPU';        % cardiac data acquired with PPU
            physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
            physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
            physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
            
            % 6. Define generic physio model parameters
            physio.model.orthogonalise = 'none';
            physio.model.output_multiple_regressors = 'multiple_regressors.txt';  %% regressor output filename
            physio.model.output_physio = 'physio.mat';  %% model output filename
            physio.model.censor_unreliable_recording_intervals = false;  %% set to true to automatically censor missing recording intervals
            
            % 7. RETROICOR phase expansion (Glover et al., 2000)
            physio.model.retroicor.include = true;  %% use RETROICOR to calculate regressors
            physio.model.retroicor.order.c = 3;
            physio.model.retroicor.order.r = 3;
            physio.model.retroicor.order.cr = 0;
            
            % 8. RVT respiratory volume per time (Birn et al., 2008)
            physio.model.rvt.include = true;
            physio.model.rvt.delays = 0;
            
            % 9. HRV heart rate response (Chang et al., 2009)
            physio.model.hrv.include = true;
            physio.model.hrv.delays = 0;
            
            % 10. ROI anatomical component-based noise extraction (aCompCor; Behzadi et al., 2007)
            physio.model.noise_rois.include = false;
            physio.model.noise_rois.thresholds = 0.9;
            physio.model.noise_rois.n_voxel_crop = 0;
            physio.model.noise_rois.n_components = 1;
            
            % 11. Motion modeling (Friston et al., 1996)
            physio.model.movement.include = false;
            physio.model.movement.file_realignment_parameters = {''};  %% No input required unless physio.model.movement.include is set to 'true'
            physio.model.movement.order = 6;
            physio.model.movement.outlier_translation_mm = 3;
            physio.model.movement.outlier_rotation_deg = Inf;
            
            % 12. Include other regressors
            physio.model.other.include = false;  %% include realignment parameters
            physio.model.other.input_multiple_regressors = {''};  %% rp*.txt file if needed
            
            % 13. Output/model options
            physio.verbose.level = 1;  %% 0 = suppress graphical output; 1 = main plots; 2 = debugging plots; 3 = all plots
            physio.verbose.process_log = cell(0, 1);
            physio.verbose.fig_handles = zeros(0, 1);
            physio.verbose.fig_output_file = [];  %% output figure filename (output format can be any matlab supported graphics format)
            physio.verbose.use_tabs = false;
            
            physio.ons_secs.c_scaling = 1;
            physio.ons_secs.r_scaling = 1;
            
            % 14. Run main script
            [physio_out,R,ons_sec]=tapas_physio_main_create_regressors(physio);
            % Save as different text files.... 
            dlmwrite('reg_retro_hr.txt',R(:,1:6));
            dlmwrite('reg_retro_resp.txt',R(:,7:12));
            dlmwrite('reg_hr.txt',[R(:,13) ons_sec.hr]); % Add the un-convolved version as well 
            dlmwrite('reg_rvt.txt',[R(:,14) ons_sec.rvt]); % add the un-convolved version as well 
            if (stop) 
                keyboard; 
                close all;
            end;
        end
    case 'ROI:csf_image'                 % Defines ROIs for brain structures
        sn=varargin{1}; 
        P = {fullfile(baseDir,'GLM_firstlevel_1','sess1',subj_name{sn},'mask.nii'),...
            fullfile(baseDir,'anatomicals',subj_name{sn},'c3anatomical.nii')}; 
        outname = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn},'csf_mask.nii'); 
        %spm_imcalc(char(P),outname,'i1 & (i2>0.3)',{0,0,2,1});
        spm_imcalc(char(P),outname,'i1 & (i2>0.3)',{0,0,1,4}); 
        
        % Before runing t
    case 'ROI:define'                 % Defines ROIs for brain structures
        % Before runing this, create masks for different structures
        sn=varargin{1}; % subjNum
        regions=varargin{2}; % Name of ROI, see cases below
        
        subjs=length(sn);
        for s=1:subjs,
            maskimg = fullfile(baseDir,'GLM_firstlevel_1','sess1',subj_name{sn(s)},'mask.nii');
            switch regions
                case 'cerebellum'
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'c_anatomical_pcereb.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['cerebellum'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                case 'cerebellum_grey'
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'c_anatomical_seg1.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['cerebellum_grey'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                case 'dentate'
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'dentate_mask.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['dentate'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                case 'pontine'
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'pontine_mask.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['pontine'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                case 'olive'
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'olive_mask.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['olive'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                case 'csf'
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'csf_mask_brainstem.nii');
                    R{1}.type = 'roi_image';
                    R{1}.mask = spm_vol(maskimg); 
                    R{1}.file= file;
                    R{1}.name = ['csf-bs'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                case 'csf_anterior'
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'rcsf_mask_anterior.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['csf-ant'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                case 'csf_posterior'
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'rcsf_mask_posterior.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['csf-post'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                
            end
            dircheck(fullfile(baseDir,regDir,'data',subj_name{sn(s)}));
            save(fullfile(baseDir,regDir,'data',subj_name{sn(s)},sprintf('regions_%s.mat',regions)),'R');
            region_saveasimg(R{1},maskimg,'name',fullfile(baseDir,regDir,'data',subj_name{sn(s)},sprintf('regions_%s.nii',regions))); 
            fprintf('ROIs have been defined for %s for %s \n',regions,subj_name{sn(s)})
        end
    case 'ROI:defall'                 % Run all types in ROI
        sn=varargin{1}; % subjNum
        bsp_imana('ROI:define',sn,'cerebellum');
        bsp_imana('ROI:define',sn,'cerebellum_grey');
        bsp_imana('ROI:define',sn,'dentate');
        bsp_imana('ROI:define',sn,'pontine');
        bsp_imana('ROI:define',sn,'olive');
        bsp_imana('ROI:define',sn,'csf');
        
    case 'TS:getRawTs'                % Get raw timeseries and save them 
        % bsp_imana('TS:getRawTs',1,1,1,'dentate');
        sn=varargin{1}; % subjNum
        glm=varargin{2}; % glmNum
        sessNum=varargin{3}; %sessNum
        regions=varargin{4}; % ROI

        glmDir =fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm),sprintf('%s',sess{sessNum}));
        subjs=length(sn);
        
        for s=1:subjs,
            glmDirSubj=fullfile(glmDir, subj_name{sn(s)});
            load(fullfile(glmDirSubj,'SPM.mat'));
            
            % load data
            load(fullfile(baseDir,regDir,'data',subj_name{sn(s)},sprintf('regions_%s.mat',regions)));
            SPM=spmj_move_rawdata(SPM,fullfile(baseDir,imagingDir,subj_name{sn(s)}));

            % Get the raw data files
            V=SPM.xY.VY;
            VresMS = spm_vol(fullfile(glmDirSubj,'ResMS.nii'));
            
            % Get time series data
            Y = region_getdata(V,R{1});  % Data is N x P
            resMS = region_getdata(VresMS,R{1});

            filename=(fullfile(baseDir,regDir,sprintf('glm%d',glm),sprintf('%s',sess{sessNum}),subj_name{sn(s)},sprintf('rawts_%s.mat',regions)));
            save(filename,'Y','resMS','-v7.3');
            fprintf('Raw ts saved for %s (%s) for %s \n',subj_name{sn(s)},sprintf('glm%d',glm),sprintf('%s/',sess{sessNum}),regions);
        end            
    case 'TS:getall'                  % Get all raw TS defined
        sn=varargin{1};
        %bsp_glm('TS:getRawTs',sn,1,'cerebellum');
        bsp_glm('TS:getRawTs',sn,1,'cerebellum_grey');
        bsp_glm('TS:getRawTs',sn,1,'dentate');
        bsp_glm('TS:getRawTs',sn,1,'csf');
        bsp_glm('TS:getRawTs',sn,1,'pontine');
        bsp_glm('TS:getRawTs',sn,1,'olive');      
        
    case 'test_GLM'                   % Get crossval R2 and R from GLM for different designs
        % example: bsp_imana('test_GLM','inK',{'Hpass','CSF'...},'inX',{'Mov',...},'ridge',0.1);
        % Input arguents :
        %    sn: Subject number
        %    inK: List of terms in the pre-filtering matrix (if any)
        %    inX: List of terms in the design matrix
        %
        sn = 5;
        glm = 1;
        roi = {'dentate'};%{'cerebellum_grey','dentate','brainstem','pontine'};
        inK = {};                   % Terms in filtering matrix - except for intercept
        inX = {{'Tasks','Instruct'}}; % Terms in the design matrix
        reg = {'OLS'};  % Regression methods
        evalX = {1,2,[1 2]}; % Evaluation on what term in inX
        runs = [1:20]; 
        D = []; % Result structure 
        % Get the posible options to test
        vararginoptions(varargin,{'roi','inX','inK','sn','reg','evalX','runs'});
        
        % Load SPM file
        glmDir = fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm));
        glmDirSubj = fullfile(glmDir, subj_name{sn});
        load(fullfile(glmDirSubj,'SPM.mat'));
        INFO = load(fullfile(glmDirSubj,'SPM_info.mat'));
        nRuns = length(SPM.nscan);
        
        % Get the Data (Y)
        for r=1:length(roi) % Loop over ROIs
            fprintf('SN: %d, ROI: %s\n',sn,roi{r});
            % Load raw time series
            load(fullfile(baseDir,regDir,sprintf('glm%d',glm),subj_name{sn},sprintf('rawts_%s.mat',roi{r})));
            % Voxel-wise prewhiten the data 
            Y = bsxfun(@rdivide,Y,sqrt(resMS));
            checksum = sum(abs(Y),1);
            badindx = isnan(checksum) | checksum==0; 
            if sum(badindx)>0
                warning('Nans or 0 in ts file'); 
                Y = Y(:,~badindx);
            end
            
            for model = 1:length(inX)
                % Add regressors of no interest to X0
                X0 = [];
                if (~isempty(inK))
                    for t=1:length(inK{model})
                        X0 = [X0 get_feature(inK{model}{t},sn,SPM,INFO,1,1,1)];
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
                    x = get_feature(inX{model}{t},sn,SPM,INFO,0,1,1);
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
                
                % Residual forming matrix
                R         = eye(N)-X0*pinv(X0);
                Xr        = R*X;    % KWX
                               
                % Whiten/Weight data and remove filter confounds
                Yr = R*Y;
                
                for method=1:length(reg) % Loop over different methods
                    fprintf('%s:',reg{method});
                    % Crossvalidated approach
                    for rn=runs
                        
                        trainI = find(row~=rn & row>0);
                        testI  = find(row==rn);
                        
                        tic; 
                        % Now estimate with the favorite regression approach
                        switch(reg{method})
                            case 'OLS'
                                Btrain = pinv(Xr(trainI,:))*Yr(trainI,:);
                                Btest  = pinv(Xr(testI, :))*Yr(testI,:);
                                theta = [];
                            case 'GLS'
                                Yw=W*Yr;
                                Xw=W*Xr;
                                Btrain = pinv(Xw(trainI,:))*Yw(trainI,:);
                                Btest  = pinv(Xw(testI, :))*Yw(testI,:);
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
                            indx = ismember(group,evalX{ev});
                            Ypred = Xr(testI,indx)*Btrain(indx,:);
                            Ytestp = Xr(testI,indx)*Btest(indx,:);
                            Ytest  = Yr(testI,:); 
                            % Record performance
                            T.roi = roi(r);
                            T.run  = rn;
                            T.sn   = sn;
                            T.method  = reg(method);
                            T.methodN = method;
                            T.evalX = ev;
                            T.model = model;
                            T.theta = nan(1,5); 
                            T.theta(1:length(theta)) = theta;
                            T.time = time; 
                            T.R     = mean(sum(Ypred.*Ytest)./sqrt(sum(Ypred.*Ypred).*sum(Ytest.*Ytest)));
                            T.Rp     = mean(sum(Ypred.*Ytestp)./sqrt(sum(Ypred.*Ypred).*sum(Ytestp.*Ytestp)));
                            T.R2    = 1-mean(sum((Ypred-Ytest).^2)./sum(Ytest.^2));
                            D = addstruct(D,T);
                        end % evalX
                    end % runs
                    fprintf('\n'); 
                end % regression methods
            end % Model terms 
            end % ROI
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
    case 'test_GLM_script'
        model = {{'Tasks','Instruct'},...
                 {'Tasks','Instruct','Retro_HR'},...
                 {'Tasks','Instruct','Retro_RESP'},...
                 {'Tasks','Instruct','HR'},...
                 {'Tasks','Instruct','RV'}}; 
        roi = {'pontine','dentate','olive','csf','cerebellum_grey','cortical_grey_left'}; 
        method = {'OLS','GLS','ridge_pcm','tikhonov_pcm'};
        
        D=bsp_glm('test_GLM','roi',roi,'reg',method,'inX',model,'evalX',{[1 2]},'runs',[1:10]);
        save(fullfile(baseDir,'results','test_GLM_5.mat'),'-struct','D'); 
        varargout={D}; 
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
            case 'Hpass'        % High pass filter
                k.HParam = 128;
                for rn = 1:nRuns
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
                    X{rn} = A(:,1:4);
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
                csf = load(fullfile(baseDir,'RegionOfInterest',sprintf('glm%d',glm),subj_name{sn},'rawts_csf.mat'));
                % mean CSF signal
                for rn = 1:nRuns
                    mcsf = mean(csf.Y(SPM.Sess(rn).row,:),2);
                    X{rn} = bsxfun(@rdivide,mcsf,sqrt(sum(mcsf.^2)));
                end
            case 'CSFPCAindiv'       % 2 Pcs of CSF
                % Get the CSF data
                csf = load(fullfile(baseDir,'RegionOfInterest',sprintf('glm%d',glm),subj_name{sn},'rawts_csf.mat'));
                % Compute the PCs
                for rn = 1:nRuns
                    runcsf = csf.Y(SPM.Sess(rn).row,:);
                    % Get the principal components
                    [~,score] = pca(runcsf);
                    X{rn} = score(:,1:2);
                end
            case 'CSFPCAall'   % 2 PCs of CSF computed over 4 runs
                % Get the CSF mask
                csf = load(fullfile(baseDir,'RegionOfInterest',sprintf('glm%d',glm),subj_name{sn},'rawts_csf.mat'));
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
                VAR(rn,:) = var(X{rn});
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
                
            end
        else
            XX = zeros(N,K);
            for rn = 1:nRuns
                XX(SPM.Sess(rn).row,:) = X{rn};
            end
        end
end
