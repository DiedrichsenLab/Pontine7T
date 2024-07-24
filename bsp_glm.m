function varargout=bsp_glm(what,varargin)
% GLM analysis for the pontine7T project, starting from preprocessed data
% (see bsp_imana).
% Define the data basedirectory 
if isdir('/Volumes/diedrichsen_data$/data')
    workdir='/Volumes/diedrichsen_data$/data';
elseif isdir('/srv/diedrichsen/data')
    workdir='/srv/diedrichsen/data';
else
    fprintf('Workdir not found. Mount or connect to server and try again.');
end
baseDir=(sprintf('%s/Cerebellum/Pontine7T',workdir));

imagingDir      ='imaging_data';
imagingDirRaw   ='imaging_data_raw';
anatomicalDir   ='anatomicals';
suitDir         ='suit';
regDir          ='RegionOfInterest';
fmapDir         ='fieldmaps';

% Load Participant information (make sure you have Dataframe/util in your
% path
pinfo = dload(fullfile(baseDir,'participants_new_format.tsv')); 
subj_name = pinfo.participant_id;
good_subj = find(pinfo.good)'; % Indices of all good subjects

%========================================================================================================================
% GLM INFO
numDummys = 3;                                                              % per run
numTRs    = pinfo.numTR;                                                             % per run (includes dummies)
runs         = {'01','02','03','04','05','06','07','08'};
runB        = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];  % Behavioural labelling of runs
sess        = [1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2];                  % session number
%========================================================================================================================
%========================================================================================================================
switch(what)
    case 'GLM:do'
        bsp_glm('GLM:glm1',[8],[1:16]);
        bsp_glm('GLM:glm2',[8],[1:16]);
        bsp_glm('GLM:estimate',[8],1);
        bsp_glm('GLM:estimate',[8],2);
        bsp_glm('GLM:Fcontrast','sn', [8], 'glm', 1, 'type', 'task');
        bsp_glm('GLM:Fcontrast','sn', [8], 'glm', 2, 'type', 'task');
    
    case 'GLM:makeMask' % Changed to include CSF
        sn=varargin{1}; % subjNum
        tissues = [1:3];
        
        P{1} = fullfile(fullfile(baseDir,imagingDir,subj_name{sn},[subj_name{sn} '_mean_bold.nii']));
        for i=1:length(tissues)
            P{i+1} = fullfile(baseDir, anatomicalDir, subj_name{sn}, sprintf('%s_T1w_c%d.nii', subj_name{sn}, tissues(i)));

        end
        out =  fullfile(fullfile(baseDir,imagingDir,subj_name{sn},'brain_mask.nii'));
        spm_imcalc(char(P),out,'i1>0 & (i2+i3+i4)>0.7');
        
        out =  fullfile(fullfile(baseDir,imagingDir,subj_name{sn},'gray_mask.nii'));
        spm_imcalc(char(P),out,'i1>0 & i2>0.4');
    
    case 'GLM:glm1'                   % FAST glm w/out hpf one regressor per task and per instruction
        % GLM with FAST and no high pass filtering
        % 'spm_get_defaults' code modified to allow for -v7.3 switch (to save>2MB FAST GLM struct)
        % EXAMPLE: bsp_imana('GLM:glm1',[1:XX],[1:XX])
        sn=varargin{1};
        runs=varargin{2}; % [1:16]
        
        announceTime=5;  % Instruction duration
        glm=1;
      %  subjs=length(sn);
        
        % load in task information
        C=dload(fullfile(baseDir,'pontine_taskConds_GLM.tsv'));
        Cc=getrow(C,C.StudyNum==1 & C.condNum==1); % Only get the first condition for each tasks
        nTask      = max(Cc.taskNum);                                 % how many tasks there are?
        
        for s=sn %IH: initially s=1:subjs
            T=[];
            A = dload(fullfile(baseDir,'data',subj_name{sn},[subj_name{sn} '_scans.tsv']));
            %A = getrow(A,A.run_num>=funcRunNum(1) & A.run_num<=funcRunNum(2)); % get only the runs we need (remove any test or Behav training)
            
            glmSubjDir =[baseDir, sprintf('/GLM_firstlevel_%d/',glm),subj_name{sn}];dircheck(glmSubjDir); % Create GLM folder
            
            % Fill up struct for glm
            J.dir = {glmSubjDir};
            J.timing.units = 'secs';
            J.timing.RT = 1.0;
            J.timing.fmri_t = 8;
            J.timing.fmri_t0 = 1;
            
            for r=1:numel(runs) % loop through runs
                P=getrow(A,A.run_num==runB(r));
                for i=1:(numTRs-numDummys) % get the filenames of the nifti volumes for each run
                    N{i} = [fullfile(baseDir,imagingDir,subj_name{sn},sprintf('run_%2.2d.nii,%d',runs(r),i+numDummys))];
                end;
                J.sess(r).scans= N; % number of scans in run
                
                % loop through tasks
                itt = 1;
                for it = 1:nTask
                    % The order of tasks are different for each run, to
                    % have a common order for the tasks, I will be reading
                    % from the Cc file for all the runs and subjects
                    ST = find(strcmp(P.task_name,Cc.taskNames{it}));
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
            
            dsave(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');  % Phase out the use of MAT files - prefer tsv files for Python-compatibility
            dsave(fullfile(J.dir{1},'SPM_info.tsv'),T);
            fprintf('glm_%d has been saved for %s \n',glm, subj_name{sn(s)});
        end
        
    case 'GLM:glm2'                   % FAST glm w/out hpf one regressor per task - common instruction regressor
        % GLM with FAST and no high pass filtering
        % 'spm_get_defaults' code modified to allow for -v7.3 switch (to save>2MB FAST GLM struct)
        % EXAMPLE: bsp_imana('GLM:glm1',[1:XX],[1:XX])
        sn=varargin{1};
        runs=varargin{2}; % 
        
        announceTime=5;  % Instruction duration
        glm=2;
        subjs=length(sn);
        
        % load in task information
        C=dload(fullfile(baseDir,'pontine_taskConds_GLM_reordered.tsv'));
        Cc=getrow(C,C.StudyNum==1 & C.condNum==1); % Only get the first condition for each tasks
        nTask      = max(Cc.taskNum);                                 % how many tasks there are?
        
        for s=1:subjs,
            T=[];
            A = dload(fullfile(baseDir,'data',subj_name{sn},[subj_name{sn} '_scans.tsv'])); % get scanning timing and order
            %A = getrow(A,A.run_num>=funcRunNum(1) & A.run_num<=funcRunNum(2)); % get only the runs we need (remove any test or Behav training)
            
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
                    N{i} = [fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('run_%2.2d.nii,%d',runs(r),i+numDummys))];
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
                
                % filling in the fields for SPM.mat for the instruction
              
            
                J.sess(r).cond(1).name     = 'Instruct';
                J.sess(r).cond(1).onset    = instruct_onset; % correct start time for numDummys and announcetime included (not for instruct)
                J.sess(r).cond(1).duration = 5;              % instructions last for 5 sec
                J.sess(r).cond(1).tmod     = 0; %parametric modulation (1 means yes, 0 means no)
                J.sess(r).cond(1).orth     = 0; %orthogonalized with respect to other conditions (1 means yes, 0 means no)
                J.sess(r).cond(1).pmod     = struct('name', {}, 'param', {}, 'poly', {}); %parametric structure name, modulation, parameter, polynomial expansion
                
                % loop through tasks
                for it = 1:nTask
                    % The order of tasks are different for each run, to
                    % have a common order for the tasks, I will be reading
                    % from the Cc file for all the runs and subjects
                    %the following line identifies the index of each
                    %task in Cc.taskNames in P.task_name. (Ex: Cc.taskNames
                    %has visual search as the first task; in P.task_name,
                    %visual search is the third task. Thus, ST = 3)
                    %This loop creates a regressor for each task 
                   
                    ST = find(strcmp(lower(P.task_name),lower(Cc.taskNames{it})));
                    
                    % get the task onset (instruction onset + announceTime)
                    onset = P.real_start_time(ST) - J.timing.RT*numDummys +announceTime;
                    
                    % filling in the fields for SPM_info.mat
                    S.task      = it;
                    S.taskName  = Cc.taskNames{it};
                    S.inst      = 0;
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
            
            % save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T'); % Phase out the use of MAT files - prefer tsv files for Python-compatibility
            dsave(fullfile(J.dir{1},'SPM_info.tsv'),T); 
            fprintf('glm_%d has been saved for %s \n',glm, subj_name{sn(s)});
        end

    case 'GLM:glm3'                   % FAST glm w/out hpf one regressor per task - common instruction regressor
        % GLM with FAST and no high pass filtering
        % 'spm_get_defaults' code modified to allow for -v7.3 switch (to save>2MB FAST GLM struct)
        % EXAMPLE: bsp_imana('GLM:glm1',[1:XX],[1:XX])

        sn=varargin{1};
        runs=varargin{2}; % 
        
        announceTime=5;  % Instruction duration
        glm=3;
        subjs=length(sn);
        
        % load in task information
        C=dload(fullfile(baseDir,'pontine_taskConds_GLM_reordered.tsv'));
        Cc=getrow(C,C.StudyNum==1 & C.condNum==1); % Only get the first condition for each tasks
        nTask      = max(Cc.taskNum);                               
        
        for s=1:subjs,
            T=[];

             A = dload(fullfile(baseDir,'data',subj_name{sn},[subj_name{sn} '_scans.tsv'])); % get scanning timing and order
            %A = getrow(A,A.run_num>=funcRunNum(1) & A.run_num<=funcRunNum(2)); % get only the runs we need (remove any test or Behav training)
            
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
                    N{i} = [fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('run_%2.2d.nii,%d',runs(r),i+numDummys))];
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
                
                % filling in the fields for SPM.mat for the instruction
              
            
                J.sess(r).cond(1).name     = 'Instruct';
                J.sess(r).cond(1).onset    = instruct_onset; % correct start time for numDummys and announcetime included (not for instruct)
                J.sess(r).cond(1).duration = 5;              % instructions last for 5 sec
                J.sess(r).cond(1).tmod     = 0; %parametric modulation (1 means yes, 0 means no)
                J.sess(r).cond(1).orth     = 0; %orthogonalized with respect to other conditions (1 means yes, 0 means no)
                J.sess(r).cond(1).pmod     = struct('name', {}, 'param', {}, 'poly', {}); %parametric structure name, modulation, parameter, polynomial expansion
                
                % loop through tasks
                for it = 1:nTask
                    % The order of tasks are different for each run, to
                    % have a common order for the tasks, I will be reading
                    % from the Cc file for all the runs and subjects
                    %the following line identifies the index of each
                    %task in Cc.taskNames in P.task_name. (Ex: Cc.taskNames
                    %has visual search as the first task; in P.task_name,
                    %visual search is the third task. Thus, ST = 3)
                    %This loop creates a regressor for each task 
                   
                    ST = find(strcmp(lower(P.task_name),lower(Cc.taskNames{it})));
                    
                    % get the task onset (instruction onset + announceTime)
                    onset = P.real_start_time(ST) - J.timing.RT*numDummys +announceTime;
                    
                    % filling in the fields for SPM_info.mat
                    S.task      = it;
                    S.taskName  = Cc.taskNames{it};
                    S.inst      = 0;
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

                % Load require Physio files 
                % struct('name', {'HR1','HR2'}, 'val', {[1,2,2],[21,2,1]})
                % Take 1-4 from Retro_HR and 1 from HR 

                J.sess(r).multi = {''};

                % Load HR file 
                A2 = load(fullfile(baseDir,'physio',subj_name{sn},sprintf('run%02d',r),'reg_hr.txt'));
                J.sess(r).regress(1).name = 'HR'; 
                J.sess(r).regress(1).val  = A2(1:(numTRs-numDummys),1); 
                
                % Load RetroIcor
                A3 = load(fullfile(baseDir,'physio',subj_name{sn},sprintf('run%02d',r),'reg_retro_hr.txt'));
                reg_name = {'sin1','cos1','sin2','cos2'};
                for i = [1:4]
                    J.sess(r).regress(1+i).name = reg_name{i}; 
                    J.sess(r).regress(1+i).val  = A3(1:1:(numTRs-numDummys),i); 
                end;
                                    % filling in the fields for SPM_info.mat
                S.task      = 100;
                S.taskName  = 'HR';
                S.inst      = 0;
                S.time      = onset;
                T  = addstruct(T, S);

                S.task = 101;
                S.taskName = 'sin1';
                S.inst = 0;
                S.time = onset;
                T = addstruct(T,S);

                S.task = 102;
                S.taskName = 'cos1';
                S.inst = 0;
                S.time = onset; 
                T = addstruct(T,S);

                S.task = 103;
                S.taskName = 'sin2';
                S.inst = 0;
                S.time = onset;
                T = addstruct(T,S);

                S.task = 104; 
                S.taskName = 'cos2';
                S.inst = 0;
                S.time = onset; 
                T = addstruct(T,S);
                
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
            
            % save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T'); % Phase out the use of MAT files - prefer tsv files for Python-compatibility
            dsave(fullfile(J.dir{1},'SPM_info.tsv'),T); 
            fprintf('glm_%d has been saved for %s \n',glm, subj_name{sn(s)});

        end

    case 'GLM:estimate'               % Estimate GLM depending on subjNum & glmNum
        % example: bsp_imana('GLM:estimate',1,1)
        sn=varargin{1};
        glm=varargin{2};
        
        subjs=length(sn);
        
        for s=1:subjs,
           % glmDir = [subj_name{sn(s)}]
            glmDir =[baseDir,sprintf('/GLM_firstlevel_%d/',glm),subj_name{sn(s)}];
            load(fullfile(glmDir,'SPM.mat'));
             SPM.swd=glmDir;        
            spm_rwls_spm(SPM);
        end

    case 'GLM:contrast'               % Create Contrast images for each task - and instruction 
        %%% Calculating contrast images.

        % 'SPM_light' is created in this step (xVi is removed as it slows
        % down code for FAST GLM).
        % Example1: bsp_glm('GLM:contrast', 'sn', 3, 'glm', 1, 'type', 'task')
        
        %sn             = 19;             %% list of subjects
        %glm            = 2;             %% The glm number
        
        vararginoptions(varargin, {'sn', 'glm'})
        
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, sprintf('GLM_firstlevel_%d', glm));
        
        for s = sn
            fprintf('******************** calculating contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'))
            num_runs = 16;

           % SPM  = rmfield(SPM,'xCon'); %this removes an xCon field 
            cd(fullfile(glmDir, subj_name{s}))
            T = readtable('SPM_info.tsv', 'FileType', 'text', 'Delimiter', '\t');

            % t contrast for tasks
            ucondition = unique(T.task);
            num_tasks = length(ucondition); 
            %SPM.xCon(1:num_tasks-1) = struct(); 
            
            idx = 1;
            tt = 1;

           %this initializes an SPM contrast 
            con  = zeros(1,size(SPM.xX.X,2));
                indx = find(T.task == tt & T.inst==0);
               
                con(:,indx) = 1/num_runs;  %this is the weight we give to each task. 
                con(:,T.task ~= tt & T.inst==0) = -1./(num_runs*(num_tasks-1)); % this is the mean of all other tasks (not tt)  
               
                name = T.taskName{indx(1)};
            
                SPM.xCon = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);

                %this fills the SPM contrast 

            for tt = 1:num_tasks-1 % 0 is "instruct" regressor
                con  = zeros(1,size(SPM.xX.X,2));
                indx = find(T.task == tt & T.inst==0);
                %indx_right = find(T.task == tt & T.inst==0 & T.hand ==1); %finds all the instances of a single task
                %indx_left = find(T.task == tt & T.inst==0 & T.hand ==2);
                con(:,indx) = 1/num_runs;  %this is the weight we give to each task. 
                con(:,T.task ~= tt & T.inst==0) = -1./(num_runs*(num_tasks-1)); % this is the mean of all other tasks (not tt)  
                %con(:,indx_right) = 1/num_runs;  %this is the weight we give to each task. 
                %con(:,indx_left) = -1./(num_runs); % this is the mean of all other tasks (not tt)
                name = T.taskName{indx(1)};
                %name = T.taskName{indx_right(1)}; 
                %name = sprintf('%s_right_vs_left', T.taskName{indx_right(1)});  % Naming the contrast
                %SPM.xCon = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);
                SPM.xCon(idx) = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);
                %creates contrast vector 
                idx=idx+1;
            end 
            % Add Instruction regressor - against mean of tasks 
            con  = zeros(1,size(SPM.xX.X,2));
            indx = find(T.inst>0); 
            con(:,indx) = 1/length(indx);
            con(:,T.task >0 & T.inst==0) = -1./(num_runs*num_tasks);                
            name = 'Instruct'; 
            SPM.xCon(idx) = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs); 
            
            
            SPM = spm_contrasts(SPM,1:length(SPM.xCon));  %makes contrast images: initially this - SPM = spm_contrasts(SPM,1:length(SPM.xCon))
             
                % save('SPM.mat', 'SPM','-v7.3');
            
            % rename contrast images and spmT images
            conName = {'con','spmT'};
            for i = 1:length(SPM.xCon)
                for n = 1:numel(conName)
                    oldName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%2.4d.nii',conName{n},i));
                    %newName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%s_RH.nii',conName{n},SPM.xCon(i).name));
                    newName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%s_16_runs.nii', conName{n}, SPM.xCon(i).name));
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
                case 'task' %Fcontrast for any between-task difference
                    numTasks = max(T.task);
                    con = zeros(numTasks,size(SPM.xX.X,2));
                    for i=1:numTasks
                        con(i,T.task==i & T.inst==0)=1-1/numTasks;
                        con(i,T.task~=i & T.inst==0)=-1/numTasks;
                    end

                case 'retroicor'
                    con = zeros(4,size(SPM.xX.X,2));
                    task=[101,102,103,104];

                    for i=1:4
                        con(i,T.task==task(i))=1;
                        con(i,T.task~=task(i))=0;
                    end
                           
                case 'HR'
                    con = zeros(1,size(SPM.xX.X,2));
                    con(1,T.task==100)=1;
                    con(1,T.task~=100)=0;

                            
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
        %   D.sn: Subject number (real)
        %   D.roi: Cell array of ROI names
        %   D.model: Specific combination in inX and inK regressors
        %   D.method: Regression method (reg)
        %   D.run: Run that served as test set
        %   USING THE ENTIRE MODEL (TASKS, INSTRUCTION, NUISANCE)
        %   THIS IS AFTER ANYTHING IN INK IS REMOVED - I.E. THE DATA IS
        %   FILTERED AND WEIGHTED (W)
        %   D.R: Predictive R of the whole model (against observed ts) **
        %   D.R2: Predictive R2 of the whole model
        %   USING TASK RELATED REGRESSORS ONLY: 
        %   D.R_Y: Predictived R on predicted vs. observed time series
        %   D.R_B: Predictive R on beta estimates (not centered)
        %   D.R_Bc: Predictive R on beta estimates (centered) **
        %   D.R_Yp: Predictive R on predicted vs. fitted time series (nc)
        %   D.R_Ypc: Predictive R on predicted vs. fitted time series (centered)
        sn = 3;
        glm = 1;
        roi = {'dentate'};%{'cerebellum_gray','dentate','brainstem','pontine'};
        inK = {};                   % Terms in filtering matrix - except for intercept
        inX = {{'Tasks','InstructC'},{'Tasks','InstructC'}}; % Terms in the design matrix
        reg = {'OLS','GLS'};  % Regression methods
        runs = [1:16];
        D = []; % Result structure
        % Get the posible options to test
        vararginoptions(varargin,{'sn','roi','inK','inX','reg','runs'});
        
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
                load(fullfile(baseDir,simDir,'data',subj_name{s},sprintf('rawts_%s.mat',roi{r})));
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
                            X0 = [X0 bsp_get_feature(inK{model}{t},s,SPM,INFO,'separate',1,'sscale',1)];
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
                        x = bsp_get_feature(inX{model}{t},s,SPM,INFO,'separate',0,'sscale',1);
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
                                R         = eye(N)-W*X0*pinv(W*X0);
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
                            % Record performance
                            T.roi = roi(r);
                            T.run  = rn;
                            T.sn   = s;
                            T.subjID = subj_name{s};
                            T.method  = reg(method);
                            T.methodN = method;
                            T.model = model;
                            T.inX = {inX{model}};  %%%%%%%%%%%%%%%%%%%%%%
                            T.modelL = length(inX{model});  %%%%%%%%%%%%%%
                            T.theta = nan(1,7);
                            T.theta(1:length(theta)) = theta;
                            T.time = time;
                            % Evaluation of the overall model: against observed time series  
                            T.R        = calc_cosang(Xr(testI,:)*Btrain,Yr(testI,:)); 
                            T.R2       = calc_R2(Xr(testI,:)*Btrain,Yr(testI,:));
                            % Evalation of the first set of regressor alone
                            % alone 
                            evindx = find(group==1);
                            q = length(evindx); 
                            C = eye(q)-ones(q)/q; 
                            Btr = Btrain(evindx,:); % Task-related beta weights from training set 
                            Bte = Btest(evindx,:);  % Task-related beta weights from test run 
                            Xte = Xr(testI,evindx); % Design martrix for the test run (task related regressors) 
                            T.R_Y      = calc_cosang(Xte*Btr,Yr(testI,:)); % R of between predicted and observed time series
                            T.R_B      = calc_cosang(Btr,Bte);             % R of between predicted and observed task betas (non-centered)                            
                            T.R_Bc     = calc_cosang(C*Btr,C*Bte);         % R of between predicted and observed task betas (centered)
                            T.R_Yp     = calc_cosang(Xte*Btr,Xte*Bte);     % R of between predicted vs. fitted time series (non-centered)                            
                            T.R_Ypc    = calc_cosang(Xte*C*Btr,Xte*C*Bte); % R of between predicted vs. fitted time series (centered)
                            D = addstruct(D,T);
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
        vararginoptions(varargin,{'roi','inX','inK','sn','reg','evalX','runs'});
        
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
                            X0 = [X0 bsp_get_feature(inK{model}{t},s,SPM,INFO,'separate',1,'sscale',1)];
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
                        x = bsp_get_feature(inX{model}{t},s,SPM,INFO,'separate',0,'sscale',1);
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
                                    T.subjID = subj_name{s};
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
    
    case 'test_GLM_lowfreq'
        % Compare different methods to deal with auto-correlated noise
        % And sporardic artifacts
        sn = [1:8];
        model = {{'Tasks','Instruct'},...
            {'Tasks','Instruct'}};
        inK   = {{'Hpass'},...
            {}};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'OLS','GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_lowfreq.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_lowfreq'
        what = 'R_Bc'; % what to plot - here correlation on 
        sn = [5:8]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_lowfreq.mat');
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1]}; 
        style={':',':','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'HpassOLS','OLS','HpassGLS','GLS'},'facecolor',color); 
            title('Different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
            subplot(num_subj,2,(s-1)*2+2); 
            lineplot(D.run,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'HpassOLS','OLS','HpassGLS','GLS'},'linecolor',color,...
                    'linestyle',style,'linewidth',2); % {'HpassOLS','HpassGLS','OLS','GLS'} 
            title('Different Runs across ROIs');
        end; 
        
    case 'test_GLM_Physio_full_model'
        sn = [1:8];
        model = {{'Tasks','InstructC'},...
            {'Retro_HR'},...
            {'Retro_RESP'},...
            {'HR'},...
            {'RV'},...
            {'Tasks','InstructC','Retro_HR'},...
            {'Tasks','InstructC','Retro_RESP'},...
            {'Tasks','InstructC','HR'},...
            {'Tasks','InstructC','RV'},...
            {'Retro_HR','Retro_RESP'},...
            {'Retro_HR','HR'},...
            {'Retro_HR','RV'},...
            {'Retro_RESP','HR'},...
            {'Retro_RESP','RV'},...
            {'HR','RV'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP'},...
            {'Tasks','InstructC','Retro_HR','HR'},...
            {'Tasks','InstructC','Retro_HR','RV'},...
            {'Tasks','InstructC','Retro_RESP','HR'},...
            {'Tasks','InstructC','Retro_RESP','RV'},...
            {'Tasks','InstructC','HR','RV'},...
            {'Retro_HR','Retro_RESP','HR'},...
            {'Retro_HR','Retro_RESP','RV'},...
            {'Retro_HR','HR','RV'},...
            {'Retro_RESP','HR','RV'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP','HR'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP','RV'},...
            {'Tasks','InstructC','Retro_HR','HR','RV'},...
            {'Tasks','InstructC','Retro_RESP','HR','RV'},...
            {'Retro_HR','Retro_RESP','HR','RV'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP','HR','RV'}
            };
        inK   = {{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}};
%         inK     = {{}};
%         roi = {'pontine','dentate','olive','csf','cerebellum_gray'};

        roi = {'simulate_snr35-12'};
        method = {'tikhonov_pcm'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_physio_simulate_snr35-12_tikhonov.mat'),'-struct','D');
        varargout={D};
        
    case 'test_GLM_Physio_full_model_task'
        sn = [1:4];
        model = {{'Tasks','InstructC'},...
%             {'Tasks','InstructC','Retro_HR'},...
%             {'Tasks','InstructC','Retro_HR','Retro_RESP'},...
%             {'Tasks','InstructC','Retro_HR','HR'},...
%             {'Tasks','InstructC','Retro_HR','RV'},...
%             {'Tasks','InstructC','Retro_HR','Retro_RESP','HR'},...
%             {'Tasks','InstructC','Retro_HR','Retro_RESP','RV'},...
%             {'Tasks','InstructC','Retro_HR','HR','RV'},...
%             {'Tasks','InstructC','Retro_RESP'},...
%             {'Tasks','InstructC','Retro_RESP','HR'},...
%             {'Tasks','InstructC','Retro_RESP','RV'},...
%             {'Tasks','InstructC','Retro_RESP','HR','RV'},...
%             {'Tasks','InstructC','HR'},...
%             {'Tasks','InstructC','HR','RV'},...
%             {'Tasks','InstructC','RV'},.....
            {'Tasks','InstructC','Retro_HR','Retro_RESP','HR','RV'}
            };
        inK   = {{},{}};
%         inK   = {{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}};
%         roi = {'simulate_signal35','simulate_snr35-3','simulate_snr35-12','simulate_signal0'};
        roi = {'cerebellum_gray','csf'};
        method = {'OLS','GLS','tikhonov_pcm'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_physio_task_model_cerebellum.mat'),'-struct','D');
        roi = {'cerebellum_gray'};
        method = {'OLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_physio_allmodels_tikhonov_S98.mat'),'-struct','D');
        varargout={D};
        
    case 'test_GLM_Physio_full_model_task_simulate'
        sn = [1];
        model = {{'Null'},...
            {'Tasks','InstructC'},...
            {'Retro_HR'},...
            {'Retro_RESP'},...
            {'HR'},...
            {'RV'},...
            {'Tasks','InstructC','Retro_HR'},...
            {'Tasks','InstructC','Retro_RESP'},...
            {'Tasks','InstructC','HR'},...
            {'Tasks','InstructC','RV'},...
            {'Retro_HR','Retro_RESP'},...
            {'Retro_HR','HR'},...
            {'Retro_HR','RV'},...
            {'Retro_RESP','HR'},...
            {'Retro_RESP','RV'},...
            {'HR','RV'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP'},...
            {'Tasks','InstructC','Retro_HR','HR'},...
            {'Tasks','InstructC','Retro_HR','RV'},...
            {'Tasks','InstructC','Retro_RESP','HR'},...
            {'Tasks','InstructC','Retro_RESP','RV'},...
            {'Tasks','InstructC','HR','RV'},...   %%%%%%
            {'Retro_HR','Retro_RESP','HR'},...
            {'Retro_HR','Retro_RESP','RV'},...
            {'Retro_HR','HR','RV'},...
            {'Retro_RESP','HR','RV'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP','HR'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP','RV'},...
            {'Tasks','InstructC','Retro_HR','HR','RV'},...
            {'Tasks','InstructC','Retro_RESP','HR','RV'},...
            {'Retro_HR','Retro_RESP','HR','RV'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP','HR','RV'}
            };
%         inK   = {{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}};
        inK = {{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}};
%         roiName = {'regress0.0a', 'regress0.5a', 'regress1.0a', 'regress1.5a'};
        roiName = {'regress2.0x2.0d'};
        k = 1;
%         k = (61:1:100);
%         k = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 29 30 31 32 33 34 35 36 37 38 39 40 42 43 44 45 46 47 48 49 50 51 52];
%         for i = 1:40
%             roi{i} = sprintf('simulate_%s_%04d',roiName{1},k(i));
%         end
        roi = {'simulate_regress2.0x2.0d_0001'};
        method = {'OLS','GLS','tikhonov_pcm'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','simulate_GLM_physio_task_allmodels_50sims_regress2.0x2.0d_S98.mat'),'-struct','D');
        varargout={D};
        
    case 'test_GLM_Physio_full_model_task_simulate_test'
        sn = [1];
        model = {
%             {'Tasks','InstructC','Retro_HR','Retro_RESP'},...
            {'Tasks','InstructC','HR','RV'},...   %%%%%%
            };
%         inK   = {{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}};
        inK = {{},{}};
%         roiName = {'regress0.0a', 'regress0.5a', 'regress1.0a', 'regress1.5a'};
        roiName = {'regress2.0x2.0d'};
        k = 1;
%         k = (61:1:100);
%         k = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 29 30 31 32 33 34 35 36 37 38 39 40 42 43 44 45 46 47 48 49 50 51 52];
%         for i = 1:40
%             roi{i} = sprintf('simulate_%s_%04d',roiName{1},k(i));
%         end
        roi = {'simulate_regress2.0x2.0d_0001'};
        method = {'tikhonov_pcm'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','simulate_GLM_physio_task_allmodels_50sims_regress2.0x2.0d_S98_test.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_Physio_full_model_task'
        what = 'R_Bc'; % what to plot - here correlation on 
        sn = [2 4 6 8]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_physio_task_model_simulate.mat');
         
        num_subj = length(sn);
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[1 0.7 0.7],[0.7 0.7 1]};
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,1,(s-1)*1+1);
            barplot(D.roi,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'OLS-Task','OLS-All','GLS-Task','GLS-All','Tikhonov-Task','Tikhonov-All'},'facecolor',color); 
            title('R-Bc for different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
        end; 
        
    case 'test_GLM_Physio_full_model_task_csfbrainstem'
        sn = [1:8];
        model = {{'Tasks','InstructC'},...
            {'Tasks','InstructC','Retro_HR'},...
            {'Tasks','InstructC','CSFmedulla'},...
            {'Tasks','InstructC','CSFpostdrain'},...
            {'Tasks','InstructC','Mov'},.....
            {'Tasks','InstructC','Retro_HR','CSFmedulla'},...
            {'Tasks','InstructC','Retro_HR','CSFpostdrain'},...
            {'Tasks','InstructC','Retro_HR','Mov'},...
            {'Tasks','InstructC','CSFmedulla','CSFpostdrain'},...
            {'Tasks','InstructC','CSFmedulla','Mov'},...
            {'Tasks','InstructC','CSFpostdrain','Mov'},...
            {'Tasks','InstructC','Retro_HR','CSFmedulla','CSFpostdrain'},...
            {'Tasks','InstructC','Retro_HR','CSFmedulla','Mov'},...
            {'Tasks','InstructC','CSFmedulla','CSFpostdrain','Mov'},...
            {'Tasks','InstructC','Retro_HR','CSFmedulla','CSFpostdrain','Mov'}
            };
        inK   = {{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_physio_full_model_task_csfbrainstem.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_Physio_full_model_task_csfbrainstem'
        what = 'R'; % what to plot - here correlation on 
        sn = [1:4]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_physio_full_model_task_csfbrainstem.mat');
         
        num_subj = length(sn);
        color={[0.7 0 0],[0 0 0.7],[0 0.7 0],[1 0.4 0.4],[0.4 0.4 1],[0.4 1 0.4],[0.8 0 0],[0 0 0.8],[0 0.8 0],[1 0.5 0.5],[0.5 0.5 1],[0.5 1 0.5],[0.9 0 0],[0 0 0.9],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,1,(s-1)*1+1);
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'Task','Task+RetHR','Task+CSFmedulla','Task+CSFpostdrain','Task+Mov','Task+RetHR+CSFmedulla','Task+RetHR+CSFpostdrain','Task+RetHR+Mov','Task+CSFmedulla+CSFpostdrain','Task+CSFmedulla+Mov','Task+CSFpostdrain+Mov','Task+RetHR+CSFmedulla+CSFpostdrain','Task+RetHR+CSFmedulla+Mov','Task+CSFmedulla+CSFpostdrain+Mov','Task+RetHR+CSFmedulla+CSFpostdrain+Mov'},'facecolor',color); 
            title('R for different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
        end;
       
    case 'test_GLM_Physio_full_model_csf'
        sn = [1:8];
        model = {{'Tasks','InstructC'},...
            {'Retro_HR'},...
            {'Retro_RESP'},...
            {'HR'},...
            {'RV'},...
            {'Tasks','InstructC','Retro_HR'},...
            {'Tasks','InstructC','Retro_RESP'},...
            {'Tasks','InstructC','HR'},...
            {'Tasks','InstructC','RV'},...
            {'Retro_HR','Retro_RESP'},...
            {'Retro_HR','HR'},...
            {'Retro_HR','RV'},...
            {'Retro_RESP','HR'},...
            {'Retro_RESP','RV'},...
            {'HR','RV'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP'},...
            {'Tasks','InstructC','Retro_HR','HR'},...
            {'Tasks','InstructC','Retro_HR','RV'},...
            {'Tasks','InstructC','Retro_RESP','HR'},...
            {'Tasks','InstructC','Retro_RESP','RV'},...
            {'Tasks','InstructC','HR','RV'},...
            {'Retro_HR','Retro_RESP','HR'},...
            {'Retro_HR','Retro_RESP','RV'},...
            {'Retro_RESP','HR','RV'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP','HR'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP','RV'},...
            {'Tasks','InstructC','Retro_RESP','HR','RV'},...
            {'Retro_HR','Retro_RESP','HR','RV'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP','HR','RV'}
            };
        inK   = {{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}};
        roi = {'galenic','medulla','midbrain','pons','postdrain','transverseL','transverseR','ventricle4'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_physio_full_model_csf_gmmask.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_Physio_full_model_csf'
        what = 'R'; % what to plot - here correlation on 
        sn = [1:4]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_physio_full_model_csf.mat');
         
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'Retro HR','Retro Resp','HR','RV','none'},'facecolor',color); 
            title('different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
            %subplot(num_subj,2,(s-1)*2+2); 
            %lineplot(D.run,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'Retro HR','Retro Resp','HR','RV','none'},'linecolor',color,...
            %        'linestyle',style,'linewidth',2);
            %title('Different Runs across ROIs');
        end; 

    case 'test_GLM_Physio_full_model_phys_csf'
        sn = [1:8];
        model = {{'Retro_HR'},...
            {'Retro_RESP'},...
            {'HR'},...
            {'RV'},...
            {'Retro_HR','Retro_RESP'},...
            {'Retro_HR','HR'},...
            {'Retro_HR','RV'},...
            {'Retro_RESP','HR'},...
            {'Retro_RESP','RV'},...
            {'HR','RV'},...
            {'Retro_HR','Retro_RESP','HR'},...
            {'Retro_HR','Retro_RESP','RV'},...
            {'Retro_RESP','HR','RV'},...
            {'Retro_HR','Retro_RESP','HR','RV'},...
            };
        inK   = {{},{},{},{},{},{},{},{},{},{},{},{},{},{}};
        roi = {'galenic','medulla','midbrain','pons','postdrain','transverseL','transverseR','ventricle4'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_physio_full_model_phys_csf.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_Physio_full_model_phys_csf'
        what = 'R'; % what to plot - here correlation on 
        sn = [1:4]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_physio_full_model_phys_csf.mat');
         
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,1,(s-1)*1+1);
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'Retro HR','Retro Resp','HR','RV','none'},'facecolor',color); 
            title('different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
        end;
        
    case 'test_GLM_Physio_task_model_csf'
        sn = 1:8;
        model = {{'Tasks','InstructC'},...
            {'Tasks','InstructC','Retro_HR'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP'},...
            {'Tasks','InstructC','Retro_HR','HR'},...
            {'Tasks','InstructC','Retro_HR','RV'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP','HR'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP','RV'},...
            {'Tasks','InstructC','Retro_HR','HR','RV'},...
            {'Tasks','InstructC','Retro_RESP'},...
            {'Tasks','InstructC','Retro_RESP','HR'},...
            {'Tasks','InstructC','Retro_RESP','RV'},...
            {'Tasks','InstructC','Retro_RESP','HR','RV'},...
            {'Tasks','InstructC','HR'},...
            {'Tasks','InstructC','HR','RV'},...
            {'Tasks','InstructC','RV'},...
            {'Tasks','InstructC','Retro_HR','Retro_RESP','HR','RV'}
            };
        inK   = {{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}};
        roi = {'galenic','medulla','midbrain','pons','postdrain','transverseL','transverseR','ventricle4'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_physio_task_model_csf_gmmask_new.mat'),'-struct','D');
        varargout={D};
        
    case 'plot_GLM_Physio_task_model_csf'
        what = 'R'; % what to plot - here correlation on 
        sn = [5:8]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_physio_task_model_csf_gmmask_new.mat');
         
        num_subj = length(sn);
        color={[0.8 0.8 0.8],[0.9 1 1],[0.7 0.8 0.9],[0.5 0.8 1],[0 0.7 1],[0 0 1],[0 0 0.5],[0.1 0.1 0.4],[0.6 0.9 0.6],[0 1 0.5],[0 0.5 0],[0 0.4 0],[1 1 0.9],[1 1 0],[1 0 0],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,1,(s-1)*1+1);
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'Task','Task+RetHR','Task+RetHR+RetRESP','Task+RetHR+HR','Task+RetHR+RV','Task+RetHR+RetRESP+HR','Task+RetHR+RetRESP+RV','Task+RetHR+HR+RV','Task+RetRESP','Task+RetRESP+HR','Task+RetRESP+RV','Task+RetRESP+HR+RV','Task+HR','Task+HR+RV','Task+RV','Task+RetHR+RetRESP+HR+RV'},'facecolor',color); 
            title('R for different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
        end;
        
    case 'test_GLM_Physio_task_motion_csf'
        sn = [1:8];
        model = {{'Tasks','InstructC'},...
            {'Tasks','InstructC','Retro_HR'},...
            {'Tasks','InstructC','Mov'},...
            {'Tasks','InstructC','MovPCA'},...
            {'Tasks','InstructC','Retro_HR','Mov'},...
            {'Tasks','InstructC','Retro_HR','MovPCA'},...
            {'Tasks','InstructC','Mov','MovPCA'},...
            {'Tasks','InstructC','Retro_HR','Mov','MovPCA'}
            };
        inK   = {{},{},{},{},{},{},{},{}};
        roi = {'galenic','medulla','midbrain','pons','postdrain','transverseL','transverseR','ventricle4'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_physio_task_motion_csf_gmmask.mat'),'-struct','D');
        varargout={D};
        
    case 'plot_GLM_Physio_task_motion_csf'
        what = 'R'; % what to plot - here correlation on 
        sn = [5:8]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_physio_task_motion_csf_gmmask.mat');
         
        num_subj = length(sn);
        color={[0.7 0 0],[0 0 0.7],[0 0.7 0],[1 0.4 0.4],[1 0.5 0.5],[0.5 0.5 1],[0.5 1 0.5],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,1,(s-1)*1+1);
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'Task','Task+RetHR','Task+Mov','Task+MovPCA','Task+RetHR+Mov','Task+RetHR+MovPCA','Task+Mov+MovPCA','Task+RetHR+Mov+MovPCA'},'facecolor',color); 
            title('R for different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
        end;
        
    case 'test_GLM_Physio_task_acompcor_csf'
        sn = [1:8];
        model = {{'Tasks','InstructC'},...
            {'Tasks','InstructC','Retro_HR'},...
            {'Tasks','InstructC','aCompCor'},...
            {'Tasks','InstructC','Retro_HR','aCompCor'}
            };
        inK   = {{},{},{},{}};
        roi = {'galenic','medulla','midbrain','pons','postdrain','transverseL','transverseR','ventricle4'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_physio_task_acompcor_csf_gmmask.mat'),'-struct','D');
        varargout={D};
        
    case 'plot_GLM_Physio_task_acompcor_csf'
        what = 'R'; % what to plot - here correlation on 
        sn = [5:8]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_physio_task_acompcor_csf_gmmask.mat');
         
        num_subj = length(sn);
        color={[0.7 0 0],[0 0.7 0],[1 0.4 0.4],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,1,(s-1)*1+1);
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'Task','Task+RetHR','Task+aCompCor','Task+RetHR+aCompCor'},'facecolor',color); 
            title('R for different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
        end;
        
    case 'test_GLM_Physio_task_csf_csf'
        sn = [1:8];
        model = {{'Tasks','InstructC'},...
            {'Tasks','InstructC','CSFmedulla'},...
            {'Tasks','InstructC','CSFmedullaPCAindiv'},...
            {'Tasks','InstructC','CSFmedullaPCAall'},...
            {'Tasks','InstructC','CSFpostdrain'},...
            {'Tasks','InstructC','CSFpostdrainPCAindiv'},...
            {'Tasks','InstructC','CSFpostdrainPCAall'}
            };
        inK   = {{},{},{},{},{},{},{}};
        roi = {'galenic','medulla','midbrain','pons','postdrain','transverseL','transverseR','ventricle4'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_physio_task_csf_csf_gmmask.mat'),'-struct','D');
        varargout={D};
        
    case 'plot_GLM_Physio_task_csf_csf'
        what = 'R'; % what to plot - here correlation on 
        sn = [5:8]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_physio_task_csf_csf_gmmask.mat');
         
        num_subj = length(sn);
        color={[0.7 0 0],[0 0 0.7],[0 0.7 0],[1 0.5 0.5],[0.5 0.5 1],[0.5 1 0.5],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,1,(s-1)*1+1);
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'Task','Task+CSFMedulla','Task+CSFMedullaPCAindiv','Task+CSFMedullaPCAall','Task+CSFpostdrain','Task+CSFpostdrainPCAindiv','Task+CSFpostdrainPCAall'},'facecolor',color); 
            title('R for different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
        end;
        
    case 'test_GLM_lowfreq_csf'
        % Compare different methods to deal with auto-correlated noise
        % And sporardic artifacts
        sn = [1:8];
        model = {{'Tasks','Instruct'},...
            {'Tasks','Instruct'}};
        inK   = {{'Hpass'},...
            {}};
        roi = {'galenic','medulla','midbrain','pons','postdrain','transverseL','transverseR','ventricle4'};
        method = {'OLS','GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_lowfreq_csf.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_lowfreq_csf'
        what = 'R_Bc'; % what to plot - here correlation on 
        sn = [1:4]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_lowfreq_csf.mat');
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1]}; 
        style={':',':','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'HpassOLS','OLS','HpassGLS','GLS'},'facecolor',color); 
            title('Different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
            subplot(num_subj,2,(s-1)*2+2); 
            lineplot(D.run,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'HpassOLS','OLS','HpassGLS','GLS'},'linecolor',color,...
                    'linestyle',style,'linewidth',2); % {'HpassOLS','HpassGLS','OLS','GLS'} 
            title('Different Runs across ROIs');
        end;
        
    case 'test_GLM_lowfreq_physio_filter_csf'
        % Compare different methods to deal with auto-correlated noise
        % And sporardic artifacts
        sn = [1:8];
        model = {{'Tasks','Instruct'},...
            {'Tasks','Instruct'},...
            {'Tasks','Instruct'},...
            {'Tasks','Instruct'},...
            {'Tasks','Instruct'}};
        inK   = {{'Retro_HR'},...
                 {'Retro_RESP'},...
                 {'HR'},...
                 {'RV'},...
            {}};
        roi = {'galenic','medulla','midbrain','pons','postdrain','transverseL','transverseR','ventricle4'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_lowfreq_physio_filter_csf.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_lowfreq_physio_filter_csf'
        what = 'R_Bc'; % what to plot - here correlation on 
        sn = [1:4]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_lowfreq_physio_filter_csf.mat');
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1]}; 
        style={':',':','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'Retro_HR','Retro_RESP','HR','RV','None'},'facecolor',color); 
            title('Different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
            subplot(num_subj,2,(s-1)*2+2); 
            lineplot(D.run,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'Retro_HR','Retro_RESP','HR','RV','None'},'linecolor',color,...
                    'linestyle',style,'linewidth',2); % {'HpassOLS','HpassGLS','OLS','GLS'} 
            title('Different Runs across ROIs');
        end;
        
    case 'test_GLM_Physio'
        sn = [1:8];
        model = {{'Tasks','InstructC','Retro_HR'},...
            {'Tasks','InstructC','Retro_RESP'},...
            {'Tasks','InstructC','HR'},...
            {'Tasks','InstructC','RV'},...
            {'Tasks','InstructC'}};
        inK   = {{},{},{},{},{}};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_physio.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_Physio'
        what = 'R_Bc'; % what to plot - here correlation on 
        sn = [5:8]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_physio.mat');
         
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'Retro HR','Retro Resp','HR','RV','none'},'facecolor',color); 
            title('different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
            subplot(num_subj,2,(s-1)*2+2); 
            lineplot(D.run,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'Retro HR','Retro Resp','HR','RV','none'},'linecolor',color,...
                    'linestyle',style,'linewidth',2);
            title('Different Runs across ROIs');
        end; 
      
    case 'test_GLM_Mov'
        sn = [1:8];
        model = {{'Tasks','InstructC','Mov'},...
            {'Tasks','InstructC','MovPCA'},...
            {'Tasks','InstructC'}};
        inK   = {{},{},{}};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_mov.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_Mov'
        what = 'R_Bc'; % what to plot - here correlation on 
        sn = [1:4]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_mov.mat');
         
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'Mov','MovPCA','none'},'facecolor',color); 
            title('different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
            subplot(num_subj,2,(s-1)*2+2); 
            lineplot(D.run,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'Mov','MovPCA','none'},'linecolor',color,...
                    'linestyle',style,'linewidth',2);
            title('Different Runs across ROIs');
        end;     
    
    case 'test_GLM_CSF'
        sn = [1:8];
        model = {{'Tasks','InstructC','CSF'},...
            {'Tasks','InstructC','CSFPCAindiv'},...
            {'Tasks','InstructC','CSFPCAall'},...
            {'Tasks','InstructC'}};
        inK   = {{},{},{},{}};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_csf.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_CSF'
        what = 'R_Bc'; % what to plot - here correlation on 
        sn = [5:8]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_csf.mat');
         
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'CSF','CSFPCAindiv','CSFPCAall','none'},'facecolor',color); 
            title('different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
            subplot(num_subj,2,(s-1)*2+2); 
            lineplot(D.run,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'CSF','CSFPCAindiv','CSFPCAall','none'},'linecolor',color,...
                    'linestyle',style,'linewidth',2);
            title('Different Runs across ROIs');
        end; 
        
    case 'test_GLM_CSFpons'
        sn = [1:8];
        model = {{'Tasks','InstructC','CSFpons'},...
            {'Tasks','InstructC','CSFponsPCAindiv'},...
            {'Tasks','InstructC','CSFponsPCAall'},...
            {'Tasks','InstructC'}};
        inK   = {{},{},{},{}};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_csfpons.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_CSFpons'
        what = 'R_Bc'; % what to plot - here correlation on 
        sn = [5:8]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_csfpons.mat');
         
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'CSFpons','CSFponsPCAindiv','CSFponsPCAall','none'},'facecolor',color); 
            title('different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
            subplot(num_subj,2,(s-1)*2+2); 
            lineplot(D.run,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'CSFpons','CSFponsPCAindiv','CSFponsPCAall','none'},'linecolor',color,...
                    'linestyle',style,'linewidth',2);
            title('Different Runs across ROIs');
        end;  
        
    case 'test_GLM_Retro_Mov_CSF'
        sn = [1:8];
        model = {{'Tasks','InstructC','Retro_HR','Retro_RESP','Mov','CSF'},...
            {'Tasks','InstructC','Retro_HR','Mov','CSF'},...
            {'Tasks','InstructC','Retro_HR','Mov'},...
            {'Tasks','InstructC','Retro_HR'},...
            {'Tasks','InstructC'}};
        inK   = {{},{},{},{},{}};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_ret_mov_csf.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_Retro_Mov_CSF'
        what = 'R_Bc'; % what to plot - here correlation on 
        sn = [1:4]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_ret_mov_csf.mat');
         
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'RetHR+RetRESP+Mov+CSF','RetHR+Mov+CSF','RetHR+Mov','RetHR','none'},'facecolor',color); 
            title('different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
            subplot(num_subj,2,(s-1)*2+2); 
            lineplot(D.run,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'RetHR+RetRESP+Mov+CSF','RetHR+Mov+CSF','RetHR+Mov','RetHR','none'},'linecolor',color,...
                    'linestyle',style,'linewidth',2);
            title('Different Runs across ROIs');
        end; 
        
    case 'test_GLM_Retro_Mov_CSF_filter'
        sn = [1:8];
        model = {{'Tasks','InstructC'},...
            {'Tasks','InstructC'},...
            {'Tasks','InstructC'},...
            {'Tasks','InstructC'},...
            {'Tasks','InstructC'}};
        inK   = {{'Retro_HR','Retro_RESP','Mov','CSF'},{'Retro_HR','Mov','CSF'},{'Retro_HR','Mov'},{'Retro_HR'},{}};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_ret_mov_csf_filter.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_Retro_Mov_CSF_filter'
        what = 'R_Bc'; % what to plot - here correlation on 
        sn = [5:8]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_ret_mov_csf_filter.mat');
         
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'RetHR+RetRESP+Mov+CSF','RetHR+Mov+CSF','RetHR+Mov','RetHR','none'},'facecolor',color); 
            title('different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
            subplot(num_subj,2,(s-1)*2+2); 
            lineplot(D.run,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'RetHR+RetRESP+Mov+CSF','RetHR+Mov+CSF','RetHR+Mov','RetHR','none'},'linecolor',color,...
                    'linestyle',style,'linewidth',2);
            title('Different Runs across ROIs');
        end; 
        
    case 'test_GLM_notask_Physio'
        sn = [1:8];
        model = {{'Retro_HR'},...
            {'Retro_RESP'},...
            {'HR'},...
            {'RV'}};
        inK   = {{},{},{},{}};
        roi = {'galenic','medulla','midbrain','pons','postdrain','transverseL','transverseR','ventricle4'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM_notask','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_notask_physio.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_notask_Physio'
        what = 'R_Y'; % what to plot - here correlation on 
        sn = [5:8]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_notask_physio.mat');
         
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'Retro HR','Retro Resp','HR','RV'},'facecolor',color); 
            title('different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
            subplot(num_subj,2,(s-1)*2+2); 
            lineplot(D.run,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'Retro HR','Retro Resp','HR','RV'},'linecolor',color,...
                    'linestyle',style,'linewidth',2);
            title('Different Runs across ROIs');
        end; 
        
    case 'test_GLM_notask_Physio_filter'
        sn = [1:8];
        model = {{'Mov'},...
            {'Mov'},...
            {'Mov'},...
            {'Mov'}};
        inK   = {{'Retro_HR'},{'Retro_RESP'},{'HR'},{'RV'}};
        roi = {'galenic','medulla','midbrain','pons','postdrain','transverseL','transverseR','ventricle4'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM_notask','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_notask_physio_filter.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_notask_Physio_filter'
        what = 'R_Y'; % what to plot - here correlation on 
        sn = [5:8]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_notask_physio_filter.mat');
         
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'Retro HR','Retro Resp','HR','RV'},'facecolor',color); 
            title('different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
            subplot(num_subj,2,(s-1)*2+2); 
            lineplot(D.run,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'Retro HR','Retro Resp','HR','RV'},'linecolor',color,...
                    'linestyle',style,'linewidth',2);
            title('Different Runs across ROIs');
        end; 
        
    case 'test_GLM_notask_Mov'
        sn = [1:8];
        model = {{'Mov'},...
            {'MovPCA'}};
        inK   = {{},{}};
        roi = {'galenic','medulla','midbrain','pons','postdrain','transverseL','transverseR','ventricle4'};
        method = {'GLS'};
        
        D=bsp_glm('test_GLM_notask','sn',sn,'roi',roi,'inK',inK,'inX',model,'reg',method,'runs',[1:16]);
        save(fullfile(baseDir,'results','test_GLM_notask_mov.mat'),'-struct','D');
        varargout={D};
    
    case 'plot_GLM_notask_Mov'
        what = 'R_Y'; % what to plot - here correlation on 
        sn = [5:8]; 
        vararginoptions(varargin,{'what','sn'});

        D=load('test_GLM_notask_mov.mat');
         
        num_subj = length(sn); 
        color={[0.7 0 0],[0 0 0.7],[1 0.4 0.4],[0.4 0.4 1],[0.5 0.5 0.5]}; 
        style={':',':','-','-','-'}; 
        for s=1:num_subj 
            subplot(num_subj,2,(s-1)*2+1); 
            barplot(D.roi,D.(what),'split',[D.model],'subset',D.sn==sn(s),'leg',{'Mov','MovPCA'},'facecolor',color); 
            title('different ROIs');
            ylabel(sprintf('%s',subj_name{sn(s)}));
            subplot(num_subj,2,(s-1)*2+2); 
            lineplot(D.run,D.(what),'split',[D.methodN D.model],'subset',D.sn==sn(s),'leg',{'Mov','MovPCA'},'linecolor',color,...
                    'linestyle',style,'linewidth',2);
            title('Different Runs across ROIs');
        end; 

    
    case 'test_GLM_script'
        model = {{'Tasks','Instruct'},...
            {'Tasks','Instruct','Retro_HR'},...
            {'Tasks','Instruct','Retro_RESP'},...
            {'Tasks','Instruct','HR'},...
            {'Tasks','Instruct','RV'}};
        roi = {'pontine','dentate','olive','csf','cerebellum_gray'};
        method = {'OLS','GLS','ridge_pcm','tikhonov_pcm'};
        
        D=bsp_glm('test_GLM','roi',roi,'reg',method,'inX',model,'runs',[1:10]);
        save(fullfile(baseDir,'results','test_GLM_5.mat'),'-struct','D');
        varargout={D};
        
    case 'physio_reg' % Examines the covariance of physiological regressors with main regressors 
        sn  = varargin{1}; 
        glm = varargin{2}; %example: bsp_glm("physio_reg",5,2), with sn=5, glm=2
        reg = {'Tasks','InstructC','Retro_HR','Retro_RESP','HR','RV'};

        glmDir = fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm));
        glmDirSubj = fullfile(glmDir, subj_name{sn});
        load(fullfile(glmDirSubj,'SPM.mat'));
        %INFO = load(fullfile(glmDirSubj,'SPM_info.mat'));
        INFO = dload(fullfile(glmDirSubj,'SPM_info.tsv'));
        nRuns = length(SPM.nscan);
        
        % find the TR and Run for each image in the entire time series
        for i=1:nRuns 
            tr(SPM.Sess(i).row,1)=[1:SPM.nscan(i)]; 
            rn(SPM.Sess(i).row,1)=i; 
        end; 
        
        % Get the corresponding feature  (regressor)
        for r = 1:length(reg)
            X{r}=bsp_get_feature(reg{r},sn,SPM,INFO,'zscale',1); 
        end 

        % Find the onset of all instruction periors 
        %onsets = vertcat(SPM.Sess(1).U([1:2:18]).ons); 
        onsets = vertcat(SPM.Sess(1).U([1:2:10]).ons);
        
        % Run-averaged 
        to_plot=[1,2,5,6]; 
        for i = 1:4 
            subplot(4,2,(i-1)*2+1);
            Y=reshape(sum(X{to_plot(i)},2),SPM.nscan(1),16); 
            traceplot([1:SPM.nscan(1)],Y','errorfcn','stderr'); 
            drawlines(onsets,'k'); 
            set(gca,'XLim',[-3 SPM.nscan(1)+3]);
            ylabel(reg{to_plot(i)});
        end; 
        
        % Run a regression of each task regressor on Physio time series 
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

            %save plots 
            
           saveas(gcf, fullfile(baseDir,'physio',subj_name{sn},sprintf('physio_%s.png', subj_name{sn})));

        end; 

        %keyboard; 
end
end

function ca=calc_cosang(A,B)
    ca  = sum(sum(A.*B))./sqrt(sum(sum(A.*A)).*sum(sum(B.*B))); % R with timeseries 
end

function R2=calc_R2(A,B)
    R2  = 1-sum(sum((A-B).^2))./sum(sum(B.^2));
end


