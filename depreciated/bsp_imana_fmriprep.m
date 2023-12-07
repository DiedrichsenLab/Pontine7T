function varargout=bsp_imana(what,varargin)

numDummys = 3;                                                              % per run
numTRs    = 321;                                                            % per run (includes dummies)
%========================================================================================================================
% PATH DEFINITIONS
baseDir         ='/srv/diedrichsen/data/Pontine7T';
%baseDir         ='/Volumes/MotorControl/data/Pontine7T';
imagingDir      ='/imaging_data_fmriprep';
imagingDirRaw   ='/imaging_data_raw_fmriprep';
anatomicalDir   ='/anatomicals_fmriprep';
suitDir         ='/suit_fmriprep';
regDir          ='/RegionOfInterest_fmriprep';
%========================================================================================================================
% PRE-PROCESSING 
subj_name = {'p04-fmriprep'};
% loc_AC={[77;116;134];[77;116;134];[82;121;134]}; % Numbers in the titlebar of mricron
%========================================================================================================================
% GLM INFO
funcRunNum  = [50,70];  % first and last behavioural run numbers
run         = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
runB        = [50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70];  % Behavioural labelling of runs
%========================================================================================================================
switch(what)

    case 'ANAT:reslice_LPI'           % Reslice anatomical image within LPI coordinate systems  
        % example: bsp_imana('ANAT:reslice_LPI',1)
        sn  = varargin{1}; % subjNum 
        subjs=length(sn);
        
        for s=1:subjs,
            
            % (1) Reslice anatomical image to set it within LPI co-ordinate frames
            source  = fullfile(baseDir,anatomicalDir,subj_name{sn(s)},'anatomical_raw.nii');
            dest    = fullfile(baseDir,anatomicalDir,subj_name{sn(s)},'anatomical.nii');
            spmj_reslice_LPI(source,'name', dest);
            
            % (2) In the resliced image, set translation to zero
            V               = spm_vol(dest);
            dat             = spm_read_vols(V);
            V.mat(1:3,4)    = [0 0 0];
            spm_write_vol(V,dat);
            display 'Manually retrieve the location of the anterior commissure (x,y,z) before continuing'
        end
    case 'ANAT:centre_AC'             % Re-centre AC
        % Before running provide coordinates in the preprocessing section
        % example: bsp_imana('ANAT:centre_AC',1)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        for s=1:subjs,
            img    = fullfile(baseDir,anatomicalDir,subj_name{sn(s)},'rinv1_MP2RAGE.nii');
            V               = spm_vol(img);
            dat             = spm_read_vols(V);
            oldOrig         = V.mat(1:3,4);
            V.mat(1:3,4)    = oldOrig-loc_AC{sn(s)};
            spm_write_vol(V,dat);
            fprintf('Done for %s \n',subj_name{sn(s)})
        end
    case 'ANAT:segmentation'          % Segmentation + Normalisation
        % example: bsp_imana('ANAT:segmentation',1)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        
        SPMhome=fileparts(which('spm.m'));
        J=[];
        for s=1:subjs,
            J.channel.vols = {fullfile(baseDir,anatomicalDir,subj_name{sn(s)},'anatomical.nii,1')};
            J.channel.biasreg = 0.001;
            J.channel.biasfwhm = 60;
            J.channel.write = [0 0];
            J.tissue(1).tpm = {fullfile(SPMhome,'tpm/TPM.nii,1')};
            J.tissue(1).ngaus = 1;
            J.tissue(1).native = [1 0];
            J.tissue(1).warped = [0 0];
            J.tissue(2).tpm = {fullfile(SPMhome,'tpm/TPM.nii,2')};
            J.tissue(2).ngaus = 1;
            J.tissue(2).native = [1 0];
            J.tissue(2).warped = [0 0];
            J.tissue(3).tpm = {fullfile(SPMhome,'tpm/TPM.nii,3')};
            J.tissue(3).ngaus = 2;
            J.tissue(3).native = [1 0];
            J.tissue(3).warped = [0 0];
            J.tissue(4).tpm = {fullfile(SPMhome,'tpm/TPM.nii,4')};
            J.tissue(4).ngaus = 3;
            J.tissue(4).native = [1 0];
            J.tissue(4).warped = [0 0];
            J.tissue(5).tpm = {fullfile(SPMhome,'tpm/TPM.nii,5')};
            J.tissue(5).ngaus = 4;
            J.tissue(5).native = [1 0];
            J.tissue(5).warped = [0 0];
            J.tissue(6).tpm = {fullfile(SPMhome,'tpm/TPM.nii,6')};
            J.tissue(6).ngaus = 2;
            J.tissue(6).native = [0 0];
            J.tissue(6).warped = [0 0];
            J.warp.mrf = 1;
            J.warp.cleanup = 1;
            J.warp.reg = [0 0.001 0.5 0.05 0.2];
            J.warp.affreg = 'mni';
            J.warp.fwhm = 0;
            J.warp.samp = 3;
            J.warp.write = [1 1];
            matlabbatch{1}.spm.spatial.preproc=J;
            spm_jobman('run',matlabbatch);
            fprintf('Check segmentation results for %s\n', subj_name{sn(s)})
        end;
        
    case 'FUNC:remDum'                % Remove the extra dummy scans from all functional runs.
        % funtional runs have to be named as run_01-r.nii, run_02-r.nii ...
        % example: bsp_imana('FUNC:remDum',1)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        for s=1:subjs,
            cd(fullfile(baseDir,imagingDirRaw,[subj_name{sn(s)},'-n']));
            funScans = dir('run-*.nii');
            for i=1:length(funScans)  
                outfilename = sprintf('run-%2.2d.nii',i);
                mkdir temp;
                spm_file_split(funScans(i).name,'temp');
                cd temp;
                list = dir('run*.nii');
                list = list(numDummys+1:end);  % Remove dummies
                V = {list(:).name};
                spm_file_merge(V,outfilename);
                movefile(outfilename,fullfile(baseDir,imagingDirRaw,[subj_name{sn(s)} '-n']))
                cd ..
                rmdir('temp','s');
                fprintf('Run %d done for %s \n',i,subj_name{sn(s)});
            end;
        end
    case 'FUNC:realign'               % Realign functional images (both sessions)
        % SPM realigns all volumes to the first volume of first run
        % example: bsp_imana('FUNC:realign',1,[1:4])
        
        sn=varargin{1}; %subjNum
        runs=varargin{2}; % runNum
        
        subjs=length(sn);
        
        for s=1:subjs,
            
            cd(fullfile(baseDir,imagingDirRaw,[subj_name{sn(s)} '-n']));
            spm_jobman % run this step first
            
            data={};
            for i = 1:length(runs),
                for j=1:numTRs-numDummys;
                    data{i}{j,1}=sprintf('run_%2.2d.nii,%d',runs(i),j);
                end;
            end;
            spmj_realign(data);
            fprintf('runs realigned for %s\n',subj_name{sn(s)});
        end
    case 'FUNC:correct deform'        % Correct Magnetic field deformations using antsRegEpi.sh
    case 'FUNC:move_data'             % Move realigned data
        % Moves image data from imaging_data_raw into imaging_data.
        % example: bsp_imana('FUNC:move_data',1,[1:4])
        sn=varargin{1}; % subjNum
        runs=varargin{2}; % runNum
        
        subjs=length(sn);
        
        for s=1:subjs,
            dircheck(fullfile(baseDir,imagingDir,subj_name{sn(s)}))
            for r=1:length(runs);
                % move realigned data for each run
                source = fullfile(baseDir,imagingDirRaw,[subj_name{sn(s)} '-n'],sprintf('rrun_%2.2d.nii',runs(r)));
                dest = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('run_%2.2d.nii',runs(r)));
                copyfile(source,dest);
                
                % move realignment parameter files for each run
                source = fullfile(baseDir,imagingDirRaw,[subj_name{sn(s)} '-n'],sprintf('rp_run_%2.2d.txt',runs(r)));
                dest = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('rp_run_%2.2d.txt',runs(r)));
                copyfile(source,dest);
            end;            
            
            fprintf('realigned epi''s moved for %s \n',subj_name{sn(s)})
        end
    case 'FUNC:coregTSE'                 % Adjust meanepi to anatomical image REQUIRES USER INPUT
        % (1) Manually seed the functional/anatomical registration
        % - Do "coregtool" on the matlab command window
        % - Select anatomical image and meanepi image to overlay
        % - Manually adjust meanepi image and save result as rmeanepi image
        % example: bsp_imana('FUNC:coreg',1)
        sn=varargin{1};% subjNum
        step=varargin{2}; % first 'manual' then 'auto'
        
        cd(fullfile(baseDir,anatomicalDir,subj_name{sn}));
        
        switch step,
            case 'manual'
                coregtool;
                keyboard();
            case 'auto'
                % do nothing
        end
        
        % (2) Automatically co-register functional and anatomical images for study 1
        J.ref = {fullfile(baseDir,anatomicalDir,subj_name{sn},'anatomical.nii')};
        J.source = {fullfile(baseDir,anatomicalDir,subj_name{sn},'tse.nii')};
        J.other = {''};
        J.eoptions.cost_fun = 'nmi';
        J.eoptions.sep = [4 2];
        J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        J.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate=J;
        spm_jobman('run',matlabbatch);
        
        % NOTE:
        % Overwrites meanepi, unless you update in step one, which saves it
        % as rmeanepi.
        % Each time you click "update" in coregtool, it saves current
        % alignment by appending the prefix 'r' to the current file
        % So if you continually update rmeanepi, you'll end up with a file
        % called r...rrrmeanepi.
    case 'FUNC:coregEPI'      % Adjust meanepi to anatomical image REQUIRES USER INPUT
        % (1) Manually seed the functional/anatomical registration
        % - Do "coregtool" on the matlab command window
        % - Select anatomical image and meanepi image to overlay
        % - Manually adjust meanepi image and save result as rmeanepi image
        % example: bsp_imana('FUNC:coregEPI',1)
        sn=varargin{1};% subjNum
        
        % (2) Automatically co-register functional and anatomical images for study 1
        J.ref = {fullfile(baseDir,anatomicalDir,subj_name{sn},'tse.nii')};
        J.source = {fullfile(baseDir,imagingDirRaw,[subj_name{sn} '-n'],'rmean_epi.nii')};
        J.other = {''};
        J.eoptions.cost_fun = 'nmi';
        J.eoptions.sep = [4 2];
        J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        J.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate=J;
        spm_jobman('run',matlabbatch);
        
    case 'FUNC:make_samealign'        % Align functional images to rmeanepi of study 1
        % Aligns all functional images from both sessions (each study done separately)
        % to rmeanepi of study 1
        % example: bsp_imana('FUNC:make_samealign',1,[1:32])
        sn=varargin{1}; % subjNum
        runs=varargin{2}; % runNum
        
        subjs=length(sn);
        
        for s=1:subjs,
            
            cd(fullfile(baseDir,imagingDir,subj_name{sn(s)}));
            
            % Select image for reference
            % For ants-registered data: TSE 
            % P{1} = fullfile(fullfile(baseDir,anatomicalDir,subj_name{sn},'tse.nii'));
            % for tradition way: rmeanepi 
            P{1} = fullfile(fullfile(baseDir,imagingDir,subj_name{sn},'rmean_epi.nii'));
            
            % Select images to be realigned
            Q={};
            for r=1:numel(runs)
              Q{end+1}    = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('run_%2.2d.nii',runs(r)));
            end;
            
            % Run spmj_makesamealign_nifti
            spmj_makesamealign_nifti(char(P),char(Q));
            fprintf('functional images realigned for %s \n',subj_name{sn(s)})
        end
        
        %spmj_checksamealign
    
    case 'SUIT:isolate'               % Segment cerebellum into grey and white matter
        
        sn=varargin{1};
        %         spm fmri
        for s=sn,
            suitSubjDir = fullfile(baseDir,suitDir,'anatomicals',subj_name{s});dircheck(suitSubjDir);
            source=fullfile(baseDir,anatomicalDir,subj_name{s},'anatomical.nii');
            dest=fullfile(suitSubjDir,'anatomical.nii');
            copyfile(source,dest);
            cd(fullfile(suitSubjDir));
            suit_isolate_seg({fullfile(suitSubjDir,'anatomical.nii')},'keeptempfiles',1);
        end
    case 'SUIT:normalise_dentate'     % Uses an ROI from the dentate nucleus to improve the overlap of the DCN
        % Create the dentate mask in the imaging folder using the tse 
        sn=varargin{1}; %subjNum
        % example: 'sc1_sc2_imana('SUIT:normalise_dentate',1)
        
        cd(fullfile(baseDir,suitDir,'anatomicals',subj_name{sn}));
        job.subjND.gray       = {'c_anatomical_seg1.nii'};
        job.subjND.white      = {'c_anatomical_seg2.nii'};
        job.subjND.dentateROI = {'dentate_mask.nii'};
        job.subjND.isolation  = {'c_anatomical_pcereb.nii'};
        
        
        suit_normalize_dentate(job);
    case 'SUIT:reslice'               % Reslice the contrast images from first-level GLM
        % example: bsm_imana('SUIT:reslice',1,4,'betas','cereb_prob_corr_grey')
        % make sure that you reslice into 2mm^3 resolution
        sn=varargin{1}; % subjNum
        glm=varargin{2}; % glmNum
        region=varargin{3}; % 'betas' or 'contrast' or 'ResMS' or 'cerebellarGrey'
        mask=varargin{4}; % 'cereb_prob_corr_grey' or 'cereb_prob_corr' or 'dentate_mask'
        
        subjs=length(sn);
        
        for s=1:subjs,
            switch region
                case 'betas'
                    glmSubjDir = fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm),subj_name{sn(s)});
                    outDir=fullfile(baseDir,suitDir,sprintf('glm%d',glm),subj_name{sn(s)});
                    images='beta_0';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                case 'contrast'
                    glmSubjDir = fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm),subj_name{sn(s)});
                    outDir=fullfile(baseDir,suitDir,sprintf('glm%d',glm),subj_name{sn(s)});
                    images='con';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                case 'ResMS'
                    glmSubjDir = fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm),subj_name{sn(s)});
                    outDir=fullfile(baseDir,suitDir,sprintf('glm%d',glm),subj_name{sn(s)});
                    images='ResMS';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                case 'cerebellarGrey'
                    source=dir(fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'c1anatomical.nii')); % image to be resliced
                    cd(fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)}));
            end
            job.subj.affineTr = {fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'Affine_c_anatomical_seg1.mat')};
            job.subj.flowfield= {fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'u_a_c_anatomical_seg1.nii')};
            job.subj.resample = {source.name};
            job.subj.mask     = {fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},sprintf('%s.nii',mask))};
            job.vox           = [1 1 1];
            suit_reslice_dartel(job);
            if ~strcmp(region,'cerebellarGrey'),
                source=fullfile(glmSubjDir,'*wd*');
                dircheck(fullfile(outDir));
                destination=fullfile(baseDir,suitDir,sprintf('glm%d',glm),subj_name{sn(s)});
                movefile(source,destination);
            end
            fprintf('%s have been resliced into suit space \n',region)
        end
        
    case 'GLM:makeMask' 
        sn=varargin{1}; % subjNum
        tissues = [1:2 7 8]; 
            
        P{1} = fullfile(fullfile(baseDir,imagingDir,subj_name{sn},'rmean_epi.nii'));
        for i=1:length(tissues) 
            P{i+1}=fullfile(baseDir,suitDir,'anatomicals',subj_name{sn},sprintf('c%danatomical.nii',tissues(i))); 
        end
        out =  fullfile(fullfile(baseDir,imagingDir,subj_name{sn},'brain_mask.nii'));
        spm_imcalc_ui(char(P),out,'i1>800 & (i2+i3+i4+i5)>0.7');
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
                    N{i} = [fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('run-%2.2d.nii,%d',runs(r),i))]; 
                end;
                J.sess(r).scans= N; % number of scans in run
                
                % loop through tasks
                itt = 1;
                for it = 1:nTask 
                    % The order of tasks are different for each run, to
                    % have a common order for the tasks, I will be reading
                    % from the Cc file for all the runs and subjects
                    ST = find(strcmp(P.taskName,Tasks{it}));
                    for taskType = 1:2 % there are two taskTypes: instruction; not instructions
                        if taskType == 1 % instructions
                            % get the isntruction onset
                            instruct_onset = P.realStartTime(ST)- J.timing.RT*numDummys; %% get the instruction start time for the first task 
                            
                            % filling in the fields for SPM_info.mat
                            S.task      = 0;              % task Number
                            S.TN        = {'Instruct'};   % task name (TN)
                            S.inst      = 1;              % is it instruction (1) or not (0)?
                            S.instOrder = ST;             % instOrder is defined by the task that comes after the instruction
                            S.time      = instruct_onset; % instruction onset time
                            % Determine taskName_after and taskName_before
                            % this instruction
                            S.taskName_after  = P.taskName(ST);
                            if ST > 1
                                S.taskName_before = P.taskName(ST-1);
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
                            onset = P.realStartTime(ST) - J.timing.RT*numDummys +announceTime;
                            
                            % filling in the fields for SPM_info.mat
                            S.task      = it;
                            S.TN        = {Tasks{it}};
                            S.inst      = 0;
                            S.instOrder = 0;
                            S.time      = onset;
                            S.taskName_after  = {'none'}; % taskName before and after are only defined for instructions
                            S.taskName_before = {'none'};
                            
                            T  = addstruct(T, S);
                            
                            % filling in the fields for SPM.mat
                            J.sess(r).cond(itt).name     = Tasks{it};
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
            J.cvi_mask = {fullfile(baseDir,imagingDir,subj_name{sn(s)},'brain_mask.nii')};
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
                        tmp(:, strcmp(T.TN, 'rest')) = 1;
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
                
                name = sprintf('%s-%s',char(unique(T.TN(T.(type) == ucondition(tt)))), con_vs);
                
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
        % Example1: bsp_imana('GLM:contrast', 'sn', 3, 'glm', 1, 'type', 'task')
        % Example2: bsp_imana('GLM:contrast', 'sn', 3, 'glm', 1, 'type', 'cond')
        
        sn             = 3;             %% list of subjects
        glm            = 1;             %% The glm number
        type           = 'task';        %% it can be set to task, .... 
        
        vararginoptions(varargin, {'sn', 'glm', 'type'})
                
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, sprintf('GLM_firstlevel_%d', glm));
        
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
                        con(i,T.task==i)=1-1/numTasks;
                        con(i,T.task>0 & T.task~=i)=-1/numTasks;
                    end
            end
            
            SPM.xCon(1) = spm_FcUtil('Set',name, 'F', 'c',con',SPM.xX.xKXs);
            SPM = spm_contrasts(SPM,1:length(SPM.xCon));
            save('SPM.mat', 'SPM','-v7.3');
            SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
            save(fullfile(glmDir, subj_name{s}, 'SPM_light.mat'), 'SPM');

        end % sn      
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
    case 'PHYS:createRegresor'        % Create Retroicor regressors using TAPAS (18 components)
        sn=varargin{1};
        
        for nrun=1:length(run)
            
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
            physio.scan_timing.sqpar.Ndummies = 0;              % Real dummy are not recorded
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
            
            % 7. RETROICOR phase expansion (Glover et al., 2000)
            physio.model.retroicor.include = true;  %% use RETROICOR to calculate regressors
            physio.model.retroicor.order.c = 3;
            physio.model.retroicor.order.r = 4;
            physio.model.retroicor.order.cr = 1;
            
            % 8. RVT respiratory volume per time (Birn et al., 2008)
            physio.model.rvt.include = false;
            physio.model.rvt.delays = 0;
            
            % 9. HRV heart rate variability response (Chang et al., 2009)
            physio.model.hrv.include = false;
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
            physio.verbose.fig_output_file = 'PhysIO.fig';  %% output figure filename (output format can be any matlab supported graphics format)
            physio.verbose.use_tabs = false;
            
            physio.ons_secs.c_scaling = 1;
            physio.ons_secs.r_scaling = 1;
            
            % 14. Run main script
            tapas_physio_main_create_regressors(physio);
        end
        
    case 'ROI:define'                 % Defines ROIs for brain structures
        % Before runing this, create masks for different structures 
        sn=varargin{1}; % subjNum
        region=varargin{2}; % Name of ROI, see cases below
        
        subjs=length(sn);
        for s=1:subjs,
            
            switch region
                case 'cerebellum'
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'rc_anatomical_pcereb.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['cerebellum'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                case 'cortical_grey'
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'rc7anatomical.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['cortical_grey'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                case 'cortical_grey_left'
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'rc7anatomical_left.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['cortical_grey'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                case 'cortical_grey_right'
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'rc7anatomical_right.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['cortical_grey'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                case 'cerebellum_grey'
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'rc_anatomical_seg1.nii');
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
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'rpontine_mask.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['pontine'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                case 'olive' 
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'rolive_mask.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['olive'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                case 'csf' 
                    file = fullfile(baseDir,suitDir,'anatomicals',subj_name{sn(s)},'rcsf_mask.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['csf-bs'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
            end
            dircheck(fullfile(baseDir,regDir,'data',subj_name{sn(s)}));
            save(fullfile(baseDir,regDir,'data',subj_name{sn(s)},sprintf('regions_%s.mat',region)),'R');
            fprintf('ROIs have been defined for %s for %s \n',region,subj_name{sn(s)})
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
        % bsp_imana('TS:getRawTs',1,1,'dentate');
        sn=varargin{1}; % subjNum
        glm=varargin{2}; % glmNum
        region=varargin{3}; % ROI

        glmDir =fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm));
        subjs=length(sn);
        
        for s=1:subjs,
            glmDirSubj=fullfile(glmDir, subj_name{sn(s)});
            load(fullfile(glmDirSubj,'SPM.mat'));
            
            % load data
            load(fullfile(baseDir,regDir,'data',subj_name{sn(s)},sprintf('regions_%s.mat',region)));
            SPM=spmj_move_rawdata(SPM,fullfile(baseDir,imagingDir,subj_name{sn(s)}));

            % Get the raw data files
            V=SPM.xY.VY;
            VresMS = spm_vol(fullfile(glmDirSubj,'ResMS.nii'));
            
            % Get time series data
            Y = region_getdata(V,R{1});  % Data is N x P
            resMS = region_getdata(VresMS,R{1});

            dircheck(fullfile(baseDir,regDir,sprintf('glm%d',glm),subj_name{sn(s)}));
            filename=(fullfile(baseDir,regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('rawts_%s.mat',region)));
            save(filename,'Y','resMS','-v7.3');
            fprintf('Raw ts saved for %s (%s) for %s \n',subj_name{sn(s)},sprintf('glm%d',glm),region);
        end            
    case 'TS:getall'                  % Get all raw TS defined
        sn=varargin{1};
        bsp_imana('TS:getRawTs',sn,1,'cerebellum');
        bsp_imana('TS:getRawTs',sn,1,'cerebellum_grey');
        %bsp_imana('TS:getRawTs',sn,1,'dentate');
        bsp_imana('TS:getRawTs',sn,1,'csf');
        %bsp_imana('TS:getRawTs',sn,1,'pontine');
        bsp_imana('TS:getRawTs',sn,1,'olive');      
        
    case 'test_GLM'                   % Get crossval R2 and R from GLM for different designs
        % example: bsp_imana('test_GLM','inK',{'Hpass','CSF'...},'inX',{'Mov',...},'hpf',128,'ridge',0.1);
        sn = 3;
        glm = 1;
        roi = {'dentate'};%{'cerebellum_grey','dentate','brainstem','pontine'};
        hpf = inf; inK = []; inX = []; ridge = [];
        glmDir = fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm));
        
        % Get the posible options to test
        vararginoptions(varargin,{'hpf','inX','inK','ridge'});
        
        % Load SPM file
        glmDirSubj = fullfile(glmDir, subj_name{sn});
        load(fullfile(glmDirSubj,'SPM.mat'));
        nRuns = length(SPM.nscan);
        
        % Create filter
        for r = 1:nRuns
            k.RT = SPM.xY.RT;
            k.row = (SPM.Sess(r).row);
            k.HParam = inf;
            K(r) = spm_filter(k);
        end
        
        % Add regressors of no interest to K.X0
        if ~isempty(inK)
            for t=1:length(inK)
                switch inK{t}
                    case 'Hpass'        % High pass filter
                        k.HParam = hpf;
                        for rn = 1:nRuns
                            HPF = spm_filter(k);
                            K(rn).X0 = [K(rn).X0 HPF.X0];
                        end
                    case 'Mov'          % Realignment parameters
                        % load the motion files of every run
                        for rn = 1:nRuns
                            mov = load(fullfile(baseDir,imagingDir,subj_name{sn},sprintf('rp_run_%02d.txt',rn)));
                            mov = bsxfun(@rdivide,mov,sqrt(sum(mov.^2)));
                            K(rn).X0 = [K(rn).X0 mov];
                        end
                    case 'MovPCA'       % 2 PCs of realigment paramenters
                        % load the motion files of every run
                        for rn = 1:nRuns
                            mov = load(fullfile(baseDir,imagingDir,subj_name{sn},sprintf('rp_run_%02d.txt',rn)));
                            % Get the principal components
                            [~,score] = pca(mov);
                            score = bsxfun(@rdivide,score,sqrt(sum(score.^2))); % Make regressors orthonormal
                            K(rn).X0 = [K(rn).X0 score(:,1:2)];
                        end
                    case 'Retroicor'    % Retroicor of cardio and resp 18 comp
                        % Load the physio regressors from TAPAS
                        for rn = 1:nRuns
                            phys = load(fullfile(baseDir,'physio',subj_name{sn},sprintf('run%02d',rn),'multiple_regressors.txt'));
                            phys = bsxfun(@rdivide,phys,sqrt(sum(phys.^2))); % Make regressors orthonormal
                            K(rn).X0 = [K(rn).X0 phys];
                        end
                    case 'CSF'          % Mean signal of the CSF around brainstem
                        % Get the CSF data
                        csf = load(fullfile(baseDir,regDir,sprintf('glm%d',glm),subj_name{sn},'rawts_csf.mat'));
                        % mean CSF signal
                        for rn = 1:nRuns
                            mcsf = mean(csf.Y(SPM.Sess(rn).row,:),2);
                            mcsf = bsxfun(@rdivide,mcsf,sqrt(sum(mcsf.^2)));
                            K(rn).X0 = [K(rn).X0 mcsf];
                        end
                    case 'CSFPCA'       % 2 Pcs of CSF
                        % Get the CSF data
                        csf = load(fullfile(baseDir,regDir,sprintf('glm%d',glm),subj_name{sn},'rawts_csf.mat'));
                        % Compute the PCs
                        for rn = 1:nRuns
                            runcsf = csf.Y(SPM.Sess(rn).row,:);
                            % Get the principal components
                            [~,score] = pca(runcsf);
                            score = bsxfun(@rdivide,score,sqrt(sum(score.^2))); % Make regressors orthonormal
                            K(rn).X0 = [K(rn).X0 score(:,1:2)];
                        end
                end
            end
        end
        
        % Get weight/whitening matrix
        W = SPM.xX.W;
        % Get the design matrix (only task related and intercepts)
        X = SPM.xX.X;
        
        % Add regressor of no interest to X
        if ~isempty(inX)
            for t=1:length(inX)
                switch inX{t}
                    case 'Hpass'        % High pass filter
                        k.HParam = hpf;
                        col = 1;
                        for rn = 1:nRuns
                            HPF = spm_filter(k);
                            F(SPM.Sess(rn).row,col:col+size(HPF.X0,2)-1) = HPF.X0;
                            col = col+size(HPF.X0,2);
                        end
                        X = [X F];
                    case 'Mov'          % realignment parameters
                        % load the motion files of every run
                        M = zeros(size(X,1),nRuns*6);
                        col = 1;
                        for rn = 1:nRuns
                            mov = load(fullfile(baseDir,imagingDir,subj_name{sn},sprintf('rp_run_%02d.txt',rn)));
                            mov = bsxfun(@rdivide,mov,sqrt(sum(mov.^2)));
                            M(SPM.Sess(rn).row,col:col+5) = mov;
                            col = col + 6;
                        end
                        X = [X M];
                    case 'MovPCA'       % 2 PCs of realignment parameters
                        % load the motion files of every run
                        M = zeros(size(X,1),nRuns*2);
                        for rn = 1:nRuns
                            mov = load(fullfile(baseDir,imagingDir,subj_name{sn},sprintf('rp_run_%02d.txt',rn)));
                            [~,score] = pca(mov);
                            score = bsxfun(@rdivide,score,sqrt(sum(score.^2)));
                            M(SPM.Sess(rn).row,rn) = score(:,1);
                            M(SPM.Sess(rn).row,rn+nRuns) = score(:,2);
                        end
                        X = [X M];
                    case 'Retroicor'    % Retroicor of cardio and resp 18 comp
                        % Load the physio regressors from TAPAS
                        P = zeros(size(X,1),nRuns*18);
                        col = 1;
                        for rn = 1:nRuns
                            phys = load(fullfile(baseDir,'physio',subj_name{sn},sprintf('run%02d',rn),'multiple_regressors.txt'));
                            phys = bsxfun(@rdivide,phys,sqrt(sum(phys.^2)));
                            P(SPM.Sess(rn).row,col:col+17) = phys;
                            col = col + 18;
                        end
                        X = [X P];
                    case 'CSF'          % Mean signal of the CSF around brainstem
                        % Get the CSF mask
                        csf = load(fullfile(baseDir,regDir,sprintf('glm%d',glm),subj_name{sn},'rawts_csf.mat'));
                        % Include mean CSF in design
                        mcsf = mean(csf.Y,2);
                        CSF = zeros(size(X,1),nRuns);
                        for rn = 1:nRuns
                            CSF(SPM.Sess(rn).row,rn) = mcsf(SPM.Sess(rn).row);
                        end
                        CSF = bsxfun(@rdivide,CSF,std(CSF))*0.1; % z-standardize the regressors
                        X = [X CSF];
                    case 'csfPCAover'   % 2 PCs of CSF computed over 4 runs
                        % Get the CSF mask
                        csf = load(fullfile(baseDir,regDir,sprintf('glm%d',glm),subj_name{sn},'rawts_csf.mat'));
                        [~,score] = pca(csf.Y);
                        % Include the first 2 PCs of CSF in design
                        CSF = zeros(size(X,1),nRuns*2);
                        for rn = 1:nRuns
                            CSF(SPM.Sess(rn).row,rn) = score(SPM.Sess(rn).row,1);
                            CSF(SPM.Sess(rn).row,rn+nRuns) = score(SPM.Sess(rn).row,2);
                        end
                        CSF = bsxfun(@rdivide,CSF,std(CSF))*0.1; % z-standardize the regressors
                        X = [X CSF];
                    case 'CSFPCA'       % 2 PCs of CSF computed independently
                        % Get the CSF mask
                        csf = load(fullfile(baseDir,regDir,sprintf('glm%d',glm),subj_name{sn},'rawts_csf.mat'));
                        % Include the first 2 PCs of CSF in design
                        CSF = zeros(size(X,1),nRuns*2);
                        for rn = 1:nRuns
                            runcsf = csf.Y(SPM.Sess(rn).row,:);
                            [~,score] = pca(runcsf);
                            score = bsxfun(@rdivide,score,sqrt(sum(score.^2)));
                            CSF(SPM.Sess(rn).row,rn) = score(:,1);
                            CSF(SPM.Sess(rn).row,rn+nRuns) = score(:,2);
                        end
                        CSF = bsxfun(@rdivide,CSF,std(CSF))*0.1; % z-standardize the regressors
                        X = [X CSF];
                end
            end
        end
        
        % Design space and projector matrix
        xKXs        = spm_sp('Set',spm_filter(K,W*X));    % KWX
        xKXs.X      = full(xKXs.X);
        pKX         = spm_sp('x-',xKXs);                  % Projector
        
        stats = zeros(length(roi),3);
        % Get the Data (Y)
        for r=1:length(roi) % Loop over ROIs
            
            % Load raw time series
            load(fullfile(baseDir,regDir,sprintf('glm%d',glm),subj_name{sn},sprintf('rawts_%s.mat',roi{r})));     
            
            % Whiten/Weight data and remove filter confounds
            KWY = spm_filter(K,W*Y);
            
            % Weighted Least Squares estimation / Ridge regression 
            if ~isempty(ridge)
                % Tikhonov Matrix
                T = eye(size(xKXs.X,2)) * ridge;
                
                % Number of T conditions on X
                Ncondx = 0;
                for i = 1:nRuns
                    Ncondx = Ncondx + length(SPM.Sess(i).U);
                end
                % Set to zero the ridge for the intercept
                T(Ncondx+1:Ncondx+nRuns,Ncondx+1:Ncondx+nRuns) = 0;
                
                % Projector matrix
                pKXr = inv(xKXs.X' * xKXs.X + T) * xKXs.X';
                B = pKXr * KWY;
            else
                B = pKX * KWY;                        % Original Parameter estimates
            end
            
            % Compute Residuals
            Bres = spm_sp('r',xKXs,KWY);              
            resMS = sum(Bres.^2)./SPM.xX.trRV*mean(diag(SPM.xX.Bcov)); % residuals MeanSquared
            
            % Univariate prewithen Betas
            B = bsxfun(@rdivide,B,sqrt(resMS));
            
            % Calculate R2
            part = kron((1:21),ones(1,18))';  % Check for each design!!!
            cond = repmat((1:18)',21,1);      % Hardcoding here!!!
            
            [R2cv,Rcv] = bsp_patternConsistency_crossval(B,part,cond);
            [R2] = rsa_patternConsistency(B,part,cond);
            
            fprintf('R2 = %.3f  R2cv = %.3f  Rcv = %.3f in %s\n',R2,R2cv,Rcv,roi{r});
            stats(r,:)=[R2 R2cv Rcv];
        end
        varargout = {stats};
    case 'graph_testGLM'
        hps = [];
        %{[
        s = bsp_imana('test_GLM');
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',0.1);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',0.2);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',0.3);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',0.4);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',0.5);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',0.6);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',0.7);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',0.8);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',0.9);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',1);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',2);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',3);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',4);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',5);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',6);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',7);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',8);
        hps = [hps;s];
        s = bsp_imana('test_GLM','ridge',9);
        hps = [hps;s];
        %{
        s = bsp_imana('test_GLM');
        hps = [hps;s];
        s = bsp_imana('test_GLM','inK',{'Hpass'},'hpf',64);
        hps = [hps;s];
        s = bsp_imana('test_GLM','inX',{'Hpass'},'hpf',64);
        hps = [hps;s];
        s = bsp_imana('test_GLM','inK',{'Hpass'},'hpf',128);
        hps = [hps;s];
        s = bsp_imana('test_GLM','inX',{'Hpass'},'hpf',128);
        hps = [hps;s];
        s = bsp_imana('test_GLM','inK',{'Hpass'},'hpf',256);
        hps = [hps;s];
        s = bsp_imana('test_GLM','inX',{'Hpass'},'hpf',256);
        hps = [hps;s];
        s = bsp_imana('test_GLM','inK',{'Mov'});
        hps = [hps;s];
        s = bsp_imana('test_GLM','inX',{'Mov'});
        hps = [hps;s];
        s = bsp_imana('test_GLM','inK',{'MovPCA'});
        hps = [hps;s];
        s = bsp_imana('test_GLM','inK',{'MovPCA'});
        hps = [hps;s];
        s = bsp_imana('test_GLM','inK',{'Retroicor'});
        hps = [hps;s];
        s = bsp_imana('test_GLM','inX',{'Retroicor'});
        hps = [hps;s];
        s = bsp_imana('test_GLM','inK',{'CSF'});
        hps = [hps;s];
        s = bsp_imana('test_GLM','inX',{'CSF'});
        hps = [hps;s];
        s = bsp_imana('test_GLM','inK',{'CSFPCA'});
        hps = [hps;s];
        s = bsp_imana('test_GLM','inX',{'CSFPCA'});
        hps = [hps;s];
        s = bsp_imana('test_GLM','inK',{'allPCA'});
        hps = [hps;s];
        s = bsp_imana('test_GLM','inX',{'allPCA'});
        hps = [hps;s];
  %}
        subplot(3,1,1);bar(hps(:,1));set(gca,'ylim',[0.1 0.27],'xlim',[0 20],'xtick',[]);ylabel('R2');hold on
        line([1,18],[hps(1,1),hps(1,1)],'color','r');hold off;
        subplot(3,1,2);bar(hps(:,2));set(gca,'ylim',[-0.5 -0.1],'xlim',[0 20],'xtick',[]);ylabel('R2cv');hold on
        line([1,18],[hps(1,2),hps(1,2)],'color','r');hold off;
        subplot(3,1,3);bar(hps(:,3));set(gca,'ylim',[-0.5 0.1],'xlim',[0 20],'xtick',(1:20));ylabel('Rcv');hold on
        line([1,18],[hps(1,3),hps(1,3)],'color','r');hold off;
        set(gca,'xticklabel',{'none','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8',...
            '0.9','1.0','2.0','3.0','4.0','5.0','6.0','7.0','8.0','9.0'},'XTickLabelRotation',45)
        %set(gca,'xticklabel',{'none','inK-hpf64','inX-hpf64','inK-hpf128','inX-hpf128','inK-hpf256','inX-hpf256','inK-Mov','inX-Mov',...
         %   'inK-MovPCA','inX-MovPCA','inK-Retro','inX-Retro','inK-CSF','inX-CSF','inK-CSFPCA','inX-CSFPCA','inK-allPCA','inX-allPCA',...
          %  },'XTickLabelRotation',45)
        
        varargout = {hps};
        
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
end

% Local functions
function dircheck(dir)
if ~exist(dir,'dir');
    warning('%s doesn''t exist. Creating one now. You''re welcome! \n',dir);
    mkdir(dir);
end

%%%  Pilot 3 log
%01. Download data from correct7T and check the run numbers and sequences
%02. Copy Mprage and tse in the anatomicals folder
%03. Copy the functional runs in the imaging_data_raw
%04. Remove dummies bsp_imana('FUNC:remDum',3);
%05. Mov Correction bsp_imana('FUNC:realign',3,[1:21]);
%06. Correct gradient deformation antsRegEpi.sh 3;
%07. Rename and move func bsp_imana('FUNC:move_data',3,[1:21]);
%08. Coregister tse and Mprage bsp_imana('FUNC:coreg',3,'auto');
%09. Apply the new aligment bsp_imana('FUNC:make_samealign',3,[1:21]);
%10. Make cerebellar mask usisng suit bsp_imana('SUIT:isolate',3);
%11. Manually create the dentate mask from tse and place it in suitDir
%12. Normalise cerebellum using dartel bsp_imana('SUIT:normalise_dentate',3)
%13. Create SPM.mat for GLM1 bsp_imana('GLM:glm1',3,[1:21]);
%14. Run the GLM bsp_imana('GLM:estimate_glm',3,1);
%15. Download from DICOM server the *_PhysioLog and cpoy to Physio folder
%16. Extract log files from dicom bsp_imana('PHYS:extract',3);
%17. Create regressors using TAPAS bsp_imana('PHYS:createRegresor',3);
%18. Manually create mask for pons, olive, CSF(thr .99,eroded) in SUIT/
%19. Create ROI definitions bsp_imana('ROI:defall',3);
%20. Extract raw Timeseries from ROIs bsp_imana('TS:getall',3);
%21. Test GLM1 with different designs
