function varargout=bsp_imana(what,varargin)
% Optional and exploratory steps for processing Pontine7T data
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
pinfo = dload(fullfile(baseDir,'participants.tsv')); 
subj_name = pinfo.participant_id;
good_subj = find(pinfo.good)'; % Indices of all good subjects

%========================================================================================================================
% GLM INFO
numDummys = 3;                                                              % per run
numTRs    = 328;                                                            % per run (includes dummies)
run         = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16'};
runB        = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];  % Behavioural labelling of runs
sess        = [1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2];                  % session number
%========================================================================================================================
switch(what)

    case 'ANAT:bet'                   % Brain extraction for anatomical.nii 
        % Run bash script /srv/diedrichsen/shell/optiBET.sh
        % Edit command variable to set path to optiBET.sh script
        % example: bsp_imana('ANAT:bet',1)
        sn=varargin{1}; % subjNum
        for s=sn
            img    = fullfile(baseDir,anatomicalDir,subj_name{s},'manatomical.nii');
            command = sprintf('bash /srv/diedrichsen/shell/optiBET.sh -i %s', img)
            system(command)
            
            in = fullfile(baseDir,anatomicalDir,subj_name{s},'manatomical_optiBET_brain.nii.gz');
            out = fullfile(baseDir,anatomicalDir,subj_name{s},'manatomical_brain.nii.gz');
            copy_command = sprintf('cp %s %s', in, out)
            system(copy_command)
            
            fprintf('optiBET completed for %s \n',subj_name{s})
            fprintf('Check the output of optiBET using FSLeyes or some other visualization software.')
        end
    
    case 'ANAT:biascorrect_tse'       % Bias correct TSE
        % example: bsp_imana('ANAT:biascorrect_tse',1)
        sn=varargin{1}; %subjNum
        for s=sn
            in_tse = fullfile(baseDir,anatomicalDir,subj_name{s},'tse.nii');
            out_tse = fullfile(baseDir,anatomicalDir,subj_name{s},'tse'); 
            command_bias = sprintf('fsl_anat --nononlinreg --strongbias --nocrop --noreg --nosubcortseg --noseg --clobber -t T2 -i %s -o %s', in_tse, out_tse)
            system(command_bias)
            
            fprintf('tse bias correction completed for %s \n',subj_name{s})
            fprintf('Check the results in FSLeyes or some other visualization software.')
        end        
            
    case 'ANAT:coregister_tse'                % Coregister TSE to anatomical
        % example: bsp_imana('ANAT:coregister_tse',1)
        sn=varargin{1}; % subjNum
        for s=sn
            in_tse = fullfile(baseDir,anatomicalDir,subj_name{s},'tse.anat','T2_biascorr.nii.gz');
            in_ref = fullfile(baseDir,anatomicalDir,subj_name{s},'manatomical.nii');
            out_mat = fullfile(baseDir,anatomicalDir,subj_name{s},'tse_to_anatomical_mi.mat');
            out_tse  = fullfile(baseDir,anatomicalDir,subj_name{s},'tse_to_anatomical_mi');
            command_mask = sprintf('flirt -in %s -ref %s -usesqform -searchrx -45 45 -searchry -45 45 -searchrz -45 45 -dof 6 -cost mutualinfo -omat %s -out %s', in_tse, in_ref, out_mat, out_tse)
            system(command_mask)
            
            fprintf('tse coregistration completed for %s \n',subj_name{s})
            fprintf('Check the results in FSLeyes or some other visualization software.')
        end        

        
    case 'SUIT:reslice'               % Reslice the contrast images from first-level GLM
        % example: bsm_imana('SUIT:reslice',1,4,'betas','cereb_prob_corr_grey')
        % make sure that you reslice into 2mm^3 resolution
        sn=2; % subjNum
        glm=1; % glmNum
        type='contrast'; % 'betas' or 'contrast' or 'ResMS' or 'cerebellarGrey'
        mask='c_anatomical_pcereb_corr'; % 'cereb_prob_corr_grey' or 'cereb_prob_corr' or 'dentate_mask'
        for s=sn
            switch type
                case 'betas'
                    glmSubjDir = fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm),subj_name{s});
                    outDir=fullfile(baseDir,suitDir,sprintf('glm%d',glm),subj_name{s});
                    images='beta_0';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                case 'contrast'
                    glmSubjDir = fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm),subj_name{s});
                    outDir=fullfile(baseDir,suitDir,sprintf('glm%d',glm),subj_name{s});
                    images='con';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                case 'ResMS'
                    glmSubjDir = fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm),subj_name{s});
                    outDir=fullfile(baseDir,suitDir,sprintf('glm%d',glm),subj_name{s});
                    images='ResMS';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                case 'cerebellarGrey'
                    source=dir(fullfile(baseDir,suitDir,'anatomicals',subj_name{s},'c1anatomical.nii')); % image to be resliced
                    cd(fullfile(baseDir,suitDir,'anatomicals',subj_name{s}));
            end
            job.subj.affineTr = {fullfile(baseDir,suitDir,'anatomicals',subj_name{s},'Affine_c_anatomical_seg1.mat')};
            job.subj.flowfield= {fullfile(baseDir,suitDir,'anatomicals',subj_name{s},'u_a_c_anatomical_seg1.nii')};
            job.subj.resample = {source.name};
            job.subj.mask     = {fullfile(baseDir,suitDir,'anatomicals',subj_name{s},sprintf('%s.nii',mask))};
            job.vox           = [1 1 1];
            % Replace Nans with zeros to avoid big holes in the the data 
            for i=1:length(source)
                V=spm_vol(source(i).name); 
                X=spm_read_vols(V); 
                X(isnan(X))=0; 
                spm_write_vol(V,X); 
            end; 
            suit_reslice_dartel(job);
            
            source=fullfile(glmSubjDir,'*wd*');
            dircheck(fullfile(outDir));
            destination=fullfile(baseDir,suitDir,sprintf('glm%d',glm),subj_name{s});
            movefile(source,destination);

            fprintf('%s have been resliced into suit space \n',type)
        end
    case 'SUIT:map_to_flat' 
        % Maps wdcon data to flatmap 
        sn = 2; 
        glm = 1; 
        vararginoptions(varargin,{'sn','glm','type','mask'});
        source_dir=fullfile(baseDir,suitDir,sprintf('glm%d',glm),subj_name{sn});
        source=dir(fullfile(source_dir,'wdcon*')); % images to be resliced
        for i=1:length(source)
            name{i} = fullfile(source_dir,source(i).name); 
        end; 
        MAP = suit_map2surf(name); 
        % G = surf_makeFuncGifti
        set(gcf,'PaperPosition',[2 2 15 7]);
        wysiwyg; 
        for i=1:10 
            subplot(2,5,i); 
            suit_plotflatmap(MAP(:,i),'cscale',[-1.5 1.5]); 
            n = source(i).name(7:end-4); 
            n(n=='_')=' '; 
            title(n); 
        end; 
     
end
        
 
  
