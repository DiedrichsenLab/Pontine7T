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

    case 'FUNC:create_mean_epis'   % Calculate mean EPIs for runs: NOT NEEDED [?]
        % Calculate mean EPIs for run(s) acquired closest to fieldmaps
        % example bsp_imana('FUNC:create_mean_epis',1,[8 16])
        sn=varargin{1}; % subjNum
        runs=varargin{2}; % runNum
        
        subjs = length(sn);
        
        for s=sn
            for r=1:length(runs);
                in_epi = fullfile(baseDir,imagingDir,subj_name{s},sprintf('run_%2.2d.nii',runs(r)));
                out_meanepi = fullfile(baseDir,imagingDirRaw,[subj_name{s} '-n'],sprintf('meanrun_%2.2d.nii.gz',runs(r)));
                command_meanepi = sprintf('fslmaths %s -Tmean %s', in_epi, out_meanepi)
                system(command_meanepi)
                fprintf('mean epi completed for run %d \n',runs(r))
                
                out = fullfile(baseDir,imagingDirRaw,[subj_name{s} '-n'],sprintf('meanrun_%2.2d.nii',runs(r)));
                command_gunzip = sprintf('gunzip -c %s > %s', out_meanepi, out)
                system(command_gunzip)
                fprintf('gunzip completed for run %d \n',runs(r))
                
                command_rm = sprintf('rm %s',out_meanepi)
                system(command_rm)
                fprintf('gzipped file removed for run %d \n',runs(r))
            end
        end          
        
    case 'FMAP:average_magnitudes'        % Average magnitude images for each session
        % Averages the two magnitude images for each session
        % example: bsp_imana('FMAP:average_magnitudes',1,1)
        sn=varargin{1}; % subjNum
        sessn=varargin{2}; %sessNum
        for s=sn
            cd(fullfile(baseDir,fmapDir,subj_name{s},sprintf('fmap_sess_%d',sessn)));            
            J.input = {sprintf('magnitude1_sess_%d.nii,1',sessn)
                       sprintf('magnitude2_sess_%d.nii,1',sessn)};
            J.output = fullfile(baseDir,fmapDir,subj_name{s},sprintf('fmap_sess_%d',sessn),sprintf('magnitudeavg_sess_%d.nii',sessn));
            J.outdir = {fullfile(baseDir,fmapDir,subj_name{s})};
            J.expression = '(i1+i2)/2';
            J.var = struct('name', {}, 'value', {});
            J.options.dmtx = 0;
            J.options.mask = 0;
            J.options.interp = 1;
            J.options.dtype = 4;
            matlabbatch{1}.spm.util.imcalc=J;
            spm_jobman('run',matlabbatch);
            fprintf('magnitude fieldmaps averaged for %s \n',subj_name{s})
        end    
        
    case 'FUNC:run_feat_coregistration'    %Run run_feat_coregistrations.sh shell script
         % example: bsp_imana('FUNC:run_feat_coregistration',1,1)
         sn=varargin{1}; %subjNum
         sessn=varargin{2}; %sessNum
         
         
         for s=sn
             
            subjID = strip(subj_name{s},'left','S') 
            command_feat = sprintf('bash /srv/diedrichsen/shell/run_feat_coregistration.sh %s %2.2d', subjID, sessn)
            system(command_feat)
            
         end
         
         fprintf('feat coregistration completed for %s \n',subj_name{s})
     
    case 'FUNC:gunzip'        % Unzip .nii.gz file to .nii
        % Run gunzip on the output file from epi_reg step
        % example: bsp_imana('FUNC:gunzip',1,8)
        sn=varargin{1}; % subjNum
        runnum=varargin{2}; %runNum
        
        
        for s=sn
            in     = fullfile(baseDir,imagingDirRaw,[subj_name{s} '-n'],sprintf('meanrun_%2.2d_func2highres.nii.gz',runnum));
            out    = fullfile(baseDir,imagingDirRaw,[subj_name{s} '-n'],sprintf('meanrun_%2.2d_func2highres.nii',runnum));
            % gunzip -c file.gz > /THERE/file
            command = sprintf('gunzip -c %s > %s', in, out)
            system(command)
            fprintf('gunzip completed for %s \n',subj_name{s})
        end
        
    case 'FUNC:coreg_meanepi'       % Coregister meanrun_01 to meanrun_01_func2struct
        % Need meanrun_01 in epi resolution coregistered to anatomical
        % example: bsp_imana('FUNC:coreg_meanepi',1,8)
        sn=varargin{1}; % subjNum
        runnum=varargin{2} %runNum
        
        J = [];
        for s=sn
            
            cd(fullfile(baseDir,imagingDirRaw,[subj_name{s} '-n']));
            
            J.ref = {fullfile(baseDir,imagingDirRaw,[subj_name{s} '-n'],sprintf('meanrun_%2.2d_func2highres.nii',runnum))};
            J.source = {fullfile(baseDir,imagingDirRaw,[subj_name{s} '-n'],sprintf('meanrun_%2.2d.nii',runnum))};
            J.other = {''};
            J.eoptions.cos_fun = 'nmi';
            J.eoptions.sep = [4 2];
            J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            J.eoptions.fwhm = [7 7];
            matlabbatch{1}.spm.spatial.coreg.estimate=J;
            spm_jobman('run',matlabbatch);
            fprintf('mean epi coregistered for %s \n',subj_name{s})
            command = sprintf('cp %s %s',fullfile(baseDir,imagingDirRaw,[subj_name{s} '-n'],sprintf('meanrun_%2.2d.nii',runnum)),fullfile(baseDir,imagingDir,subj_name{s},sprintf('rmeanrun_%2.2d.nii',runnum)))
            system(command)
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
        
 
  
