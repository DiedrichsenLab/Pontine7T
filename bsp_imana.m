function varargout=bsp_imana(what,varargin)
% Function for minimal preprocessing of the Pontine7T data
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
bidsDir         ='BIDS';

loc_AC = [
    0, 0, 0;   
    0, 0, 0;   % IH: put this here to run centering code 
    0, 0, 0   
];


% Load Participant information (make sure you have Dataframe/util in your
% path
pinfo = dload(fullfile(baseDir,'participants1.tsv')); 
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

    case 'ANAT:reslice_LPI'           % Reslice anatomical image within LPI coordinate systems  
        % example: bsp_imana('ANAT:reslice_LPI',1)
        sn  = varargin{1}; % subjNum 
        for s=sn
            % (1) Reslice anatomical image to set it within LPI co-ordinate frames
            source  = fullfile(baseDir,anatomicalDir,subj_name{s},'anatomical_raw.nii'); %originally: baseDir, anatomicalDir, subj_name{s}, 'anatomial.nii'
            dest    = fullfile(baseDir,anatomicalDir,subj_name{s},'anatomical.nii');
            spmj_reslice_LPI(source,'name', dest);
            
            % (2) In the resliced image, set translation to zero
            V               = spm_vol(dest);
            dat             = spm_read_vols(V);
            V.mat(1:3,4)    = [0 0 0];
            spm_write_vol(V,dat);
            disp 'Manually retrieve the location of the anterior commissure (x,y,z) from anatomical.nii before continuing'
        end

    case 'ANAT:center_AC'             % Re-centre AC
        % Description:
        % Recenters the anatomical data to the Anterior Commissure
        % coordiantes. Doing that, the [0,0,0] coordiante of subject's
        % anatomical image will be the Anterior Commissure.

        % You should manually find the voxel coordinates of AC 
        % for each from their anatomical scans and add it to the
        % participants.tsv file under the loc_ACx loc_ACy loc_ACz columns.

        % This function runs for all subjects and sessions.

        % handling input args:
        sn = varargin{1};
        for s=sn
        % path to the raw anatomical:
            anat_file = fullfile(baseDir,anatomicalDir,subj_name{s},[subj_name{sn} '_T1w.nii,1']);
        
            % Get header info for the image:
            V = spm_vol(anat_file);
            % Read the volume:
            dat = spm_read_vols(V);
        
            % changing the transform matrix translations to put AC near [0,0,0]
            % coordinates:
            R = V.mat(1:3,1:3);
            t = -1 * R * [pinfo.locACx(s),pinfo.locACy(s),pinfo.locACz(s)]';
            V.mat(1:3,4) = t;

            % writing the image with the changed header:
            spm_write_vol(V,dat);
        end
        
    case 'ANAT:segmentation'          % Segmentation + Normalisation %IH: Need to check orientations of images so that they match; for now, did this through the SPM window..
        % example: bsp_imana('ANAT:segmentation',1)
        sn=varargin{1}; % subjNum
        SPMhome=fileparts(which('spm.m'));
        J=[];
        for s=sn

            %T1w scan 
            J.channel(1).vols = {fullfile(baseDir,anatomicalDir,subj_name{s},[subj_name{sn} '_T1w.nii,1'])};
            J.channel(1).biasreg = 0.001;
            J.channel(1).biasfwhm = 60;
            J.channel(1).write = [0 1];

            %T2w scan

          %  J.channel(2).vols = {fullfile(baseDir,anatomicalDir,subj_name{s},[subj_name{sn} '_whole_T2w.nii,1'])};
          %  J.channel(2).biasreg = 0.001;
          %  J.channel(2).biasfwhm = 60;
          %  J.channel(2).write = [0 1];
            
            %tissues

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
            fprintf('Check segmentation results for %s\n', subj_name{s})
        end;

    case 'ANAT:T2T1coreg'                                                      
        % coregister the whole brain T2 to the T1 image 
        
        % handling input args:
        sn=varargin{1}; % subjNum
        image = '_TSE_T2w.nii'; 
        % loop on sessions:
        J.source = {fullfile(baseDir,anatomicalDir,subj_name{sn},[subj_name{sn} image])};
        J.ref = {fullfile(baseDir,anatomicalDir,subj_name{sn},[subj_name{sn} '_T1w.nii'])}; 
        J.other = {''};
        J.eoptions.cost_fun = 'nmi';
        J.eoptions.sep = [4 2];
        J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        J.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate=J; 
        spm_jobman('run',matlabbatch);

    case 'FUNC:realign'      

        % realign functional images
        % SPM realigns all volumes to the mean volume of first run

       sn=varargin{1}; % subjNum
       for s = sn
            spm_jobman('initcfg')
            
            data = {};
                % initialize data cell array which will contain file names for runs/TR images
                %func_ses_subj_dir = fullfile(imagingDirRaw ,subj_name{s});

                                
                for r_cell = run(1:min(numel(run)))
                    current_run = str2double(r_cell{1});
                    % Obtain the number of TRs for the current run
                    for j = 1:numTRs
                        data{current_run}{j,1} = fullfile(imagingDirRaw,strcat(subj_name{s},'-n'),sprintf('uint16_run-%02d.nii,%d', current_run, j));
                    end % j (TRs/images)
                end % r (runs)            
            spmj_realign(data);
            fprintf('- runs realigned for %s  ',subj_name{s});

        end % s (sn)

    case 'FUNC:move_data'             % Move realigned data
        % Moves image data from imaging_data_raw into imaging_data.
        % example: bsp_imana('FUNC:move_data',1,[1:16])
        sn=varargin{1}; % subjNum
        runs=varargin{2}; % runNum
        for s=sn
            dircheck(fullfile(baseDir,imagingDir,subj_name{s}))
            for r=1:length(runs);
                % move realigned data for each run
                source = fullfile(baseDir,imagingDirRaw,[subj_name{s} '-n'],sprintf('ruint16_run-%02d.nii',runs(r)));
                dest = fullfile(baseDir,imagingDir,subj_name{s},sprintf('run_%2.2d.nii',runs(r)));
                copyfile(source,dest);
                
                % move realignment parameter files for each run
                source = fullfile(baseDir,imagingDirRaw,[subj_name{s} '-n'],sprintf('rp_uint16_run-%02d.txt',runs(r)));
                dest = fullfile(baseDir,imagingDir,subj_name{s},sprintf('rp_run_%2.2d.txt',runs(r)));
                copyfile(source,dest);
            end;            
            
            fprintf('realigned epi''s moved for %s \n',subj_name{s})
        end
   
    case 'FUNC:make_samealign'      % Align functional images to rmeanepi of run 1, session 1
        % Aligns all functional images from both sessions
        % to rmeanepi of run 1 of session 1
        % example: bsp_imana('FUNC:make_samealign',1,[1:8])
        sn=varargin{1}; % subjNum
        %runnum=varargin{2}; % runNum used for coregistration
        runs=varargin{2}; % runNum to coregister
        
        for s=sn
            
            cd(fullfile(baseDir,imagingDir,subj_name{s}));
            
            % Select image for reference
            P{1} = fullfile(fullfile(baseDir,imagingDir,subj_name{s},sprintf('S10_mean_bold.nii'))); %IH: original was 'rmeanrun_%2.2d.nii', runnum
            
            % Select images to be realigned
            Q={};
            for r=1:numel(runs)
              Q{end+1}    = fullfile(baseDir,imagingDir,subj_name{s},sprintf('run_%2.2d.nii',runs(r)));
            end;

             % Coregistration options
            coreg_options = spm_get_defaults('coreg');
            coreg_options.sep = [4 2]; % Separation between sampling points (mm)
            coreg_options.params = [0 0 0 0 0 0]; % Translation and rotation
            coreg_options.cost_fun = 'nmi'; % Normalized Mutual Information
        
            % Run spmj_makesamealign_nifti with specified coregistration options
            spmj_makesamealign_nifti(char(P), char(Q));
           % spmj_makesamealign_nifti(char(P), char(Q), coreg_options);
            
            fprintf('functional images realigned for %s \n',subj_name{s})
        end

    case 'FUNC:coreg'                                                      
        % coregister rbumean image to anatomical image for each session
        
        % handling input args:
        sn=varargin{1}; % subjNum

        for s=sn
     
        % loop on sessions:
            for r_cell = run(1:min(numel(run),1))
            
                mean_file_name = fullfile(fullfile(baseDir,imagingDir,subj_name{s},[subj_name{sn} '_whole_sbref.nii']));
                J.source = {mean_file_name};
                J.ref = {fullfile(baseDir,anatomicalDir,subj_name{s},[subj_name{sn} '_T1w.nii'])}; 
                J.other = {''};
                J.eoptions.cost_fun = 'nmi';
                J.eoptions.sep = [4 2];
                J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                J.eoptions.fwhm = [7 7];
                matlabbatch{1}.spm.spatial.coreg.estimate=J;
                spm_jobman('run',matlabbatch);
            
            % (3) Check alignment manually by using fsleyes similar to step
            % one.
            end
        end 
   
    
    case 'SUIT:isolate'             % Segment cerebellum into grey and white matter
        % LAUNCH SPM FMRI BEFORE RUNNING!!!!!
        % example: 'bsp_imana('SUIT:isolate',1)'
        sn=varargin{1};
        %         spm fmri
        for s=sn
            suitSubjDir = fullfile(baseDir,suitDir,'anatomicals',subj_name{s});dircheck(suitSubjDir);
            % Copy over the T1map image 
            T1name = fullfile(baseDir,anatomicalDir,subj_name{s},[subj_name{sn} '_T1map.nii']); 
            dest=fullfile(suitSubjDir,[subj_name{sn} '_T1w.nii']);
            copyfile(T1name,dest);
            % If available, reslice the T2w image into the T1 voxel
            % resolution 
            if (pinfo.T2_whole(s) >0)
                source_vol=spm_vol(fullfile(baseDir,anatomicalDir,subj_name{s},[subj_name{sn} '_whole_T2w.nii']));
                target_vol=spm_vol(T1name); 
                outname = fullfile(suitSubjDir,[subj_name{sn} '_reslice_T2w.nii']); 
                spmj_reslice_vol(source_vol,target_vol.dim,target_vol.mat,outname);
                input_files = {T1name,outname}; 
            else
                input_files = {T1name};
            end
            cd(fullfile(suitSubjDir));
            suit_isolate_seg(input_files,'keeptempfiles',1);
        end
        
    case 'SUIT:normalise_dartel'    % Uses an ROI from the dentate nucleus to improve the overlap of the DCN
        % Create the dentate mask in the imaging folder using the tse
        % LAUNCH SPM FMRI BEFORE RUNNING!!!!!
        sn=varargin{1}; %subjNum
        % example: 'bsp_imana('SUIT:normalise_dartel',1)'
        
        cd(fullfile(baseDir,suitDir,'anatomicals',subj_name{sn}));
        job.subjND.gray       = {'c_S12_T1w_seg1.nii'};
        job.subjND.white      = {'c_S12_T1w_seg2.nii'};
        job.subjND.isolation  = {'c_S12_T1w_pcereb_corr.nii'};
        suit_normalize_dartel(job) 

    case 'SUIT:save_dartel_def'    
        % Saves the dartel flow field as a deformation file. 
        for sn = 12 %IH: original was for sn = [1:length(subj_name)]
            cd(fullfile(baseDir,suitDir,'anatomicals',subj_name{sn}));
            anat_name = 'S10_T1w';
            suit_save_darteldef(anat_name);
        end; 
    case 'SUIT:normalise_dentate'   % Uses an ROI from the dentate nucleus to improve the overlap of the DCN
        % Create the dentate mask in the imaging folder using the tse 
        sn=varargin{1}; %subjNum
        % example: 'sc1_sc2_imana('SUIT:normalise_dentate',1)
        
        cd(fullfile(baseDir,suitDir,'anatomicals',subj_name{sn}));
        job.subjND.gray       = {'c_anatomical_seg1.nii'};
        job.subjND.white      = {'c_anatomical_seg2.nii'};
        job.subjND.dentateROI = {'dentate_mask.nii'};
        job.subjND.isolation  = {'c_anatomical_pcereb_corr.nii'};
        suit_normalize_dentate(job);
    
    case 'SUIT:reslice_ana'         % Reslice the contrast images from first-level GLM
        % example: bsm_imana('SUIT:reslice',1,'anatomical')
        % make sure that you reslice into 2mm^3 resolution
        sn=varargin{1}; % subjNum
        image=varargin{2}; % 'betas' or 'contrast' or 'ResMS' or 'cerebellarGrey'
        
        for s=sn
            suitSubjDir = fullfile(baseDir,suitDir,'anatomicals',subj_name{s});             
            job.subj.affineTr = {fullfile(suitSubjDir ,sprintf('Affine_c_%s_T1w_seg1.mat', subj_name{s}))};
            job.subj.flowfield= {fullfile(suitSubjDir ,sprintf('u_a_c_%s_T1w_seg1.nii', subj_name{s}))};
            job.subj.mask     = {fullfile(suitSubjDir ,sprintf('c_%s_T1w_pcereb_corr.nii', subj_name{s}))};
            switch image
                case 'anatomical'
                    sourceDir = suitSubjDir; 
                    source = fullfile(sourceDir,[subj_name{sn} '_T1w.nii']); 
                    job.subj.resample = {source};
                    outDir = suitSubjDir; 
                    %outDir = baseDir;
                    job.vox           = [1 1 1];
                case 'TSE'
                    sourceDir = fullfile(baseDir,'anatomicals',subj_name{s}); 
                    source = fullfile(sourceDir,'rtse.nii');
                    job.subj.resample = {source};
                    outDir = suitSubjDir; 
                    job.vox           = [1 1 1];
                case 'csf'
                    sourceDir = fullfile(baseDir,'anatomicals', subj_name{s}, 'mp2rage - T1w'); 
                    source = fullfile(sourceDir,'c3anatomical.nii');
                    job.subj.resample = {source};
                    job.subj.mask     =  {}; 
                    outDir = suitSubjDir; 
                    job.vox           = [1 1 1];
                case 'functional'
                    taskNames = {"finger_sequence"};
                    for taskIdx = 1:numel(taskNames)
                        taskName = taskNames{taskIdx};
                        sourceDir = fullfile(baseDir,'GLM_firstlevel_2',subj_name{s});
                        source = fullfile(sourceDir, sprintf('con_%s_16_runs.nii', taskName));
                        job.subj.resample = {source};
                        outDir = fullfile(baseDir,suitDir,'glm2',subj_name{s}); 
                    end 
            end
            suit_reslice_dartel(job);   
            source=fullfile(sourceDir,'*wd*');
            movefile(source,outDir);
        end

    case 'SURF:reconall'       % Freesurfer reconall routine
        % Calls recon-all, which performs, all of the
        % FreeSurfer cortical reconstruction process
        % Example usage: bsp_imana('SURF:reconall')
        
        sn   = good_subj; % subject list
        vararginoptions(varargin, {'sn'});
           
        % check if freesurfer directory already exists
        freesurferDir   ='/surfaceFreesurfer';
        dircheck(fullfile(baseDir, freesurferDir));
        
        for s = sn
            % get the subject id folder name
            freesurfer_reconall(fullfile(baseDir, freesurferDir), subj_name{s}, ...
                                    fullfile(baseDir, anatomicalDir, subj_name, 'anatomical.nii'));
        end % s (sn)
    case 'SURF:fs2wb'          % Resampling subject from freesurfer fsaverage to fs_LR
        % Example usage: bsp_imana('SURF:fs2wb', 'sn', 2 , 'res', 32)
        
        sn   = good_subj; % subject list
        res  = 32;          % resolution of the atlas. options are: 32, 164
        hemi = [1, 2];      % list of hemispheres
        
        vararginoptions(varargin, {'sn', 'res', 'hemi'});
        
        % setting some directories:
%         fs_LRDir = fullfile(baseDir, sprintf('%s_%d', atlas, res));
        
        freesurferDir   ='/surfaceFreesurfer';
        
        for s = sn 
            % get the subject id folder name
            outDir   = fullfile(baseDir, 'surfaceWB', 'data'); dircheck(outDir);
            surf_resliceFS2WB(subj_name{s}, fullfile(baseDir, freesurferDir), outDir, 'hemisphere', hemi, 'resolution', sprintf('%dk', res))
        end % s (sn)    
    
    case 'ROI:group_define'         % Defines the group regions from the group-space images
        regions={'cerebellum_gray','dentate','pontine','olive','rednucleus'};
        outfilename = 'regions.mat'; 
        vararginoptions(varargin,{'regions','outfilename'}); 
        regGroupDir = fullfile(baseDir,'RegionOfInterest','regdef','group');
        for r = 1:length(regions)
            file = fullfile(regGroupDir,sprintf('%s_mask.nii',regions{r}));
            R{r}.type = 'roi_image';
            R{r}.file= file;
            R{r}.name = regions{r};
            R{r}.value = 1;
        end
        R=region_calcregions(R);     % Calculate the ROI coordinates             
        save(fullfile(regGroupDir,outfilename),'R');

    case 'ROI:group_define_csf' 
        regCsfDir = fullfile(baseDir, 'RegionOfInterest', 'regdef', 'group', 'csf_separate_masks');
        regions={'csf_','csf_'};
        outfilename = 'regions_csf.mat'; 
        vararginoptions(varargin,{'regions', 'outfilename'});
        for r = 1:length(regions)
            file = fullfile(regCsfDir, sprintf('%s_mask.nii', regions{r}));
            R{r}.type = 'csf_image';
            R{r}.file = file;
            R{r}.name = regions{r};
            R{r}.value = 1;
        end 
        R=region_calcregions(R);
        save(fullfile(regCsfDir, outfilename), 'R')  %IH: edited this to make it run; not sure if the output is correct (contains no data)
        
      %  bsp_imana('ROI:group_define','regions',regions,'outfilename',outfilename); 
   
    case 'ROI:group_exclude'        % OPTIONAL: Checks for overlapping voxels in group space 
        regions={'cerebellum_gray','dentate','pontine','olive','rednucleus'};
        
        vararginoptions(varargin,{'regions'}); 
        regGroupDir = fullfile(baseDir,'RegionOfInterest','regdef','group');
        for r = 1:length(regions)
            V(r)=spm_vol(fullfile(regGroupDir,sprintf('%s_mask.nii',regions{r})));
            X(:,:,:,r)=spm_read_vols(V(r));
            
        end; 
        M=sum(X,4);
        indx=find(M>1);
        [i,j,k]=ind2sub(V(1).dim,indx)
        keyboard; 
    case 'ROI:group_cifti'          % Generate cifti file for ROI labels 
        regGroupDir = fullfile(baseDir,'RegionOfInterest','regdef','group');
        load(fullfile(regGroupDir,'regions.mat'));
        % Make labels 
        for r=1:length(R)
            data{r}=ones(size(R{r}.data,1),1)*r; 
        end; 
        dname = {'roi_label'}; 
        V=spm_vol(fullfile(baseDir,'RegionOfInterest','regdef','group','SUIT.nii')); % space defining image
        C=region_make_cifti(R,V,'data',data,'dtype','scalars');
        cifti_write(C, fullfile(baseDir, 'RegionOfInterest', 'regdef', 'group', 'regions.dscalar.nii'));
    
    case 'ROI:make_mask'            % Generates masks to determine available voxels in individual space 
        sn = 11;
        vararginoptions(varargin,{'sn'}); %example: bsp_imana("ROI:make_mask",'sn',11)
        for s=sn 
            glm_mask = fullfile(baseDir,'GLM_firstlevel_2',subj_name{s},'mask.nii');
            pcorr = fullfile(baseDir,'suit','anatomicals',subj_name{s},['c_' subj_name{sn} '_T1w_pcereb_corr.nii']); 
            outfile = fullfile(baseDir,'RegionOfInterest','regdef',subj_name{s},'pcereb_mask.nii'); 
            Vi(1)=spm_vol(glm_mask); 
            Vi(2)=spm_vol(pcorr); 
            spm_imcalc(Vi,outfile,'i1>0 & i2>0'); 
        end
   
    case 'ROI:deformation'          % Deform ROIs into individual space and retain mapping. 
        sn = 11; 
        saveasimg = 1;
        region_file = 'regions.mat';   % File with group ROI definitions 
        def_dir = 'suit/anatomicals/S08';  % This is where the deformation can be found 
        def_img = 'c_S08_T1w_seg1';  
        vararginoptions(varargin,{'sn','saveasimg','region_file','def_dir','def_img'}); 
        
        % Load the group regions 
        groupDir = fullfile(baseDir,'RegionOfInterest','regdef','group');
        groupR = load(fullfile(groupDir,region_file)); 
        % For all subjects, deform those regions and save as regions 
        for s = sn
            mask = fullfile(baseDir,'RegionOfInterest','regdef',subj_name{s},'pcereb_mask_T1w.nii');
            Vmask = spm_vol(mask); 
            Def = fullfile(baseDir,def_dir,['u_a_' def_img '.nii']);  %            Def = fullfile(baseDir,def_dir,subj_name{s},['u_a_' def_img '.nii']);

            mat = fullfile(baseDir,def_dir,['Affine_' def_img '.mat']);
            R=region_deformation(groupR.R,{Def,mat},'mask', mask);
            outdir = fullfile(baseDir,'RegionOfInterest','regdef',subj_name{s});
            save(fullfile(outdir,[region_file]),'R'); 
            % For testing purposes, we can also save these ROIs as images 
            % in the original ROI space 
            if (saveasimg)
                for r=1:length(R)
                    region_saveasimg(R{r},Vmask,'name',fullfile(outdir,[R{r}.name '_mask.nii'])); 
                end
            end
        end
        
    case 'ROI:getRawTs'                % Get raw timeseries and save them
        % bsp_glm('ROI:getRawTs',1,1);
        sn=varargin{1}; % subjNum
        glm=varargin{2}; % glmNum
        
        glmDir =fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm));
        
        for s=sn,
            glmDirSubj=fullfile(glmDir, subj_name{s});
            load(fullfile(glmDirSubj,'SPM.mat'));
            
            % load data
            R = load(fullfile(baseDir,regDir,'regdef',subj_name{s},sprintf('regions.mat')), 'R');
            
            % SPM=spmj_move_rawdata(SPM,fullfile(baseDir,imagingDir,subj_name{s}));
            
            % Get the raw data files
            V=SPM.xY.VY; %originally: SPM.xY.VY; IH
            VresMS = spm_vol(fullfile(glmDirSubj,'ResMS.nii'));
            
            % Get time series data
            for r = 1:length(R)
                Y = region_getdata(V,R.R{1,r});  % Data is N x P; IH: initially was R(r). What needs to be read is regions{1,1}.R{1,4}.name
                resMS = region_getdata(VresMS,R.R{1,r}); %NOTE. 
                filename=(fullfile(baseDir,regDir,'data',subj_name{s},sprintf('rawts_%s.mat',R.R{1,r}.name)));
                save(filename,'Y','resMS','-v7.3');
                fprintf('Raw ts saved for %s for %s \n',subj_name{s},R.R{1,r}.name);
            end
        end
    case 'ROI:getBetas'       % Get Betas from a specific GLM save as as cifti files 
        % bsp_glm('ROI:getRawTs',1,1);
        sn=varargin{1}; % subjNum
        glm=varargin{2}; % glmNum
        
        glmDir =fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm));
        
        for s=sn,
            fprintf('SN: %d\n',s);
            glmDirSubj=fullfile(glmDir, subj_name{s});
            D = readtable('SPM_info.tsv', 'FileType', 'text', 'Delimiter', '\t');
            %D=load(fullfile(glmDirSubj,'SPM_info.mat'));     %IH: changes SPM_info.tsv to .mat
            % load 
            load(fullfile(baseDir,regDir,'regdef',subj_name{s},'regions.mat'));
            
            
            for i=1:length(D.run)
                Vbeta(i) = spm_vol(fullfile(glmDirSubj,sprintf('beta_%04d.nii',i)));
                vname{i} = sprintf('%03d_run-%02d-cond%s',i,D.run(i),D.taskName{i});
            end;
            Vbeta(i+1) = spm_vol(fullfile(glmDirSubj,'ResMS.nii'));
            
            % Get betas 
            for r = 1:length(R)
                Y = region_getdata(Vbeta,R{r});  % Data is N x P
                Y = Y(1:end-1,:)./sqrt(Y(end,:));        % Prewhiten beta estimates 
                C = region_make_cifti(R(r),Vbeta(1),'data',{Y'},'struct',{'OTHER'},'dtype','scalars','dnames',vname);
                filename=(fullfile(baseDir,regDir,'data',subj_name{s},sprintf('%s_beta_glm%d_%s.dscalar.nii',subj_name{s},glm,R{r}.name)));
                cifti_write(C,filename); 
            end
        end
        
    case 'ROI:data_subj_to_group'
        sn=varargin{1}; % subjNums 
        dname=varargin{2}; % Name of data file ('beta_glm1', 'beta_glm2')
        reg_name = 'regions.mat';  % Name of regions
        vol = spm_vol(fullfile(baseDir,regDir,'regdef','group','dentate_mask.nii')); % Space defining image
        for s=sn,
            fprintf('SN: %d\n',s);
            load(fullfile(baseDir,regDir,'regdef',subj_name{s},reg_name));
            % Load and transform the CIFTI files 
            for r = 1:length(R)
                filename=(fullfile(baseDir,regDir,'data',subj_name{s},sprintf('%s_beta_glm%d_%s.dscalar.nii',subj_name{s},dname,R{r}.name)));
                C = cifti_read(filename);
                Cnew = region_deform_cifti(R{r},C,'vol',vol); 
                filename=(fullfile(baseDir,regDir,'data','group',sprintf('beta_glm%d_%s_%s.dscalar.nii',dname,R{r}.name,subj_name{s})));
                cifti_write(Cnew,filename); 
            end
        end

    case 'ROI:make_gmwm_mask'        % Create GM-WM exclusion mask for defining CSF 
        % example: bsp_imana('ROI:make_gmwm_mask',1)
        sn=varargin{1}; % subjNum
        
        for s=sn
            J.input = {fullfile(baseDir,anatomicalDir,subj_name{s},'c1anatomical.nii')
                       fullfile(baseDir,anatomicalDir,subj_name{s},'c2anatomical.nii')};
            J.output = fullfile(baseDir,anatomicalDir,subj_name{s},'gm_wm_exclusion_mask.nii');
            J.outdir = {fullfile(baseDir,anatomicalDir,subj_name{s})};
            J.expression = '(i1+i2)<1e-6';
            J.var = struct('name', {}, 'value', {});
            J.options.dmtx = 0;
            J.options.mask = 0;
            J.options.interp = 1;
            J.options.dtype = 4;
            matlabbatch{1}.spm.util.imcalc=J;
            spm_jobman('run',matlabbatch);
            
            fprintf('GM-WM exlusion mask created for %s \n',subj_name{s})
        end 
        
    case 'ROI:reslice_gmwm_mask'       %Resample gm-wm exclusion mask into epi resolution
        % usage 'bsp_imana('ROI:reslice_gmwm_mask',1)'
        sn=varargin{1};
        J = [];
        
        
        for s=sn
            J.ref = {fullfile(baseDir,'GLM_firstlevel_1',subj_name{s},'mask.nii')};
            J.source = {fullfile(baseDir,anatomicalDir,subj_name{s},'gm_wm_exclusion_mask.nii')};
            J.roptions.interp = 0;
            J.roptions.wrap = [0 0 0];
            J.roptions.mask = 0;
            J.roptions.prefix = 'r';
        
            matlabbatch{1}.spm.spatial.coreg.write = J;
            spm_jobman('run',matlabbatch);
            fprintf('gm-wm exclusion mask resliced for %s \n',subj_name{s})
        end  
    
        
    case 'ROI:mask_rois'        % Mask CSF rois with subject-specific GM-WM exclusion mask
        % example: bsp_imana('ROI:mask_rois',1)
        sn=varargin{1}; % subjNum
        images = {'csf_galenic','csf_medulla','csf_midbrain','csf_pons','csf_postdrain','csf_transverseL','csf_transverseR','csf_ventricle4'};
        
        
        
        for s=sn
            regCsfDir = fullfile(baseDir,'RegionOfInterest', 'regdef', 'group', 'csf_separate_masks');
            regSubjDir = fullfile(baseDir, 'RegionOfInterest', 'data', subj_name{s});

            for im=1:length(images)
                J.input = {fullfile(regCsfDir,sprintf('%s_mask.nii',images{im}))
                           fullfile(baseDir,anatomicalDir,subj_name{s},'rgm_wm_exclusion_mask.nii')};
                J.output = fullfile(regSubjDir,sprintf('mask_%s.nii',images{im}));
                J.outdir = {fullfile(regSubjDir)};
                J.expression = 'i1.*i2';
                J.var = struct('name', {}, 'value', {});
                J.options.dmtx = 0;
                J.options.mask = 0;
                J.options.interp = 1;
                J.options.dtype = 4;
                matlabbatch{1}.spm.util.imcalc=J;
                spm_jobman('run',matlabbatch);
            end
            fprintf('csf ROIs masked for %s \n',subj_name{s})
        end    
        
    case 'ROI:cerebellar_gray'  %IH: I might remove this (it's mask is worse than that created in ROI: deformation) 
        sn=varargin{1};
        for s = sn 
            suitSubjDir = fullfile(baseDir,suitDir,'anatomicals',subj_name{s});             
            regSubjDir = fullfile(baseDir,'RegionOfInterest','data',subj_name{s});
            P = {fullfile(baseDir,'GLM_firstlevel_1',subj_name{sn},'mask.nii'),...
                fullfile(suitSubjDir,'c_anatomical_seg1.nii'),...
                fullfile(suitSubjDir,'c_anatomical_pcereb_corr.nii')}; 
            outname = fullfile(regSubjDir,'cerebellum_gray_mask.nii'); 
            spm_imcalc(char(P),outname,'i1 & (i2>0.1) & (i3>0.3)',{0,0,1,4}); 
        end
        
        
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
    case 'PHYS:createRegressor'       % Create Retroicor regressors using TAPAS (18 components)

      %  sn=15; 
       % run = [1:8]; 
        %stop = true; 
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
            physio.preproc.cardiac.initial_cpulse_select.max_heart_rate_bpm = 110; % Allow higher heart rate to start
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

        
    otherwise
        error('Unknown option: %s',what);

end
        
 
  
