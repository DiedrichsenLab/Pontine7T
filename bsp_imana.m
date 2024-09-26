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
pinfo = dload(fullfile(baseDir,'participants.tsv')); 
subj_name = pinfo.participant_id;
good_subj = find(pinfo.good)'; % Indices of all good subjects

%========================================================================================================================
% GLM INFO
numDummys = 3;                                                              % per run
numTRs    = pinfo.numTR;                                                            % per run (includes dummies)
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
        image = '_whole_T2w.nii'; 
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

    case 'FUNC:split' % for instances where realign doesn't work becuase of 'offending images'; split 4D functional image into individual 3D volumes 
   
        source_image = '/Volumes/Diedrichsen_data$/data/Cerebellum/Pontine7T/imaging_data/S16/run_16.nii'; 
    
        % Load the header information of the 4D image
        hdr = spm_vol(source_image);
    
        % Read the 4D volume data
        data = spm_read_vols(hdr);
    
        % Create a directory to store the split 3D volumes
        output_dir = fullfile(fileparts(source_image), 'split_volumes');
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
    
        % Initialize a cell array to store file paths of split 3D volumes
        split_files = cell(size(data, 4), 1);
    
        % Loop through each volume in the 4D data
        for i = 1:size(data, 4)
            % Extract each 3D volume
            volume = data(:,:,:,i);
        
            % Create a new header for the 3D volume
            hdr_3D = hdr(1);
        
            % Construct file path for the split 3D volume
            volume_filename = sprintf('volume_%d.nii', i);
            output_path = fullfile(output_dir, volume_filename);
        
            % Update the header with the new file path
            hdr_3D.fname = output_path;
        
            % Write the 3D volume to a new NIfTI file
            spm_write_vol(hdr_3D, volume);
        
            % Store the file path of the split 3D volume
            split_files{i} = output_path;
        end
    
        % Display confirmation message
        fprintf('- 4D image split into %d 3D volumes\n', numel(split_files));
    
        % Optionally return the file paths of split 3D volumes
        varargout{1} = split_files;

    case 'FUNC:concatenate' %concatenate 3D images back to 4D after estimate + reslice 

        resliced_dir = fullfile(baseDir,imagingDir,'S11','split_volumes'); % Adjust to your directory

        % Define the output directory for the concatenated 4D image
        output_dir = fullfile(baseDir,imagingDir,'S11'); % Adjust to your directory

        % Get the list of resliced 3D volumes using the wildcard pattern
        volumes = dir(fullfile(resliced_dir, 'rvolume_*.nii'));
        numVolumes = length(volumes);

        % Initialize the header for the 4D NIfTI file using the first volume
        hdr = spm_vol(fullfile(resliced_dir, volumes(1).name));
        data = spm_read_vols(hdr);

        % Initialize the 4D data array
        dim = [hdr.dim, numVolumes];
        data4D = zeros(dim);

    % Read each 3D volume and store it in the 4D array
        for i = 1:numVolumes
            hdr = spm_vol(fullfile(resliced_dir, volumes(i).name));
            data = spm_read_vols(hdr);
            data4D(:,:,:,i) = data;
        end

    % Create a new header for each 3D volume and write to the same NIfTI file
        hdr4D = hdr;  % Use the header from the last volume
        hdr4D.fname = fullfile(output_dir, 'concatenated_4D_image.nii');
        hdr4D.dim = [hdr4D.dim, numVolumes];
      %  hdr4D.n = [1 numVolumes];  % Indicate the number of volumes
        hdr4D.private.dat.fname = hdr4D.fname;
        hdr4D.private.dat.dim = [hdr4D.dim, 1];

        %errors out here 
    % Write each 3D volume to the 4D NIfTI file
        for i = 1:numVolumes
            hdr4D.n(1) = i;  % Set the volume index
            spm_write_vol(hdr4D, data4D(:,:,:,i));
        end

        fprintf('4D image saved as %s\n', hdr4D.fname);



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
            P{1} = fullfile(fullfile(baseDir,imagingDir,subj_name{s},[subj_name{sn},'_mean_bold.nii'])); %IH: original was 'rmeanrun_%2.2d.nii', runnum    
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
            T1name = fullfile(baseDir,anatomicalDir,subj_name{s},[subj_name{sn} '_T1w.nii']); 
            dest=fullfile(suitSubjDir,[subj_name{sn} '_T1w.nii']);
            copyfile(T1name,dest);
            % If available, reslice the T2w image into the T1 voxel
            % resolution 
          
            if (pinfo.T2_whole(s) == 1)
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
        job.subjND.gray       = {['c_', subj_name{sn}, '_T1w_seg1.nii']};
        job.subjND.white      = {['c_', subj_name{sn}, '_T1w_seg2.nii']};
        job.subjND.isolation  = {['c_', subj_name{sn}, '_T1w_pcereb_corr.nii']};
        suit_normalize_dartel(job) 

    case 'SUIT:save_dartel_def'     % Save deformation into SUIT space as deformation ('y_,,') file
        sn = good_subj;
        target_dir = fullfile(baseDir,[regDir '_SUIT'],'deform');
        vararginoptions(varargin,{'sn'}); 

        % Saves the dartel flow field as a deformation file. 
        for s = sn %
            cd(fullfile(baseDir,suitDir,'anatomicals',subj_name{s}));
            anat_name = [subj_name{s},'_T1w'];
            suit_save_darteldef(anat_name);
            copyfile(['y_' anat_name '_suitdef.nii'],fullfile(target_dir,['y_' subj_name{s} '.nii'])); 
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
                    job.vox           = [1 1 1];
                case 'TSE'
                    sourceDir = fullfile(baseDir,'anatomicals',subj_name{s}); 
                    source = fullfile(sourceDir,'rtse.nii');
                    job.subj.resample = {source};
                    outDir = suitSubjDir; 
                    job.vox           = [1 1 1];
                case 'csf'
                    sourceDir = fullfile(baseDir,'anatomicals', subj_name{s}); 
                    source = fullfile(sourceDir,'c3anatomical.nii');
                    job.subj.resample = {source};
                    job.subj.mask     =  {}; 
                    outDir = suitSubjDir; 
                    job.vox           = [1 1 1];
                case 'functional'
                    taskNames = {"flexion_extension", "romance_movie", "finger_sequence", "n_back", "action_observation", "theory_of_mind", "visual_search", "Instruct", "semantic_prediction", "rest"};
                    for taskIdx = 1:numel(taskNames)
                        taskName = taskNames{taskIdx};
                        sourceDir = fullfile(baseDir,'GLM_firstlevel_2',subj_name{s});
                        source = fullfile(sourceDir, sprintf('con_%s_16_runs.nii', taskName));
                        job.subj.resample = {source};
                        outDir = fullfile(baseDir,suitDir,'glm2',subj_name{s}); 
                        suit_reslice_dartel(job);   
                        source=fullfile(sourceDir,'*wd*');
                        movefile(source,outDir);
                    end 
            end
            suit_reslice_dartel(job);   
            source=fullfile(sourceDir,'*wd*');
            movefile(source,outDir);
        end

    case 'SUIT:flatmap'
        %creates flatmaps for each task 

       sn=varargin{1};
       taskNames = {"flexion_extension", "romance_movie", "finger_sequence", "n_back", "action_observation", "theory_of_mind", "visual_search", "Instruct", "semantic_prediction", "rest"};
       for taskIdx = 1:numel(taskNames)
           taskName = taskNames{taskIdx};
           sourceDir = fullfile(baseDir,suitDir, 'glm2', subj_name{sn});
           source = fullfile(sourceDir,sprintf('wdcon_%s.nii',taskName));
           Data = suit_map2surf(source,'space','SUIT');

           %save plots 

           figure; 
           suit_plotflatmap(Data);
           saveas(gcf, fullfile(sourceDir, sprintf('flatmap_%s.png', taskName)));

           %close figure 
           close(gcf);

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
        % EXAMPLE: bsp_imana('ROI:group_define','reg_type','BOLD','regions',{'dentate'}); 
        reg_type = 'SUIT'; % Use 'SUIT' or 'BOLD'
        regions={'cerebellum_gray','dentate','pontine','olive','rednucleus'};
        outfilename = 'regions.mat'; 
        vararginoptions(varargin,{'reg_type','regions','outfilename'}); 
        regGroupDir = fullfile(baseDir,[regDir '_' reg_type],'regdef','group');
        for r = 1:length(regions)
            file = fullfile(regGroupDir,sprintf('%s_mask.nii',regions{r}));
            R{r}.type = 'roi_image';
            R{r}.file= file;
            R{r}.name = regions{r};
            R{r}.value = 1;
        end
        R=region_calcregions(R);     % Calculate the ROI coordinates             
        save(fullfile(regGroupDir,outfilename),'R');
   
    case 'ROI:group_exclude'        
        % OPTIONAL: Checks for overlapping voxels in group space 
        reg_type = 'SUIT'; % Use 'SUIT' or 'BOLD'
        regions={'cerebellum_gray','dentate','pontine','olive','rednucleus'};        
        vararginoptions(varargin,{'reg_type','regions'}); 
        regGroupDir = fullfile(baseDir,[regDir '_' reg_type],'regdef','group');
        for r = 1:length(regions)
            V(r)=spm_vol(fullfile(regGroupDir,sprintf('%s_mask.nii',regions{r})));
            X(:,:,:,r)=spm_read_vols(V(r));
        end; 
        M=sum(X,4);
        indx=find(M>1);
        [i,j,k]=ind2sub(V(1).dim,indx)
        keyboard; 
    case 'ROI:group_cifti'          % Generate cifti file for ROI labels 
        % OPTIONAL: bsp_imana('ROI:group_cifti','reg_type','BOLD'); 
        reg_type = 'SUIT'; % Use 'SUIT' or 'BOLD'
        vararginoptions(varargin,{'reg_type'}); 
        regGroupDir = fullfile(baseDir,[regDir '_' reg_type],'regdef','group');
        load(fullfile(regGroupDir,'regions.mat'));
        % Make labels 
        for r=1:length(R)
            data{r}=ones(size(R{r}.data,1),1)*r; 
        end; 
        dname = {'roi_label'}; 
        V=spm_vol(fullfile(baseDir,[regDir '_' reg_type],'template.nii')); % space defining image
        C=region_make_cifti(R,V,'data',data,'dtype','scalars');
        cifti_write(C, fullfile(regGroupDir, 'regions.dscalar.nii'));
    
    case 'ROI:make_mask'            
       % Generates masks to determine available cerebellar voxels in individual space 
       % This is being saved into the imaging_data directory - as it's
       % indenpendent of the type of ROI 
       % example: bsp_imana("ROI:make_mask",'sn',11)
       sn  = good_subj; 
       vararginoptions(varargin,{'sn'}); 
       for s=sn 
            glm_mask = fullfile(baseDir,'GLM_firstlevel_3',subj_name{s},'mask.nii');
            pcorr = fullfile(baseDir,'suit','anatomicals',subj_name{s},['c_' subj_name{s} '_T1w_pcereb_corr.nii']); 
            outfile = fullfile(baseDir,'imaging_data',subj_name{s},'pcereb_mask.nii'); 
            Vi(1)=spm_vol(glm_mask); 
            Vi(2)=spm_vol(pcorr); 
            % Calculate mask in functional space - datatype to unit8 (2)
            Vo = spm_imcalc(Vi,outfile,'i1>0 & i2>0',{0,0,0,2,'mask'}); 
       end; 
    case 'ROI:deformation'          % Deform ROIs into individual space and retain mapping. 
        % Example bsp_imana('ROI:deformation','reg_type','BOLD','sn',[18 19]);
        sn  = good_subj; 
        saveasimg = 1;
        reg_type = 'SUIT'; % Use 'SUIT' or 'BOLD'
        region_file = 'regions.mat'; 
        vararginoptions(varargin,{'saveasimg','reg_type','sn','region_file'}); 
        groupDir = fullfile(baseDir,[regDir '_' reg_type],'regdef','group');
        defDir = fullfile(baseDir,[regDir '_' reg_type],'deform');
        groupR = load(fullfile(groupDir,region_file)); 
        for s = sn
            mask = fullfile(baseDir,'imaging_data',subj_name{s},'pcereb_mask.nii'); % This is the cerebellar mask in fMRI space
            Vmask = spm_vol(mask); 
            def_file = fullfile(defDir,['y_' subj_name{s} '.nii']); 
            R=region_deformation(groupR.R,def_file,'mask', mask);
            outdir = fullfile(baseDir,[regDir '_' reg_type],'regdef',subj_name{s});
            if ~exist(outdir)
                mkdir(outdir);
            end
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
        % bsp_glm('ROI:getRawTs','sn',[17],'glm',2,'reg_type','BOLD');
        sn=good_subj; % subjNum
        glm=2; % glmNum
        reg_type = 'SUIT'; % Use 'SUIT' or 'BOLD'
        vararginoptions(varargin,{'sn','reg_type','glm'}); 
        
        glm_dir =fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm));
        reg_dir = fullfile(baseDir,[regDir '_' reg_type]); % Specific reg_dir 

        
        for s=sn,
            glmDirSubj=fullfile(glm_dir, subj_name{s});
            outdir = fullfile(reg_dir,'data',subj_name{s});
            if ~exist(outdir)
                mkdir(outdir);
            end

            load(fullfile(glmDirSubj,'SPM.mat'));
            
            % load data
            load(fullfile(reg_dir,'regdef',subj_name{s},sprintf('regions.mat')), 'R');
            
            % SPM=spmj_move_rawdata(SPM,fullfile(baseDir,imagingDir,subj_name{s}));
            
            % Get the raw data files
            V=SPM.xY.VY; 
            
            % Get time series data
            for r = 1:length(R) % 1:length(R)
                Y = region_getdata(V,R{r});  % Data is N x P
                C = region_make_cifti(R(r),V(1),'data',{Y'},'struct',{'OTHER'},'dtype','series');
                filename=(fullfile(outdir,sprintf('%s_%s.dtseries.nii',subj_name{s},R{r}.name)));
                cifti_write(C,filename); 
            end
        end

    case 'ROI:getBetas'     % Get Betas from a specific GLM save as as cifti files 
        % bsp_imana('ROI:getBetas','sn',[18],'glm',2,'reg_type','BOLD');
        sn=good_subj; % subjNum
        glm=2; % glmNum
        reg_type = 'SUIT'; % Use 'SUIT' or 'BOLD'
        vararginoptions(varargin,{'sn','reg_type','glm'}); 
        
        glm_dir =fullfile(baseDir,sprintf('GLM_firstlevel_%d',glm));
        reg_dir = fullfile(baseDir,[regDir '_' reg_type]); % Specific reg_dir 

        
        for s=sn,
            fprintf('SN: %d\n',s);
            glmDirSubj=fullfile(glm_dir, subj_name{s});
            outdir = fullfile(reg_dir,'data',subj_name{s});
            if ~exist(outdir)
                mkdir(outdir);
            end
            
            D = dload(fullfile(glmDirSubj,'SPM_info.tsv'));
            load(fullfile(reg_dir,'regdef',subj_name{s},'regions.mat'));
            
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
                filename=(fullfile(outdir,sprintf('%s_beta_glm%d_%s.dscalar.nii',subj_name{s},glm,R{r}.name)));
                cifti_write(C,filename); 
            end
        end
        
    case 'ROI:data_subj_to_group'
        % Transforms any individual region data file in individual space
        % into group space. Use this mainly on betas 
        % Example: 
        % bsp_imana('ROI:data_subj_to_group','reg_type','BOLD','cifti_name','beta_glm2','sn',[18 19]);
        sn=good_subj; % subjNum
        cifti_name = 'beta_glm2'; % glmNum
        reg_type = 'SUIT'; % Use 'SUIT' or 'BOLD'
        reg_name = 'regions.mat';

        vararginoptions(varargin,{'sn','reg_type','cifti_name'}); 
        reg_dir = fullfile(baseDir,[regDir '_' reg_type]); % Specific reg_dir 
        vol = spm_vol(fullfile(reg_dir,'template.nii')); % Group Space defining image
        for s=sn,
            fprintf('SN: %d\n',s);
            load(fullfile(reg_dir,'regdef',subj_name{s},reg_name));
            % Load and transform the CIFTI files 
            for r = 1:length(R)
                filename=(fullfile(reg_dir,'data',subj_name{s},sprintf('%s_%s_%s.dscalar.nii',subj_name{s},cifti_name,R{r}.name)));
                C = cifti_read(filename);
                Cnew = region_deform_cifti(R{r},C,'vol',vol); 
                filename=(fullfile(reg_dir,'data','group',sprintf('%s_%s_%s.dscalar.nii',cifti_name,R{r}.name,subj_name{s})));
                cifti_write(Cnew,filename); 
            end
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

%manual peak selection 

           % c = [2043	2043	2014	2014	2014	2014	1965	1965	1965	1965	1895	1895	1895	1895	1804	1804	1804	1804	1700	1700	1700	1700	1590	1590	1590	1590	1483	1483	1483	1483	1380	1380	1380	1380	1290	1290	1290	1290	1210	1210	1210	1210	1145	1145	1145	1145	1094	1094	1094	1094	1053	1053	1053	1053	1023	1023	1023	1023	998	998	998	998	979	979	979	979	965	965	965	965	963	963	963	963	974	974	974	974	990	990	990	990	1013	1013	1013	1013	1026	1026	1026	1026	1043	1043	1043	1043	1100	1100	1100	1100	1266	1266	1266	1266	1578	1578	1578	1578	2014	2014	2014	2014	2511	2511	2511	2511	2990	2990	2990	2990	3381	3381	3381	3381	3651	3651	3651	3651	3785	3785	3785	3785	3788	3788	3788	3788	3694	3694	3694	3694	3539	3539	3539	3539	3350	3350	3350	3350	3154	3154	3154	3154	2951	2951	2951	2951	2748	2748	2748	2748	2555	2555	2555	2555	2389	2389	2389	2389	2264	2264	2264	2264	2181	2181	2181	2181	2144	2144	2144	2144	2148	2148	2148	2148	2161	2161	2161	2161	2166	2166	2166	2166	2154	2154	2154	2154	2128	2128	2128	2128	2083	2083	2083	2083	2010	2010	2010	2010	1921	1921	1921	1921	1828	1828	1828	1828	1733	1733	1733	1733	1634	1634	1634	1634	1531	1531	1531	1531	1440	1440	1440	1440	1358	1358	1358	1358	1286	1286	1286	1286	1229	1229	1229	1229	1180	1180	1180	1180	1146	1146	1146	1146	1128	1128	1128	1128	1121	1121	1121	1121	1111	1111	1111	1111	1105	1105	1105	1105	1105	1105	1105	1105	1111	1111	1111	1111	1120	1120	1120	1120	1125	1125	1125	1125	1148	1148	1148	1148	1243	1243	1243	1243	1473	1473	1473	1473	1841	1841];          
            
            %t = [39106.6080000000 39106.6120000000 39106.6160000000 39106.6200000000 39106.6240000000 39106.6280000000 39106.6320000000 39106.6360000000 39106.6400000000 39106.6440000000 39106.6480000000 39106.6520000000 39106.6560000000 39106.6600000000 39106.6640000000 39106.6680000000 39106.6720000000 39106.6760000000 39106.6800000000 39106.6840000000 39106.6860000000 39106.6900000000 39106.6940000000 39106.6980000000 39106.7020000000 39106.7060000000 39106.7100000000 39106.7140000000 39106.7180000000 39106.7220000000 39106.7260000000 39106.7300000000 39106.7340000000 39106.7380000000 39106.7420000000 39106.7460000000 39106.7500000000 39106.7540000000 39106.7580000000 39106.7620000000 39106.7660000000 39106.7700000000 39106.7740000000 39106.7780000000 39106.7820000000 39106.7860000000 39106.7900000000 39106.7940000000 39106.7980000000 39106.8020000000 39106.8060000000 39106.8100000000 39106.8140000000 39106.8180000000 39106.8220000000 39106.8260000000 39106.8300000000 39106.8340000000 39106.8380000000 39106.8420000000 39106.8460000000 39106.8500000000 39106.8540000000 39106.8580000000 39106.8620000000 39106.8660000000 39106.8700000000 39106.8740000000 39106.8780000000 39106.8820000000 39106.8860000000 39106.8900000000 39106.8940000000 39106.8980000000 39106.9020000000 39106.9060000000 39106.9100000000 39106.9140000000 39106.9180000000 39106.9220000000 39106.9260000000 39106.9300000000 39106.9340000000 39106.9380000000 39106.9420000000 39106.9460000000 39106.9500000000 39106.9540000000 39106.9580000000 39106.9620000000 39106.9660000000 39106.9700000000 39106.9740000000 39106.9780000000 39106.9820000000 39106.9860000000 39106.9900000000 39106.9940000000 39106.9980000000 39107.0020000000 39107.0060000000 39107.0100000000 39107.0140000000 39107.0180000000 39107.0220000000 39107.0260000000 39107.0300000000 39107.0340000000 39107.0380000000 39107.0420000000 39107.0460000000 39107.0500000000 39107.0540000000 39107.0580000000 39107.0620000000 39107.0660000000 39107.0700000000 39107.0740000000 39107.0780000000 39107.0820000000 39107.0860000000 39107.0900000000 39107.0940000000 39107.0980000000 39107.1020000000 39107.1060000000 39107.1100000000 39107.1140000000 39107.1180000000 39107.1220000000 39107.1260000000 39107.1300000000 39107.1340000000 39107.1380000000 39107.1420000000 39107.1460000000 39107.1500000000 39107.1540000000 39107.1580000000 39107.1620000000 39107.1660000000 39107.1700000000 39107.1740000000 39107.1780000000 39107.1820000000 39107.1860000000 39107.1900000000 39107.1940000000 39107.1980000000 39107.2020000000 39107.2060000000 39107.2100000000 39107.2140000000 39107.2180000000 39107.2220000000 39107.2260000000 39107.2300000000 39107.2340000000 39107.2380000000 39107.2420000000 39107.2460000000 39107.2500000000 39107.2540000000 39107.2580000000 39107.2620000000 39107.2660000000 39107.2700000000 39107.2740000000 39107.2780000000 39107.2820000000 39107.2860000000 39107.2900000000 39107.2940000000 39107.2980000000 39107.3020000000 39107.3060000000 39107.3100000000 39107.3140000000 39107.3180000000 39107.3220000000 39107.3260000000 39107.3300000000 39107.3340000000 39107.3380000000 39107.3420000000 39107.3460000000 39107.3500000000 39107.3540000000 39107.3580000000 39107.3620000000 39107.3660000000 39107.3700000000 39107.3740000000 39107.3780000000 39107.3820000000 39107.3860000000 39107.3900000000 39107.3940000000 39107.3980000000 39107.4020000000 39107.4060000000 39107.4100000000 39107.4140000000 39107.4180000000 39107.4220000000 39107.4260000000 39107.4300000000 39107.4340000000 39107.4380000000 39107.4420000000 39107.4460000000 39107.4500000000 39107.4540000000 39107.4580000000 39107.4620000000 39107.4660000000 39107.4700000000 39107.4740000000 39107.4780000000 39107.4820000000 39107.4860000000 39107.4900000000 39107.4940000000 39107.4980000000 39107.5020000000 39107.5060000000 39107.5100000000 39107.5140000000 39107.5180000000 39107.5220000000 39107.5260000000 39107.5300000000 39107.5340000000 39107.5380000000 39107.5420000000 39107.5460000000 39107.5500000000 39107.5540000000 39107.5580000000 39107.5620000000 39107.5660000000 39107.5700000000 39107.5740000000 39107.5780000000 39107.5820000000 39107.5860000000 39107.5900000000 39107.5940000000 39107.5980000000 39107.6020000000 39107.6060000000 39107.6100000000 39107.6140000000 39107.6180000000 39107.6220000000 39107.6260000000 39107.6300000000 39107.6340000000 39107.6380000000 39107.6420000000 39107.6460000000 39107.6500000000 39107.6540000000 39107.6580000000 39107.6620000000 39107.6660000000 39107.6700000000 39107.6740000000 39107.6780000000 39107.6820000000 39107.6860000000 39107.6900000000 39107.6940000000 39107.6980000000 39107.7020000000 39107.7060000000 39107.7100000000 39107.7140000000 39107.7180000000 39107.7220000000 39107.7260000000 39107.7300000000 39107.7340000000 39107.7380000000 39107.7420000000 39107.7460000000 39107.7500000000 39107.7540000000 39107.7580000000 39107.7620000000 39107.7660000000 39107.7700000000 39107.7740000000 39107.7780000000 39107.7820000000 39107.7860000000];

            %pulse_detect_options.method = 'manual'; % Use 'manual' for manual selection or 'load' for loading template
            %pulse_detect_options.min = 0.5; % Threshold for correlation with QRS-wave
            %pulse_detect_options.file = 'posthoc_cpulse.mat'; % Path to save/load template

            % Define verbose options
            %verbose.level = 2; % Level of verbosity (0 for silent, >=2 for debugging plots)
            %verbose.fig_handles = [];

            %[cpulse, verbose] = tapas_physio_get_cardiac_pulses_manual_template(c, t, pulse_detect_options, verbose);
 

            %physio.preproc.cardiac.posthoc_cpulse_select.method = 'manual';
            %physio.preproc.cardiac.posthoc_cpulse_select.output_posthoc_cpulse = 'posthoc_cpulse.mat';
            %physio.preproc.cardiac.posthoc_cpulse_select.load.file = 'posthoc_cpulse.mat';
            
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
        
 
  
