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
        
        