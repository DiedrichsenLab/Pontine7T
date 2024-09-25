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
