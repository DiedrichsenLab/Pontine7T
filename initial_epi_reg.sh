sub=96
output_dir=/Users/CN/Documents/Projects/7T_Pontine/derivatives/sub-${sub}

# --- Brain extract T1 ---
# Bet the T1 (needs AFNI & FSL - did this on the fmrib cluster for now, but could set it up on the robarts server)
SCRIPTDIR=/Users/CN/Documents/Projects/7T_Pontine/code/analysis/Pontine7T
cd /Volumes/diedrichsen_data$/data/Pontine7T/anatomicals/S${sub}/
sh ${SCRIPTDIR}/optiBET.sh -i anatomical.nii


# --- Register functional data to structural using epi_reg ---
# Set paramters
output_dir=/Users/CN/Documents/Projects/7T_Pontine/derivatives/sub-${sub}
epi=/Volumes/diedrichsen_data$/data/Pontine7T/imaging_data/S${sub}/rmeanrun_01.nii
struct=/Volumes/diedrichsen_data$/data/Pontine7T/anatomicals/S${sub}/anatomical.nii
struct_betted=${output_dir}/anatomical_optiBET_brain.nii.gz
wmseg=/Volumes/diedrichsen_data$/data/Pontine7T/anatomicals/S${sub}/c2anatomical.nii

epi_reg \
--wmseg=${wmseg} \
--epi=${epi} \
--t1=${struct} \
--t1brain=${struct_betted} \
--out=${output_dir}/func_to_struct

# --- Inspect result ---
fsleyes ${struct} -dr 0 3000 \
${epi} -dr 0 15000 \
${output_dir}/func_to_struct.nii.gz -dr 0 15000 &

