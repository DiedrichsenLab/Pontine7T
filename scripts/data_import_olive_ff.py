"""
Imports SPM GLM output into the Functional Fusion folder structure.

For every subject this copies
    beta_XXXX.nii  ->  <sub>_<ses>_run-XX_reg-XX_beta.nii
    ResMS.nii      ->  <sub>_<ses>_resms.nii
    mask.nii       ->  <sub>_<ses>_mask.nii
    SPM_info.tsv   ->  <sub>_<ses>_reginfo.tsv   (with standard column names)
    deformation    ->  <sub>_space-MNI152NLin2009cSym_xfm.nii
participants.tsv is not touched.

Run:   python import_olive.py      (only needs numpy and pandas)
"""
import shutil
from pathlib import Path
import numpy as np
import pandas as pd

# ----------------------------- SETTINGS ------------------------------------
base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new'   # Functional Fusion folder (output)
source_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Olive7T'    # your own data (input)
dataset = 'Olive'                           # dataset folder inside base_dir
ses = 'ses-s1'                              # session name used by Functional Fusion
space = 'MNI152NLin2009cSym'                # template of the deformation files

# Functional Fusion subject id -> name of that subject in your own folders
subjects = {
    'sub-07': 'S07',
}

# Where your files are. {src} is replaced by your own subject name (S01, ...)
# when each subject is imported. These are plain strings, NOT f-strings.
glm_dir = source_dir + '/GLM_firstlevel_2/{src}'   # has beta_XXXX.nii, ResMS.nii, mask.nii, SPM_info.tsv
xfm_file = source_dir + '/anatomicals/{src}/{src}_space-MNI152NLin2009cSym_xfm.nii'

# Column in SPM_info.tsv -> standard Functional Fusion column.
# 'run' and 'task_code' are required. 'cond_code', 'half' and 'reg_id' are
# optional: delete the line if you do not have such a column.
col = {
    'run': 'run',
    'taskName': 'task_code',
}
# Column that flags instruction regressors (non-zero = instruction). Those rows
# get task_code 'instrct', which Functional Fusion leaves out when extracting.
# Set to None if you have no such column.
instr_col = 'inst'
keep_other_col = False                  # True: carry the remaining SPM_info columns along
# ---------------------------------------------------------------------------


def make_reginfo(spm_info_file):
    """Reads SPM_info.tsv and returns the reginfo table.
    Row i must describe beta_<i+1>.nii."""
    info = pd.read_csv(spm_info_file, sep='\t')
    missing = [c for c in col if c not in info.columns]
    if missing:
        raise KeyError(f'{spm_info_file} has no column(s) {missing}. '
                       f'Available: {list(info.columns)}')
    # Run numbers must not restart in a second session of the same GLM
    if 'sess' in info.columns and 'run' in col and \
            len(info[['sess', 'run']].drop_duplicates()) > info['run'].nunique():
        raise ValueError(f'{spm_info_file}: run numbers repeat across "sess". '
                         'Runs from different sessions would be merged.')
    D = info.rename(columns=col)
    is_instr = np.zeros(len(D), dtype=bool)
    if instr_col is not None:
        if instr_col not in info.columns:
            raise KeyError(f'{spm_info_file} has no column "{instr_col}"')
        is_instr = info[instr_col].fillna(0).astype(float).values != 0
    if not keep_other_col:
        D = D[list(col.values())]
    for required in ['run', 'task_code']:
        if required not in D.columns:
            raise KeyError(f'col must map one of your columns to "{required}"')

    D['run'] = D['run'].astype(float).astype(int)
    if 'reg_id' not in D.columns:
        # position of the regressor within its run: 0, 1, 2, ...
        D['reg_id'] = D.groupby('run').cumcount()
    D['reg_id'] = D['reg_id'].astype(float).astype(int)
    if 'half' not in D.columns:
        # odd runs = 1, even runs = 2
        D['half'] = 2 - (D['run'] % 2)
    if 'cond_code' not in D.columns:
        # required by the extraction; empty entries are treated as 'task'
        D['cond_code'] = np.nan
    for c in ['task_code', 'cond_code']:
        D[c] = D[c].astype('string').str.strip()
    D.loc[is_instr, 'task_code'] = 'instrct'

    # The same task / condition twice in one run would be averaged on extraction
    task = D[~is_instr]
    dup = task[task.duplicated(['run', 'task_code', 'cond_code'], keep=False)]
    if len(dup) > 0:
        print(f'WARNING {spm_info_file}: these appear more than once per run and '
              f'will be averaged: {sorted(dup.task_code.unique())}')

    if D.duplicated(['run', 'reg_id']).any():
        raise ValueError('run / reg_id combinations are not unique')
    first = ['run', 'half', 'reg_id', 'task_code', 'cond_code']
    return D[first + [c for c in D.columns if c not in first]]


def copy(src, dest):
    if not Path(src).exists():
        raise FileNotFoundError(f'Missing file: {src}')
    shutil.copyfile(src, dest)


def import_subject(sub, src_name):
    src_dir = Path(glm_dir.format(src=src_name))
    ff_dir = Path(base_dir) / dataset / 'derivatives' / 'ffimport' / sub
    func_dir = ff_dir / 'func' / ses
    anat_dir = ff_dir / 'anat'
    func_dir.mkdir(parents=True, exist_ok=True)
    anat_dir.mkdir(parents=True, exist_ok=True)

    # Regressor information
    D = make_reginfo(src_dir / 'SPM_info.tsv')
    D.to_csv(func_dir / f'{sub}_{ses}_reginfo.tsv', sep='\t', index=False)

    # Betas: row i of SPM_info <-> beta_<i+1>.nii
    for i, row in enumerate(D.itertuples()):
        copy(src_dir / f'beta_{i + 1:04d}.nii',
             func_dir / f'{sub}_{ses}_run-{row.run:02d}_reg-{row.reg_id:02d}_beta.nii')

    # Residual variance, mask and deformation
    copy(src_dir / 'ResMS.nii', func_dir / f'{sub}_{ses}_resms.nii')
    copy(src_dir / 'mask.nii', func_dir / f'{sub}_{ses}_mask.nii')
    copy(xfm_file.format(src=src_name), anat_dir / f'{sub}_space-{space}_xfm.nii')

    # Report
    n_beta = len(list(src_dir.glob('beta_*.nii')))
    print(f'{sub}: {len(D)} betas imported, {D.run.nunique()} runs, '
          f'{n_beta - len(D)} beta files in the GLM folder not listed in SPM_info '
          f'(intercepts / nuisance regressors)')
    if n_beta < len(D):
        raise ValueError(f'{sub}: SPM_info has more rows than there are beta files')


if __name__ == '__main__':
    for sub, src_name in subjects.items():
        import_subject(sub, src_name)