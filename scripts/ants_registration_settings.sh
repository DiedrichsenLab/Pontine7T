SYNCC: 

antsRegistration 
-d 3 # dimension
-r [0x35bdc87b0,0x35bdd1680,1]  #
-m mattes[0x35bdc87b0,0x35bdd1680,1,32,regular,0.2] 
-t Rigid[1] 
-c 2100x1200x1200x0 # convergence
-s 3x2x1x0  # smoothing factors
-f 4x4x2x1  # shrink factors
-x [NA,NA] # Mask   
-m mattes[0x35bdc87b0,0x35bdd1680,1,32,regular,0.2] 
-t Affine[1] 
-c 1200x1200x100 
-s 2x1x0 
-f 4x2x1 
-x [NA,NA] # Mask   
-m CC[0x35bdc87b0,0x35bdd1680,1,4] 
-t SyN[0.15,3,0] 
-c [2100x1200x1200x20,1e-7,8] 
-s 3x2x1x0 
-f 4x3x2x1 
-u 1 
-z 1 
-o [/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/bold_normalization/transforms/MNI2009c_T2_CC/xfm_S98_,0x35bde6010,0x35bdffff0] 
-x [NA,NA] 
--float 1 
--write-composite-transform 0 
-v 1

SyNBold
antsRegistration 
-d 3 
-r [0x3519988a0,0x3519c08f0,1] 
-m mattes[0x3519988a0,0x3519c08f0,1,32,regular,0.2] 
-t Rigid[0.25] 
-c [1200x1200x100,1e-6,5] 
-s 2x1x0 
-f 4x2x1 -
x [NA,NA] -m mattes[0x3519988a0,0x3519c08f0,1,32] 
-t SyN[0.200000,3.000000,0.000000] 
-c [40x20x0,1e-7,8] 
-s 2x1x0 
-f 4x2x1 
-u 1 
-z 1 
-o [/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/bold_normalization/transforms/MNI2009c_T1bold/xfm_S95_,0x35199ecf0,0x3519a6080] -x [NA,NA] --float 1 --write-composite-transform 0 -v 1
