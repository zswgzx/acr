#  ACR Phantom QA
MATLAB codes for [ACR phantom](https://www.acraccreditation.org/-/media/ACRAccreditation/Documents/MRI/LargePhantomGuidance.pdf) QA, work in progress

Tests are:
- geometry accuracy 			(geom.m, done)
- high contrast spatial resolution	(spatial_res.m)
- slice thickness accuracy 		(slice_th.m)
- slice position accuracy 		(slice_pos.m)
- image intensity uniformity		(iiu_psg.m)
- percentage signal ghosting		(iiu_psg.m)
- low contrast object detectability 	(lcod.m)

otsu.m is used for binary segmentation