This is a Matlab code for verifying the JPEG copy-move constraint, implementing the method proposed in 
 
*Adam Novozámský, Michal Šorel, "Detection of image modification using JPEG compression model", Forensic Science International, 2017*

*Adam Novozámský, Michal Šorel, "JPEG Compression Model in Copy-move Forgery Detection", 2017 (conference version)*
 
This code is mainly for review purposes. A more detailed version will follow soon. 
 
**Files**
 
*copymove_constraint.p* - main function [](see help included in the file)
 
*test_jpegcopymove_random_selection.m* - The code showing how to test one random patch
	against several other random patches.

*test_jpegcopymove_tampering.m* - The code tampers the input image by copying one patch
to a different position, saving the image with given JPG quality and after reading the tapered 
image, the patch is compared with respect to the original position. 



