run:


for compiling mex file for edgelets detection
>> mex mx_lines


for compiling mex file for diamond space accumulation
>> mex mx_raster_space


for diamond space accumulation and detection    
>> R = diamond_vanish(InputImg, Normalization, SpaceSize, PatchSize, VanishNumber)
 input:
   InputImg - input image 
   Normalization - float number; normalization of the image (1 means normalization form -1 to 1)
   SpaceSize - int number; resolution of the accumulation space (final space has dims SpaceSize x SpaceSize)
   PatchSize - int row vector; radius of a patch, from which edge pixels are extracted and use for ellipse fitting for detection of orientation of the edge point
   VanishNumber - number of vanishing points to be detected
 output:
   R - structure with fields
    R.Space - accumulated diamond space (further is used for orthogonalization)
    R.PC_VanP - positions of the maxima in R.Space
    R.PC_VanP_Norm - normalized position of the maxima (R.Space bounded from -1 to 1)
    R.CC_VanP - position of the vanishing point in the input image coordinates


for evaluation on a dataset	
>> R = run_on_dataset(DB, DB_path, Normalization, SpaceSize, PatchSize)
 input:
   DB - cell array with images names (e.g. Manhattan_Image_DB_Names)
   DB_path - path to dataset data (with camera parameters and image data)
   Normalization - same as above (if empty, default paramters are used, Normalization = 0.4)
   SpaceSize - same as above (if empty, default paramters are used, SpaceSize = 321)
   PatchSize - same as above (if empty, default paramters are used, PatchSize = [6:4:22])
   (default parameters are from training on a ECCV_TrainingAndTestImageNumbers -> trainingSetIndex)	
 output:
   R - structure with fields:
    R.Image - image name
   R.Detected - 3x3 matrix with three detected vanishing points, each in one column
     R.Orthogonal - 3x3 matrix with three vanishing points after orthogonalization, each in one column