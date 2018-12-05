# MP_MRI
Miscellaneous tools for working with multiparametric MRI

create_imagefolder.py
command line program that puts all your files in Pytorch ImageFolder.  Expects .jpeg images with label after final underscore
flags:
--input_path (str)# full path to directory with images
--output_path (str) # full path to directory you want your ImageFolder
--percent(float) #what percent of total patients allocated for validation

Example Usage:
1) Download file and navigate to folder in command line
2) type:
python3 create_imagefolder.py --input_path '/path_to_folder_with_images' --output_path '/path_to_output_folder' --percent=0.2
