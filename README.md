Miscellaneous tools 

# create_imagefolder.py
command line program that puts all your files in Pytorch ImageFolder.  Expects .jpeg images with label after final underscore

flags:

--input_path # full path to directory with images

--output_path # full path to directory you want your ImageFolder 

--percent(float) #what percent of total patients allocated for validation 

Example Usage:
1) Download file and navigate to folder in command line
2) type:
python3 create_imagefolder.py --input_path '/path_to_folder_with_images' --output_path '/path_to_output_folder' --percent=0.2

# torchvision_train.py
command line program that allows you to train neural networks using transfer learning in pytorch.  Basically just takes the 
examples from pytorch website and allows them to run with command line

flags:

--image_dir        image directory in ImageFolder format (can use create_imagefolder script above to create

--out_dir          filepath to where you want your model saved

--mn (path)        name of model.  Options include resnet, resnet18 - resnet 152, alexnet, vgg, squeezenet.  Default resnet18

--nc (int)         num classes.  default 2

--num_epochs (int) number of training iterations, default 1

--input_size (iut) pixel size of input image.  Default is 224, larger images will be cropped to this size

--feat_ext (bool)  if false, finetune the entire model.  Otherwise, only reshaped layers

--pretrain(bool)   if false, no transfer learning used.  May the force be with you.
