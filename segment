#author @t_sanf


import os
import pydicom
import numpy as np
from PIL import Image

#local file imports
from parsing_VOI import *
from utils import *


class Segment(ParseVOI):
    '''Creates bounding boxes based on polygon segmentation, and forms a three channel .jpg and also T2 for reference
    '''

    def __init__(self,filetype='PIRADS'):
        self.datapath=r'/home/*user*/Desktop/prostateX_processed'

    def segment_allpatients(self):
        '''loop over all patients in root directory, perform segmentation on all of them
        '''

        path_list=[os.path.join(self.datapath,file) for file in os.listdir(self.datapath)]
        for path in path_list:
            __, file = os.path.split(path)
            print("segmenting file for patient {}".format(file))
            try:
                self.savefiles(path)
            except Exception as e:
                print('unable to perform segmentation for patient {}'.format(file))
                if not os.path.exists(os.path.join(path, 'exception_log')):
                    os.mkdir(os.path.join(path, 'exception_log'))
                os.chdir(os.path.join(path, 'exception_log'))
                f = open('segmentation_log.txt', 'w')
                f.write('file %s not able to be downloaded' % e)
                f.close()

    def savefiles(self,patient_dir):
        '''
        alignes arrays for each sequence as 3D arrays, then saves these as .jpg, also saves T2 to make sure segmentation worked
        :param patient_dir: full path to root directomy for each patient
        :return: none - just save files as
        '''

        #get filename for file saving later
        path, file = os.path.split(patient_dir)

        # make a directory if one doesn't already exist to save all images/dicom files
        if not os.path.exists(os.path.join(patient_dir, 'jpg')):
            os.mkdir(os.path.join(patient_dir, 'jpg'))

        if not os.path.exists(os.path.join(patient_dir, 'segmented_t2')):
            os.mkdir(os.path.join(patient_dir, 'segmented_t2'))


        #get list of all paths to voi files
        voi_paths = self.voi_path(patient_dir)

        for voi_path in voi_paths:
            array_dict=self.segment(patient_dir,voi_path)

            for index in array_dict.keys():
                index_dict=array_dict[index]

                # extract each sequance array and combine into numpy array
                t2_array = index_dict['t2'][0].pixel_array
                adc_array = index_dict['adc'][0].pixel_array
                highb_array = index_dict['highb'][0].pixel_array
                stacked_image = np.dstack((t2_array, adc_array, highb_array))

                #extract label from tuple
                label=index_dict['t2'][1]

                #rescale using function from utils
                for i in range(0,3):
                    stacked_image[:,:,i]=rescale_array(stacked_image[:,:,i])

                #save stack as jpg
                img = Image.fromarray(stacked_image, 'RGB')
                img.save(os.path.join(patient_dir, 'jpg', file+ '_' + str(index) + '_' + label + '.jpg'),'JPEG')

                #save t2 as dicom images to check segemtation was correct
                index_dict['t2'][0].save_as(os.path.join(patient_dir,'segmented_t2',file + '_' + str(index) + '_' + label + '.dcm'))


    def segment(self,patient_dir,voi_path, buffer=10):
        '''extract bounding boxes based on polygon segmentation
        :param patient_dir patient directroy
        :param voi_path --> path to voi file, using function in this class
        :param buffer (int) --> number of voxels to buffer around the class
        return dictionaries within dictionaries with the following format {index:{t2:(array,label),adc:(array,label),highb:(array,label)}}
        '''

        bbox_dict = self.BBox_from_position(voi_path) #dictionary of bbox coords, BBox_from_position inherited
        dicom_dict=self.create_dicom_dict(patient_dir) #dictionary of slices

        segmented_index_dict={}
        for index in sorted(int(k) for k in bbox_dict.keys()):
            coords=bbox_dict[str(index)]  #get the coordinates and lablel for this location
            segmented_series_dict = {}
            for series in dicom_dict.keys():
                series_dict=dicom_dict[series]
                ds = pydicom.dcmread(series_dict[str(index)])
                data = ds.pixel_array
                data_downsampled = data[coords[2] - buffer:coords[4] + buffer, coords[1] - buffer:coords[3] + buffer]
                ds.PixelData = data_downsampled.tobytes()
                ds.Rows, ds.Columns = data_downsampled.shape
                segmented_series_dict[series]=(ds,coords[0]) #tuple that includes (pixel array, label)
            segmented_index_dict[index]=segmented_series_dict
        return segmented_index_dict


    def create_dicom_dict(self,patient_dir):
        '''
        creates a dictionary that maps slice number to a filepath for each series
        :param patient_dir full path to patient directory
        :return: all segmented files
        '''
        #obtain path to the directory containing the .voi files and all 3 dicom file series
        dicom_paths={'t2':order_dicom(os.path.join(patient_dir,'dicoms','t2')),
                     'adc':order_dicom(os.path.join(patient_dir,'dicoms','adc','aligned')),
                     'highb':order_dicom(os.path.join(patient_dir,'dicoms','highb','aligned'))}

        #create dict of slices for which
        dicom_dict={}
        for series in dicom_paths.keys():
            series_dict={}
            pathlist=dicom_paths[series]
            for index,path in enumerate(pathlist):
                series_dict[str(index)]=path
            dicom_dict[str(series)]=series_dict
        return(dicom_dict)


    def voi_path(self, patient_dir,filetype='PIRADS'):  #consier converting to regular expression to be more robust
        '''
        function to find all .voi files within a specific patient directory based on
        :param patient_dir (str) full path to patient directory
        :param filetype (str): type of segementation to perform
               PIRADS - tumor
               wp- whole prostate
               u -urethra
        :return: list of files that have the search term of interest (i.e. 'PIRADS')
        '''
        files=os.listdir(os.path.join(patient_dir,'voi'))

        file_path=[]
        for file in files:
            split_file=file.replace(' ','_').split('_')
            if filetype not in split_file:
                continue
            elif filetype in split_file:
                file_path+=[str(os.path.join(patient_dir,'voi',file))]
        return file_path


if __name__=='__main__':
    c=Segment()
    #datapath=os.path.join(c.datapath,'7880984_20180917')
    c.segment_allpatients()

