#!/usr/bin/env python3
import os
import pandas as pd
import shutil
import argparse
import re


class CreateImageFolder:
    '''
    note - this class takes an input path and sorts the
    expects files with .jpeg output and label after final underscore
    '''

    def __init__(self,path_to_input='',path_to_output='',percent=0.2,neg_lab='neg',pos_lab='pos'):
        '''
        Args:
            path_to_input: path to the folder containing images.  Expect them to have image extension (i.e. .jpg) and label after final underscore
            path_to_output: the path where you want your image folder saved
            percent: percent of patients in validation (if you have 100 patients and this value is 0.2, you will have 20 validation patients)
        '''

        self.inpath=path_to_input
        self.outpath=path_to_output
        self.percent=percent
        self.neg_label=neg_lab
        self.pos_label=pos_lab

    def create_image_folder(self):
        '''creates an image folder'''

        #use function to below to make sure all files present in output directory
        self.make_out_dir(neg_lab=self.neg_label,pos_lab=self.pos_label)


        #get all file paths
        all_file_path=pd.Series([os.path.join(self.inpath,file) for file in os.listdir(self.inpath)])

        #split into training and test tests, place into dictionary
        valid=all_file_path.sample(frac=self.percent, random_state=200)
        train=all_file_path.drop(valid.index)
        split={'train':train.tolist(),'valid':valid.tolist()}

        #iteratve over train/valid data and copy files into negative and positive based on filename
        for state in split.keys():
            state_item=split[state]
            for path in state_item:
                neg=re.compile(self.neg_label); pos=re.compile(self.pos_label)
                if neg.search(os.path.basename(path)):
                    shutil.copy2(path,os.path.join(self.outpath,state,self.neg_label))
                if pos.search(os.path.basename(path)):
                    shutil.copy2(path, os.path.join(self.outpath, state, self.pos_label))


    def make_out_dir(self,neg_lab,pos_lab):
        '''
        looks through all the values in the output directory and makes the correct data structure and makes folds
        for all relevant files.
        Returns:
        '''

        #define all directories
        states=['train','valid']

        for state in states:
            path_state_dir=os.path.join(self.outpath,state)
            if not os.path.exists(path_state_dir):
                os.mkdir(path_state_dir)
                os.chdir(path_state_dir)
                if not os.path.exists(os.path.join(path_state_dir,neg_lab)):
                    os.mkdir(os.path.join(path_state_dir,neg_lab))
                if not os.path.exists(os.path.join(path_state_dir,pos_lab)):
                    os.mkdir(os.path.join(path_state_dir,pos_lab))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Create ImageFolder')
    parser.add_argument('--input_path')
    parser.add_argument('--output_path')
    parser.add_argument('--percent', default=0.2, type=float)
    parser.add_argument('--neg_lab', default='neg', type=str)
    parser.add_argument('--pos_lab',default='pos', type=str)
    args = parser.parse_args()

    c=CreateImageFolder(args.input_path,args.output_path,args.percent,args.neg_lab,args.pos_lab)
    c.create_image_folder()
