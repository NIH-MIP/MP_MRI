#author: @t_sanf


import os
import pandas as pd
from pathlib import Path
import subprocess

class Align:

    def __init__(self):
        self.alignmentPATH=r'path to exacutable'
        self.basePATH=Path(r'')
        self.failsaveFILE = r''

    def align_files(self):
        '''
        file to iterate over all files in robots_only folder
        :return: nothing specific returned, but files are aligned
        '''
        list_of_files = self.check_empty_files()
        files=self.check_aligned_files(list_of_files)
        align_fails=[]
        for file in files:
            print("aligning file {}".format(file))
            path=str(self.basePATH)+r"\\"+file+r"\\"+r"dicoms"+r"\\"
            print(self.alignmentPATH+path)
            try:
                os.system(self.alignmentPATH+path)
            except:
                print("cannot do this file")
                align_fails+=[file]
        try:
            align_fails=pd.Series(align_fails)
            align_fails.to_csv(os.path.join(self.failsaveFILE,'alignment_fails.csv'))
        except:
            print("TOM! Same Fregin error. FIX IT! ")
            pass

    def check_aligned_files(self,list_of_files):
        '''
        check for folder called 'align' in adc/highb
        :return: aligned files
        '''
        #for development only
        aligned_files=[]
        for file in list_of_files:
            adc_output=os.listdir(os.path.join(self.basePATH,file,'dicoms','adc'))
            highb_output=os.listdir(os.path.join(self.basePATH,file,'dicoms','highb'))
            if not 'aligned' in adc_output and not 'aligned' in highb_output:
                aligned_files+=[file]
        print("there are a total of {} files in need of alignment".format(len(aligned_files)))
        return(aligned_files)


    def check_empty_files(self):
        '''
        check for file completion in robert_huang and return only completed files
        :return: list of patient names
        '''

        path = self.basePATH
        files = os.listdir(path)
        print("searching total of {} files in robots_only folder to see which ones need aligning".format(len(files)))
        completed=[]
        for file in files:
            if len(file.split('_')) == 2:
                voi_dir = os.path.join(path, file, 'voi')
                adc_dir = os.path.join(path, file, r'dicoms\adc')
                highb_dir = os.path.join(path, file, r'dicoms\highb')
                t2_dir = os.path.join(path, file, r'dicoms\t2')

                voi_files = os.listdir(voi_dir)
                adc_files = os.listdir(adc_dir)
                highb_files = os.listdir(highb_dir)
                t2_files = os.listdir(t2_dir)

                if len(voi_files) > 0 and len(adc_files) > 5 and len(highb_files) > 5 and len(t2_files) > 5:
                    completed += [file]

        print('total of {} complete files'.format(len(completed)))
        return(completed)


if __name__=="__main__":
    c=Align()
    c.align_files()




