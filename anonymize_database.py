import pandas as pd
import os
import random
import shutil
import pydicom
import win32com.client
from io import StringIO
from docx import Document
import string

#local import
from utils import *


class Anonymize:

    def __init__(self):
        self.phi_database=r''
        self.anonymize_database=r''
        self.transfer_database=r''
        self.anonymize_database_key=r''
        self.database=

    def anonymize_ID_transfer(self):
        '''
        copy files from network drive to new folder, anonymize the name of the patient
        :param dir: directory in which all files that need anonymization
        '''

        #read in key:value dataframe, create new one if there isn't already one available
        if os.path.exists(os.path.join(self.anonymize_database_key,'key_file.csv')):
            keys=pd.read_csv(os.path.join(self.anonymize_database_key,'key_file.csv'))
            keys=keys.drop(columns=['Unnamed: 0'])
        else:
            keys=pd.DataFrame(columns=['anonymize_id','original_id'])
            keys.to_csv(os.path.join(self.anonymize_database_key,'key_file.csv'))

        databases=self.database

        dictionary={}
        random_num_list=[]
        for database in databases:

            dir=os.path.join(self.phi_database,database)
            for file in os.listdir(dir):
                print(file)
                random_num=random.randint(1000000000,9999999999)
                print(random_num)
                if str(random_num) not in os.listdir(os.path.join(self.anonymize_database,database)) and file not in keys['original_id'].tolist():
                    dictionary[random_num]=file
                    random_num_list+=[str(random_num)]
                    shutil.copytree(os.path.join(self.phi_database,database,file),os.path.join(self.anonymize_database,database,str(random_num)))
                else:
                    print("file {} already in this location".format(file))
        database=pd.DataFrame.from_dict(dictionary,orient='index')
        if database.empty == False:
            database.reset_index(inplace=True)
            database.columns=['anonymize_id','original_id']
            database=pd.concat([keys,database])
            database.to_csv(os.path.join(self.anonymize_database_key,'key_file.csv'))


    def anonomyze_word_docs(self):
        '''
        from word document, creates .txt file with raddiology report
        :param dir: directory in which all files that need anonymization
        :return:
        '''

        databases=self.database

        exception_logger=[]
        for database in databases:
            for dir in self.check_for_anonymized_word_docs():
                print('anonymizing files for patient {}'.format(dir))
                #define directory, word file path, and wordx file path
                reports_dir=os.path.join(self.anonymize_database,database,dir,'radiologyreports')

                try:

                    #define length files
                    word_file = os.path.join(reports_dir, os.listdir(reports_dir)[0])
                    docx_file = word_file + 'x'

                    # first convert to docx to play more
                    save_doc_as_docx(word_file,docx_file) #note - this function is from the utils folder
                    output=getText_without_first_line(docx_file)  #note - this function is from the utils folder

                    #save output as .txt file
                    os.chdir(reports_dir)
                    text_file = open("radiology_report.txt", "w")
                    text_file.write(output)
                    text_file.close()

                    #remove dates from .txt file
                    remove_dates("radiology_report.txt") #from .utils

                    #remove all files that are not the anonymized file
                    for file in os.listdir(reports_dir):
                        if file !='radiology_report.txt':
                            os.remove(file)

                except:
                    print("patient {} not able to be processed".format(dir))
                    exception_logger+=[dir]
        print('error for the following patient {}'.format(exception_logger))
        return exception_logger


    def check_for_anonymized_word_docs(self):
        '''
        iterate through all databases and check to see which have already been anonymized
        :return: list of all patients that need updating
        '''

        databases = self.database

        needs_update=[]
        for database in databases:
            for patient_dir in os.listdir(os.path.join(self.anonymize_database, database)):
                reports_dir = os.path.join(self.anonymize_database, database, patient_dir, 'radiologyreports')
                if os.path.exists(reports_dir):
                    reports_list=os.listdir(reports_dir)

                    if "radiology_report.txt" in reports_list:
                        print('patient {} all G!!!!!!!!!!!!!!!!!!'.format(patient_dir))

                    if not "radiology_report.txt" in reports_list:
                        print('patient {} needs updating'.format(patient_dir))
                        needs_update+=[patient_dir]
        return sorted(needs_update)

    def transfer_files(self):
        '''transfer anonymized files to new database.  Also remove dates from transferred file so completely anonymized'''
        
        databases = self.database
        
        for database in databases:
            for patient_dir in os.listdir(os.path.join(self.anonymize_database, database)):
                reports_dir = os.path.join(self.anonymize_database, database, patient_dir, 'radiologyreports')
                if os.path.exists(reports_dir):
                    reports_list=os.listdir(reports_dir)

                    if "radiology_report.txt" in reports_list:
                        original_file=os.path.join(reports_dir,'radiology_report.txt')
                        new_file=os.path.join(self.transfer_database,patient_dir+'.txt')
                        shutil.copy2(original_file,new_file)
                        remove_dates(new_file)   #remove dates function from utils
                    else:
                        print("patient {} radiology report not transferred".format(patient_dir))


if __name__=='__main__':
    c=Anonymize()
    c.transfer_files()


