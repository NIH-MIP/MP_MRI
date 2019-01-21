#author @t_sanf
import os
import pydicom
from sklearn.preprocessing import MinMaxScaler
from collections import OrderedDict
import docx
import win32com.client
from io import StringIO
import re



def rescale_array(array):
    '''rescales the array from 0-250'''
    scaler =MinMaxScaler(feature_range=(0 ,250))
    scaler =scaler.fit(array)
    X_scaled =scaler.transform(array)
    return (X_scaled)


def order_dicom(dicom_dir):
    '''
    sorts the list of filenames based on z position of dicom files
    :param dicom_dir full path to directory with dicom files
    :return list of files in correct order
    '''

    dicoms={}
    for path in [os.path.join(dicom_dir,file) for file in os.listdir(dicom_dir) if file!='VERSION']:
        ds=pydicom.read_file(path)
        dicoms[path] = float(ds.SliceLocation)
    updated_imagelist=[key for (key, value) in sorted(dicoms.items(), key=lambda x: x[1])]
    return(updated_imagelist)

#########################
# function to manipulate word documents

def getText_without_first_line(filename):
    doc = docx.Document(filename)
    fullText = []
    for para in doc.paragraphs[1:]: #starting on 9th line because dates live in the header
        fullText.append(para.text)
    return '\n'.join(fullText)

def save_doc_as_docx(doc_path,docx_path):
    '''
    hackery to deal with stoopd .doc file.  First convert to .docx
    :param doc_path:
    :param docx_path:
    :return:
    '''
    word = win32com.client.gencache.EnsureDispatch('Word.Application')
    word.Visible = False
    wb = word.Documents.Open(doc_path)
    doc = word.ActiveDocument
    doc.SaveAs2(docx_path, FileFormat=16)
    doc.Close()


######################
#regular expression methods

def remove_dates(path_to_txt_file):
    '''uses regular expression to remove all the dates from a text file'''
    f = open(path_to_txt_file, "r")
    lines = f.readlines()
    f.close()
    f = open(path_to_txt_file, "w")
    for line in lines:
        line=re.sub('[0-9]{1,2}[\/][0-9]{1,2}[\/][0-9]{2,4}','date removed', line)
        line = re.sub('Baris Turkbey', ' name removed', line)
        line = re.sub('Choyke, Peter L.', ' name removed', line)
        line = re.sub('15-C-0124', '', line)
        f.write(line)
    f.close
