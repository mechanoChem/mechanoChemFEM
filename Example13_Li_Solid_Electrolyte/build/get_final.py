#!/usr/bin/env python3

import os

def exe_cmd(cmd, output=False):
    import subprocess, os
    if output:
        output_info = os.popen(cmd).read()
        return output_info
    else:
        os.system(cmd)

def get_frame_number(str0):
    str1 = str0.split('.vtk')[0]
    str2 = str1.split('-')[1]
    #print(str2)
    return int(str2)

cmd = 'find ./ | grep \.vtk'
ofile = exe_cmd(cmd, True)

file_list = ofile.strip().split()



folder_list = {}
for f0 in file_list:
    #print(f0)
    name_list = f0.split('/')
    _folder = ('/').join(name_list[0:-1]) + '/'
    _file = name_list[-1]
    #print(_folder, _file)
    try:

        if get_frame_number(folder_list[_folder]) < get_frame_number(_file):
            folder_list[_folder] = _file
    except:
        folder_list[_folder] = _file
    print('dict: folder=', _folder,  ' val=', folder_list[_folder])
    #folder_list
#     if f0.find('0.vtk') > 0:
#         print('keep', f0)
#     else:
#         if f0.find('.vtk') > 0:
#             cmd = 'rm ' + f0
#             print(cmd)
#             exe_cmd(cmd)
cmd = 'rm -rf final/*'
exe_cmd(cmd)
cmd = 'mkdir -p final'
exe_cmd(cmd)
for key, val in folder_list.items():
    cmd = 'cp ' + key + '/' + val + ' final/' + key.replace('./', '').replace('/','_') + val
    exe_cmd(cmd)
    #print(cmd)
    #print(key, val)
