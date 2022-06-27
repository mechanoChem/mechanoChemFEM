#!/usr/bin/env python3

import os

def exe_cmd(cmd, output=False):
    import subprocess, os
    if output:
        output_info = os.popen(cmd).read()
        return output_info
    else:
        os.system(cmd)


cmd = 'find ./ | grep \.vtk'
ofile = exe_cmd(cmd, True)

file_list = ofile.strip().split()

for f0 in file_list:
    if f0.find('0.vtk') > 0:
        print('keep', f0)
    else:
        if f0.find('.vtk') > 0:
            cmd = 'rm ' + f0
            print(cmd)
            exe_cmd(cmd)