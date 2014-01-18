#!/usr/bin/env python

import os, subprocess

if __name__ == '__main__':
    hadoop_dir = '/usr/class/cs149/hadoop-1.1.1'
    hadoop_bin_dir = os.path.join(hadoop_dir, 'bin')
    hadoop_exe = os.path.join(hadoop_bin_dir, 'hadoop')

    hadoop_conf_dir = os.path.abspath('conf')

    subprocess.check_call(
        ['./hadoop-daemon.sh', '--config', hadoop_conf_dir, 'stop', 'namenode'],
        cwd = hadoop_bin_dir)
    subprocess.check_call(
        ['./hadoop-daemon.sh', '--config', hadoop_conf_dir, 'stop', 'datanode'],
        cwd = hadoop_bin_dir)
