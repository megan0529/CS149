#!/usr/bin/env python

import getpass, os, random, shutil, socket, subprocess, time

def copy_with_substitutions(from_filename, to_filename, substitutions):
    with open(to_filename, 'wb') as to_file:
        with open(from_filename, 'rb') as from_file:
            to_file.write(from_file.read().format(**substitutions))

def copy_dir_with_substitutions(from_dir, to_dir, substitutions):
    filenames = os.listdir(from_dir)
    for filename in filenames:
        copy_with_substitutions(
            os.path.join(from_dir, filename),
            os.path.join(to_dir, filename),
            substitutions)

if __name__ == '__main__':
    user = getpass.getuser()

    java_home_dir = '/usr/lib/jvm/java-1.6.0-openjdk-amd64'

    hadoop_dir = '/usr/class/cs149/hadoop-1.1.1'
    hadoop_bin_dir = os.path.join(hadoop_dir, 'bin')
    hadoop_exe = os.path.join(hadoop_bin_dir, 'hadoop')

    pristine_conf_dir = os.path.join(hadoop_dir, 'conf')
    template_conf_dir = '/usr/class/cs149/assignments/pa3/local-hadoop/conf'
    user_conf_dir = os.path.abspath('conf')

    hadoop_tmp_dir = '/tmp/hadoop-%s' % user
    hadoop_namenode_pid = '/tmp/hadoop-%s-namenode.pid' % user
    hadoop_datanode_pid = '/tmp/hadoop-%s-datanode.pid' % user

    # Stop any running daemons
    if os.path.exists(hadoop_namenode_pid):
        subprocess.check_call(
            ['./hadoop-daemon.sh', '--config', user_conf_dir, 'stop', 'namenode'],
            cwd = hadoop_bin_dir)
    if os.path.exists(hadoop_datanode_pid):
        subprocess.check_call(
            ['./hadoop-daemon.sh', '--config', user_conf_dir, 'stop', 'datanode'],
            cwd = hadoop_bin_dir)

    # Generate new configuration
    hdfs_port = random.randint(20000, 50000)
    substitutions = {
        'host': socket.gethostname(),
        'hdfs_port': hdfs_port,
        }
    if os.path.exists(user_conf_dir):
        shutil.rmtree(user_conf_dir)
    shutil.copytree(pristine_conf_dir, user_conf_dir)
    copy_dir_with_substitutions(
        template_conf_dir,
        user_conf_dir,
        substitutions)

    # Start daemons
    os.environ['USER_CONF_DIR'] = user_conf_dir
    os.environ['JAVA_HOME'] = java_home_dir
    if os.path.exists(hadoop_tmp_dir):
        shutil.rmtree(hadoop_tmp_dir)

    subprocess.check_call([hadoop_exe, 'namenode', '-format'])
    subprocess.check_call(
        ['./hadoop-daemon.sh', '--config', user_conf_dir, 'start', 'namenode'],
        cwd = hadoop_bin_dir)
    subprocess.check_call(
        ['./hadoop-daemon.sh', '--config', user_conf_dir, 'start', 'datanode'],
        cwd = hadoop_bin_dir)
    subprocess.check_call(
        [hadoop_exe, '--config', user_conf_dir, 'fs', '-mkdir', '/user'])
    subprocess.check_call(
        [hadoop_exe, '--config', user_conf_dir, 'fs', '-mkdir', '/user/%s' % user])
