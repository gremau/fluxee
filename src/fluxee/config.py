import os
import yaml

# Get the project configuration
yamlfile = os.path.join("fluxee_config", "fluxee_conf.yaml")
if os.path.isfile(yamlfile):
    print("Loading fluxee configuration file {0}".format(yamlfile))
    stream = open(yamlfile, 'r')
    conf = yaml.load(stream)
else:
    raise ValueError('Project YAML file missing from default location')

# Project name
projectname = conf['projectname']
print('Configuration for project {0}'.format(projectname))

# Get valid paths from configuration
# NOTE - If base_path is unset in conf, all other paths must be complete
base_path = conf['paths']['base']

sitedata_file = os.path.join(base_path,
        conf['paths']['sitedata'])
rawdata_incoming_path = os.path.join(base_path,
        conf['paths']['rawdata_incoming'])
rawdata_backup_path = os.path.join(base_path,
        conf['paths']['rawdata_backup'])
qadata_path = os.path.join(base_path,
        conf['paths']['qadata'])
fluxee_py_path = os.path.join(base_path,
        conf['paths']['fluxee_py'])

# Get name and subdirectory for fluxee data types
datasubdirs = {'rawdata_incoming':'',  
    'rawdata_backup': 'raw_bak/',
    'rawdata_standardized': 'raw_std/',
    'quality_assured': 'qa/'}

datadirs = list(datasubdirs.keys())

print('Each site has {0} data directories available: \n {1} \n'
        'Use "iodat.site_datadir(sitename, datadir=datadir_name)"'
        'to get proper path'.format(len(datadirs), ', '.join(datadirs)))

