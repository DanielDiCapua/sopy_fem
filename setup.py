from distutils.core import setup
import setuptools

setup(
  name = 'sopy_fem',
  #packages = ['sopy_fem'], # this must be the same as the name above
  packages= setuptools.find_packages(),
  version = '0.1.2',
  description = 'An academic program for finite elments analysis of 2D solids and structures',
  author = 'Daniel Di Capua',
  author_email = 'dicapua67@gmail.com',
  url = 'https://github.com/DanielDiCapua/sopy_fem.git', # use the URL to the github repo
  download_url = 'https://github.com/DanielDiCapua/sopy_fem.git/tarball/0.1.2',
  keywords=['developing', 'example'],
  # data_files=[('Examples',
  #   ['sopy_fem/Examples/dynamics_TRUSS02/data.json',
  #   'sopy_fem/Examples/electrical_BR02/data.json',
  #   'sopy_fem/Examples/mechanics_BR02/data.json',
  #   'sopy_fem/Examples/mechanics_QU04/data.json',
  #   'sopy_fem/Examples/mechanics_TR03/data.json',
  #   'sopy_fem/Examples/structural_TRUSS02/data.json',
  #   'sopy_fem/Examples/thermal_BR02/data.json'])],
  include_package_data=True,
  package_dir={'sopy_fem': 'sopy_fem/sopy_fem'},
  package_data={'sopy_fem': [
    'Examples/dynamics_TRUSS02/data.json',
    'Examples/electrical_BR02/data.json',
    'Examples/mechanics_BR02/data.json',
    'Examples/mechanics_QU04/data.json',
    'Examples/mechanics_TR03/data.json',
    'Examples/structural_TRUSS02/data.json',
    'Examples/thermal_BR02/data.json'
    ]},
  classifiers = [],
)