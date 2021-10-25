from setuptools import setup

setup(
   name='adler',
   version='1.0',
   description='module for studying noisy adler statistics',
   author='fiorefabris',
   author_email='fiorefabris@gmail.com.com',
   packages=['adler'],  #same as name
   install_requires=['matplotlib', 'pandas','numpy','itertools','functools','pickle','math','time'], #external packages as dependencies
   scripts=[]
)
#multiprocessing
