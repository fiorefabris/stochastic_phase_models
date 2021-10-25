from setuptools import setup

setup(
   name='adler',
   version='1.0',
   description='module for studying noisy adler statistics',
   author='fiorefabris',
   author_email='fiorefabris@gmail.com.com',
   packages=['adler'],  #same as name
   install_requires=['matplotlib', 'pandas','numpy'], #external packages as dependencies
   packages=['adler', 'adler.get_time_series','adler.plotting','adler.pulse_detection'],
)
#multiprocessing,,'itertools','functools','pickle','math','time'
#scripts
