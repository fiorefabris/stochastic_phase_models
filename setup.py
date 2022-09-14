from setuptools import setup

setup(
   name='adler',
   version='1.22968221223121111111111111',
   description='module for studying noisy adler statistics',
   author='fiorefabris',
   author_email='fiorefabris@gmail.com.com',
   packages=['adler', 'adler.get_time_series','adler.plotting','adler.pulse_detection'],
   install_requires=['matplotlib', 'pandas','numpy'], #external packages as dependencies
)
#multiprocessing,,'itertools','functools','pickle','math','time'
#scripts
