from setuptools import setup

setup(
    name='rtergpy',
    version='0.2.2',
    description='A package for real-time processing of Earthquake Energy and duration',
    url='https://github.com/avnewman/',
    author='Andrew Newman',
    author_email='anewman@gatech.edu',
    license='MIT',
    packages=['rtergpy'],
    install_requirements=['pygmt>=0.3.1', 
                        'numpy>=1.0',
                        'matplotlib>=3.3.0',
                        'tqdm',
                        'obspy>=1.2.0',
                        'pandas>=1.2.0'
                        ],
    classifiers=[ 
        'Development Status :: 2 - Beta',
        'Intended Audience :: Science',
        'Operating System :: OS Independent',
        'Programming Language :: Python :3.9',
        ],
)
