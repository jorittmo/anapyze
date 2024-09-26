from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
        name='anapyze_processing',
        version='0.1.0',
        author='Jesus Silva',
        author_email='jesus.bubuchis@gmail.com',
        description='A bunch of functions to analyze Nifti data in Python',
        long_description=open('README.md').read(),
        long_description_content_type='text/markdown',
        url='https://github.com/txusser/anapyze',
        packages=find_packages(exclude = ['pipelines', 'resources', 'matlab_scripts']),
        install_requires=requirements,
        classifiers=[
            'Programming Language :: Python :: 3',
            'License :: OSI Approved :: MIT License',
            'Operating System :: OS Independent',
            ],
        python_requires='>=3.6',
        )