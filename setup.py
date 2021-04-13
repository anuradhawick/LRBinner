#!/usr/bin/env python3

import setuptools
import subprocess
from distutils.command.build import build
from setuptools.command.install import install as SetuptoolsInstall
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

class SHBuild(build):
    def run(self):
        try:
            subprocess.check_call(["sh", "build.sh"])
        except subprocess.CalledProcessError as e:
            sys.exit("Compilation error: ", e)    

packages = setuptools.find_packages()
package_data = {
    "mbcclr_utils": [
        "mbcclr_utils/*",
        "mbcclr_utils/bin/*",
    ]
}
data_files = [(".", ["LICENSE", "README.md"])]

setuptools.setup(
    name="MetaBCC-LR-2", # Replace with your own username
    version="0.0.1",
    zip_safe=True,
    author="Anuradha Wickramarachchi",
    author_email="anuradhawick@gmail.com",
    description="MetaBCC-LR-2",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/anuradhawick/MetaBCC-LR",
    packages=packages,
    package_data=package_data,
    data_files=data_files,
    include_package_data=True,
    scripts=['MetaBCC-LR'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "biopython",
        "tqdm",
        "tabulate",
        "seaborn",
        "setuptools",
        "pytorch"],
    python_requires='>=3.7',
    cmdclass={  'build': SHBuild }
)