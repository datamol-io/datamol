from setuptools import setup
from setuptools import find_packages

setup(
    name="datamol",
    version="0.2.7",
    author="Hadrien Mary",
    author_email="hadrien.mary@gmail.com",
    url="https://github.com/datamol-org/datamol",
    description="A python library to work with molecules. Built on top of RDKit.",
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "License :: OSI Approved :: Apache-2.0",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    entry_points={"console_scripts": []},
)
