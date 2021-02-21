from setuptools import setup
from setuptools import find_packages

# Sync the env.yml file here
install_requires = [
    "tqdm",
    "loguru",
    "joblib",
    "fsspec>=0.8.0",
    "pandas",
    "numpy",
    "scipy",
    "pillow",
    "ete3",
    "selfies",
]

setup(
    name="datamol",
    version="0.2.7",
    author="Hadrien Mary",
    author_email="hadrien.mary@gmail.com",
    url="https://github.com/datamol-org/datamol",
    description="A python library to work with molecules. Built on top of RDKit.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    project_urls={
        "Bug Tracker": "https://github.com/datamol-org/datamol/issues",
        "Documentation": "https://datamol.readthedocs.io/en/stable/",
        "Source Code": "https://github.com/datamol-org/datamol",
    },
    python_requires=">=3.7",
    install_requires=install_requires,
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
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
