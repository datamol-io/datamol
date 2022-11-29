from setuptools import setup
from setuptools import find_packages

# Sync the env.yml file here
install_requires = [
    "tqdm",
    "loguru",
    "joblib",
    "fsspec>=2021.9",
    "pandas",
    "numpy",
    "scipy",
    "matplotlib",
    "pillow",
    "selfies",
    "appdirs",
    "scikit-learn",
    "packaging",
    # NOTE(hadim): can't add rdkit because of `pip` will always override
    # the conda package at the moment.
    # See:
    # - https://github.com/rdkit/rdkit/issues/5378
    # - https://github.com/conda-forge/rdkit-feedstock/issues/104
    # "rdkit",
]

setup(
    name="datamol",
    version="0.8.7",
    author="Valence Discovery",
    author_email="hadrien@valencediscovery.com",
    url="https://github.com/datamol-org/datamol",
    description="A python library to work with molecules. Built on top of RDKit.",
    long_description=open("README.md", encoding="utf8").read(),
    long_description_content_type="text/markdown",
    project_urls={
        "Bug Tracker": "https://github.com/datamol-org/datamol/issues",
        "Documentation": "https://doc.datamol.io",
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
