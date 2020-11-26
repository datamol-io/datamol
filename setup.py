from setuptools import setup
from setuptools import find_packages

setup(
    name="datamol",
    version="0.2.5",
    author="Hadrien Mary",
    author_email="hadrien.mary@gmail.com",
    url="https://github.com/invivoai-platform/datamol",
    description="Python library to create and submit workflow on an argo server.",
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
    ],
    entry_points={
        "console_scripts": []
    },
)
