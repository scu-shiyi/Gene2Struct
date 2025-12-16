from setuptools import setup, find_packages
import os


setup(
    name="phyloselect",
    version="1.1.0",
    description="PhyloSelect: A comprehensive toolkit for phylogenetic analysis and evolutionary selection.",
    author="Shishi",
    author_email="shi@stu.scu.edu.cn",
    url="https://github.com/Shishi/gene2struct",  
    packages=find_packages(include=["scripts", "scripts.*"]),
    python_requires='>=3.8',
    install_requires=[
        "biopython",
        "matplotlib",
        "numpy",
        "pandas",
        "scipy",
        "seaborn",
        "psutil",
        "ete3",
        "requests"],
    entry_points={
        'console_scripts': [
            'phyloselect = phyloselect.cli:main',
        ]
    },
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
)
