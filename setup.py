from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="genome3danalysis",
    version="0.3",
    author="Francesco Musella, Ye West",
    author_email="fmusella@g.ucla.edu",
    description="3D Genome Analysis Tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/alberlab/genome3danalysis",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
    ],
    python_requires='>=3.7.3',
    install_requires=[
        "numpy>=1.20.3",
        "alabtools>=1.1.13",
        "scipy>=1.7.0",
        "markov_clustering>=0.0.6.dev0",
        "mrcfile>=1.4.3",
    ],
    entry_points={
        'console_scripts': [
              'structfeat-run=genome3danalysis.structfeat.run_feature_extractor:main',
        ],
    }
)
