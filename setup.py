from setuptools import setup, find_packages

setup(
    name="fibrilsite",
    version="1.0",
    description="A package for aiding fibril site definition",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Ahmed Sadek",
    author_email="ahmed.sadek@epfl.ch",
    python_requires="==3.6",
    license="MIT",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "numpy==1.19.2",
        "pandas==0.25.3",
        "seaborn==0.11.2",
        "matplotlib==3.2.2",
        "biopython==1.79",
        "scipy==1.5.4",
        "StrBioInfo==0.9a0.dev1",
        "open3d== 0.15.1"
    ],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
