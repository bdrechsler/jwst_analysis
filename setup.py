from setuptools import setup, find_packages

setup(
    name="jwst_analysis",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        # List your dependencies here
        # e.g., 'requests>=2.25.1'
    ],
    author="Blake Drechsler",
    author_email="wbd814@gmail.com",
    description="Package to help analyze jwst data",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/bdrechsler/jwst_analysis.git",  # Optional
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # Choose appropriate license
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)