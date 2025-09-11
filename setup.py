from setuptools import setup, find_packages

setup(
    name="analysistools",
    version="0.1.0",
    packages=find_packages(),
    license="GPL-3.0",
    description="Various scripts to process astrophysical simulations",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Chris Power (chris.power@uwa.edu.au)",
    install_requires=["h5py",
                      "numpy",
		      "scipy",
		      "matplotlib",
		      "plotly",
                     ],  

)
