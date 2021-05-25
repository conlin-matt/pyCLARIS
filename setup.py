import setuptools

setuptools.setup(
    name="pyCLARIS",
    version="0.1.2",
    url="https://github.com/conlin-matt/pyCLARIS",
    author="Matthew P. Conlin",
    author_email="conlinm@ufl.edu",
    description="Analysis of CLARIS data in Python- tools and projects.",
    packages=setuptools.find_packages(),
    install_requires=['alphashape','matplotlib','numpy','numpy_groupies','pptk','scipy'],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
)
