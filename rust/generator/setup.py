import os
import glob
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "posidonius",
    version = "0.0.1",
    author = "Sergi Blanco-Cuaresma",
    author_email = "marblestation@users.noreply.github.com",
    description = ("Case generator for Posidonius."),
    license = "GNU Affero General Public License",
    keywords = "N-Body simulations exoplanets tides",
    url = "http://www.blancocuaresma.com/s/",
    packages=['posidonius', ],
    data_files=[(basedir[3:], [filename for filename in glob.iglob('{}/*.*'.format(basedir))]) for basedir in glob.iglob("../input/*")],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Rust",
        "Programming Language :: Python",
        "Intended Audience :: Science/Research",
        "Environment :: Console",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: GNU Affero General Public License v3",
    ],
    install_requires=[
        'numpy',
    ],
)
