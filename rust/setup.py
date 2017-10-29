import os
import glob
import shutil
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
    packages=['posidonius', 'posidonius.analysis' ],
    data_files=[(basedir, [filename for filename in glob.iglob('{}/*.*'.format(basedir))]) for basedir in glob.iglob("input/*")],
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
        'pandas',
        'matplotlib',
    ],
)

for dirname in ["build/", "dist/", "posidonius.egg-info/"]:
    print("Removing {}".format(dirname))
    try:
        shutil.rmtree(dirname)
    except OSError:
        pass
