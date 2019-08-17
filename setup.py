import os
import glob
import shutil
from subprocess import Popen, PIPE

try:
    from setuptools import setup
except:
    from distutils.core import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except (IOError, ImportError):
    long_description = read('README.md')

with open('requirements.txt') as f:
    required = f.read().splitlines()

def get_git_version(default="v0.0.1"):
    try:
        p = Popen(['git', 'describe', '--tags'], stdout=PIPE, stderr=PIPE)
        p.stderr.close()
        line = p.stdout.readlines()[0]
        line = line.strip()
        return line.decode('utf-8')
    except:
        return default

setup(
    name = "posidonius",
    version = get_git_version(default="v0.0.1"),
    author = "Sergi Blanco-Cuaresma",
    author_email = "marblestation@users.noreply.github.com",
    description = ("Case generator for Posidonius."),
    license = "GNU Affero General Public License",
    keywords = "N-Body simulations exoplanets tides",
    url = "http://www.blancocuaresma.com/s/",
    packages=['posidonius', 'posidonius.analysis', 'posidonius.effects', 'posidonius.particles', 'posidonius.integrator' ],
    data_files=[(basedir, [filename for filename in glob.iglob('{}/*.*'.format(basedir))]) for basedir in glob.iglob("input/*")],
    long_description=long_description,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Rust",
        "Programming Language :: Python",
        "Intended Audience :: Science/Research",
        "Environment :: Console",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: GNU Affero General Public License v3",
    ],
    install_requires=required,
)

for dirname in ["build/", "dist/", "posidonius.egg-info/"]:
    print("Removing {}".format(dirname))
    try:
        shutil.rmtree(dirname)
    except OSError:
        pass
