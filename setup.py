import os
import glob
import shutil
from subprocess import Popen, PIPE
from setuptools import setup
from setuptools.command.install import install

class CustomInstallCommand(install):
    """
    Copy input/ into the package. This is necessary because `package_data` or
    `include_package_data` (with `MANIFEST.in`) expect the data to be within the
    posidonius python package and not in the root. On the other hand, data_files
    can copy files outside the posidonius python package but they are copied
    directly at the "venv/" and not within the installed posidonius package.
    """
    def run(self):
        install.run(self)
        # Copy the input directory to the desired location
        src_dir = 'input'
        dst_dir = os.path.join(self.install_lib, 'posidonius', 'input')
        shutil.copytree(src_dir, dst_dir)

def remove_building_dirs(building_dirs):
    for dirname in ["build/", "posidonius.egg-info/"]:
        print("Removing {}".format(dirname))
        try:
            shutil.rmtree(dirname)
        except OSError:
            pass

def remove_older_wheels(dist_dir="dist"):
    # Get the list of wheel files and source files in the dist directory
    for pattern in ("*.whl", "*.tar.gz",):
        wheel_files = glob.glob(os.path.join(dist_dir, pattern))

        # Sort the wheel files by modification time, newest first
        wheel_files.sort(key=os.path.getmtime, reverse=True)

        # Keep the latest wheel and remove the older ones
        for wheel_file in wheel_files[1:]:
            print(f"Removing older wheel: {wheel_file}")
            os.remove(wheel_file)

# Setup configuration
setup(
    cmdclass={'install': CustomInstallCommand}, # <= necessary to include "input/" given that it is outside the "posidonius/" package
)

# Cleanup (it will only work with `pip install .`, since `python -m build` uses an isolated environment)
remove_building_dirs(["build/", "posidonius.egg-info/"])
remove_older_wheels("dist/")
