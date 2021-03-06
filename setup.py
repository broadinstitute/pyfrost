import os
import re
import sys
import platform
import subprocess
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

from typing import Union, List

import versioneer


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir='', target: Union[List[str], str]=None):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)
        self.target = target


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(),
                                                                           extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''), self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)

        if ext.target:
            if isinstance(ext.target, list):
                for target in ext.target:
                    subprocess.check_call(['cmake', '--build', '.', '--target', target] + build_args,
                                          cwd=self.build_temp)
            else:
                subprocess.check_call(['cmake', '--build', '.', '--target', ext.target] + build_args,
                                      cwd=self.build_temp)
        else:
            subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)


setup(
    name='pyfrost',
    version=versioneer.get_version(),
    author='Lucas van Dijk',
    author_email='lvandijk@broadinstitute.org',
    description='A Python interface to the Bifrost colored compacted De Bruijn graph '
                'library with a NetworkX like API',
    long_description='',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,

    setup_requires=['cmake>=3.10'],
    install_requires=['networkx>=2.4'],

    ext_modules=[CMakeExtension('pyfrostcpp', target=['pyfrostcpp'])],
    cmdclass=dict(build_ext=CMakeBuild, **versioneer.get_cmdclass()),
    zip_safe=False,
)
