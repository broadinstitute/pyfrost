import os
import sys
from pathlib import Path
import setuptools
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.sysconfig import customize_compiler

import versioneer


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __str__(self):
        import pybind11
        return pybind11.get_include()


def bifrost_sources():
    file_exts = ["*.c", "*.cpp"]
    bifrost_path = Path(__file__).parent / "vendor/bifrost/src"

    for ext in file_exts:
        for fname in bifrost_path.glob(ext):
            if fname.name in {"Bifrost.cpp", "xxhash.c"}:
                continue

            yield str(fname)


ext_modules = [
    Extension(
        'bifrost_python',
        # Sort input source files to ensure bit-for-bit reproducible builds
        # (https://github.com/pybind/python_example/pull/53)
        sorted([
            *bifrost_sources(),
            "src/bifrost_python/Kmer.cpp",
            "src/bifrost_python/KmerCounter.cpp",
            "src/bifrost_python/UnitigColors.cpp",
            "src/bifrost_python/UnitigDataProxy.cpp",
            "src/bifrost_python/AdjacencyProxy.cpp",
            "src/bifrost_python/NodeView.cpp",
            "src/bifrost_python/EdgeView.cpp",
            "src/bifrost_python/BifrostDiGraph.cpp",
            "src/bifrost_python/bifrost_python.cpp",
        ]),
        define_macros=[
            ('XXH_NAMESPACE', 'BIFROST_HASH_'),
            ('MAX_KMER_SIZE', os.environ.get("MAX_KMER_SIZE", "32")),
        ],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            'vendor/bifrost/src',
            'vendor/cereal/include/',
            'vendor/robin-hood-hashing/include/'
        ],
        language="c++"
    ),
]


# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    import os
    with tempfile.NamedTemporaryFile('w', suffix='.cpp', delete=False) as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        fname = f.name
    try:
        compiler.compile([fname], extra_postargs=[flagname])
    except setuptools.distutils.errors.CompileError:
        return False
    finally:
        try:
            os.remove(fname)
        except OSError:
            pass
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14/17] compiler flag.

    The newer version is prefered over c++11 (when it is available).
    """
    flags = ['-std=c++14']

    for flag in flags:
        if has_flag(compiler, flag):
            return flag

    raise RuntimeError('Unsupported compiler -- at least C++14 support '
                       'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': ['-Wno-unused-private-field', '-Wno-unused-function',
                 '-Wno-unused-variable', '-Wno-sign-compare'],
    }
    l_opts = {
        'msvc': [],
        'unix': ['-lz', '-lpthread'],
    }

    if sys.platform == 'darwin':
        darwin_opts = ['-stdlib=libc++', '-mmacosx-version-min=10.7']
        c_opts['unix'] += darwin_opts
        l_opts['unix'] += darwin_opts

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])
        if ct == 'unix':
            extra_cpp_flags = cpp_flag(self.compiler)

            os.environ['CPPFLAGS'] = os.environ.get('CPPFLAGS', "") + f" {extra_cpp_flags}"
            customize_compiler(self.compiler)

            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')

        for ext in self.extensions:
            if not ext.define_macros:
                ext.define_macros = []

            ext.define_macros.append(('VERSION_INFO', '"{}"'.format(self.distribution.get_version())))
            ext.extra_compile_args = opts
            ext.extra_link_args = link_opts

        super().build_extensions()

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

    ext_modules=ext_modules,
    cmdclass=dict(build_ext=BuildExt, **versioneer.get_cmdclass()),
    zip_safe=False,
)
