# Author: Xiaotao Wang

"""
This is a free software under GPLv3. Therefore, you can modify, redistribute
or even mix it with other GPL-compatible codes. See the file LICENSE
included with the distribution for more details.

"""
import os, sys, raichu, glob
import setuptools

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if (sys.version_info.major!=3) or (sys.version_info.minor<7):
    print('PYTHON 3.7+ IS REQUIRED. YOU ARE CURRENTLY USING PYTHON {}'.format(sys.version.split()[0]))
    sys.exit(2)

# Guarantee Unix Format
for src in glob.glob('scripts/*'):
    text = open(src, 'r').read().replace('\r\n', '\n')
    open(src, 'w').write(text)

setuptools.setup(
    name = 'RaichuNorm',
    version = raichu.__version__,
    author = raichu.__author__,
    author_email = 'wangxiaotao@fudan.edu.cn',
    url = 'https://github.com/XiaoTaoWang/Raichu',
    description = 'A cross-platform method for chromatin contact normalization',
    keywords = 'Hi-C ChIA-PET HiChIP PLAC-Seq single-cell normalization',
    packages = setuptools.find_packages(),
    scripts = glob.glob('scripts/*'),
    long_description = read('README.rst'),
    long_description_content_type='text/x-rst',
    classifiers = [
        'Programming Language :: Python :: 3 :: Only',
        'Operating System :: POSIX',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ]
    )

