import os
from setuptools import setup, find_packages

pkg = 'rfinder'
build_root = os.path.dirname(__file__)

scripts = ["bin/" + i for i in os.listdir("bin")]


def readme():
    """Get readme content for package long description"""
    with open(os.path.join(build_root, 'README.md')) as f:
        return f.read()


def requirements():
    """Get package requirements"""
    with open(os.path.join(build_root, 'requirements.txt')) as f:
        return [pname.strip() for pname in f.readlines()]


setup(name=pkg,
      version="1.0.0",
      description="A tool to identify the presence of RFI in an observation and visualize it",
      long_description=readme(),
      author="Filippo Marcello Maccagni",
      author_email="filippo.maccagni@gmail.com",
      packages=find_packages(),
      url="https://github.com/Fil8/RFInder",
      license="GNU GPL 3",
      classifiers=["Intended Audience :: Developers",
                   "Programming Language :: Python :: 2",
                   "Topic :: Scientific/Engineering :: Astronomy",
                   "Topic :: Software Development :: Libraries :: Python Modules"],
      platforms=["OS Independent"],
      install_requires=requirements(),
      scripts=scripts)
