from setuptools import setup, find_packages
 
setup(name="ImPhSh",
      version="0.1",
      license="MIT",
      author="Gigi Sayfan",
      description="Add static script_dir() method to Path",
      packages=find_packages(exclude=["tests"]),
      long_description=open("README.md").read(),
      zip_safe=False)