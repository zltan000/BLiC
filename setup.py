import setuptools

with open("README.md", "r") as f:
  long_description = f.read()

setuptools.setup(
  name="BLiC",
  version="2023.5.0",
  author="Zhenlin Tan, Jiaxin Han",
  author_email="zltan999@sjtu.edu.cn",
  description="A python package to Build Light-cones and output mock Catalogues",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/zltan000/BLiC",
  packages=setuptools.find_packages(),
  classifiers=[
  "Programming Language :: Python :: 3",
  "License :: OSI Approved :: MIT License",
  "Operating System :: OS Independent",
  ],
)
