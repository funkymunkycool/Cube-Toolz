----
Cube Toolz
----

A python library and tool to read in and manipulate Gaussian cube files. This code allows you to:
- Read and write Gaussian cube files
- Translate and rotate cube data
- Integrate around a particular atom
- Integrate around a sphere
- Integrate around the whole cube file
- Take the planar average

Requires numpy and scipy. If you have pip, you can just run
```
pip install -r requirements.txt
```
to install all the requirements.

This is designed to be used either as a library, or to be run from the command line. If you decide to run it from the command line, here are the required arguments:

positional arguments:
  Files                 Cube files used in program

optional arguments: 
* -h, --help            show this help message and exit
* -a, --add             Add two or more cube files together
* -s, --subtract        Subtract two or more cube files together
* -p [POWER], --power [POWER]  Raise the cube file to a certain power. Any number of cube files can be specified and they will all be raised to the power defined. Default is to square the cube file(s).
* -t TRANSLATE TRANSLATE TRANSLATE, --translate TRANSLATE TRANSLATE TRANSLATE  Translate a cube file. Requires a translation vector as an argument.

usage: cube_tools.py [-h] [-a] [-s] [-p [POWER]] [-t TRANSLATE TRANSLATE TRANSLATE] Files [Files ...]
