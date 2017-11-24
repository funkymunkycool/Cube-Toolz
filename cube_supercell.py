#!/usr/bin/env python

import cube_tools
from sys import argv


rCube = cube_tools.cube(argv[1])
rCube.super_cube([2,2,2])
rCube.write_cube('test.cube')
