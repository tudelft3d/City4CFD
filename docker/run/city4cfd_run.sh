#!/bin/bash

docker run --rm -v `pwd`:/data tudelft3d/city4cfd:latest $@
