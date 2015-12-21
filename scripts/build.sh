#!/bin/bash

cd build/gnu/ &&
cmake ../../qsim/ &&
make &&
cd ../../