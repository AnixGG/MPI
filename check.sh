#!/bin/bash

clang-format -i lab1/$1
clang-tidy lab1/$1 -- -std=c17 -I/usr/lib/x86_64-linux-gnu/openmpi/include