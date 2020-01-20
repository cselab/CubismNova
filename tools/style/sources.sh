#!/usr/bin/env bash
# File       : sources.sh
# Created    : Mon Jan 20 2020 11:29:28 AM (+0100)
# Author     : Fabian Wermelinger
# Description: Source file regex patterns
# Copyright 2020 ETH Zurich. All Rights Reserved.

HEADER='\.\./\.\./include/Cubism/.*\.h'
SOURCE='\.\./\.\./src/.*\.(c|cpp)'
TEST_H='\.\./\.\./test/[A-Z].*\.h'
TEST_S='\.\./\.\./test/[A-Z].*\.(c|cpp)'
