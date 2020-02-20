#!/usr/bin/env python
# File       : post.py
# Created    : Thu Feb 20 2020 02:27:15 PM (+0100)
# Author     : Fabian Wermelinger
# Description: Post-processing benchmark data
# Copyright 2020 ETH Zurich. All Rights Reserved.
import numpy as np
import pandas as pd
import argparse

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--files', nargs='+', type=str, help='Data files',
        required=True)
    parser.add_argument('-o', '--output', type=str, help='Output file',
        required=True)
    return parser.parse_known_args()

def main():
    args, _ = parseArgs()
    cols = ["blocks", "t_init", "t_comp", "t_dump", "t_tot"]
    da = pd.DataFrame(columns = cols)
    for f in sorted(args.files):
        df = pd.read_csv(f, sep='\s+', names = cols)
        ds = df.describe()
        da = da.append(ds.loc['mean'])
        with open(f, 'a') as out:
            out.write(ds.to_string())
    with open(args.output, 'w') as out:
        out.write(da.to_string(index=False))

if __name__ == "__main__":
    main()
