#!/usr/bin/env python3
# File       : gen_version.py
# Created    : Mon Mar 15 2021 12:20:47 PM (+0100)
# Author     : Fabian Wermelinger
# Description: Generate CubismNova build version strings
# Copyright 2021 ETH Zurich. All Rights Reserved.
import os
import argparse


def parse_args(*, partial=False):
    parser = argparse.ArgumentParser()
    # yapf: disable
    parser.add_argument('-i', '--input', type=str, help="Input file template", required=True)
    parser.add_argument('-o', '--output', type=str, help="Output source file", required=True)
    parser.add_argument('-p', '--project_root', default=".", type=str, help="Root directory of project")
    # yapf: enable
    if partial:
        return parser.parse_known_args()
    else:
        return parser.parse_args()


def main(args):
    fversion = os.path.join(args.project_root, 'VERSION')
    with open(fversion, 'r') as f:
        version = f.read().strip()
    major, minor, patch = version.split('.')

    HEAD = ''
    SHA1 = ''
    fHEAD = os.path.join(args.project_root, '.git', 'HEAD')
    if os.path.isfile(fHEAD):
        with open(fHEAD, 'r') as git_head:
            HEAD = git_head.read().strip()
            HEAD = HEAD.split()
            if len(HEAD) > 1:
                HEAD = HEAD[1]
                fSHA1 = os.path.join(args.project_root, '.git', HEAD)
                with open(fSHA1, 'r') as git_sha1:
                    SHA1 = git_sha1.read().strip()
            else:
                SHA1 = HEAD

    with open(args.input, 'r') as f:
        version_template = f.read()

    version_template.replace('@VERSION@', version)
    version_template.replace('@HEAD@', HEAD)
    version_template.replace('@SHA1@', SHA1)
    with open(args.output, 'r') as out:
        out.write(version_template)


if __name__ == "__main__":
    args = parse_args()
    main(args)
