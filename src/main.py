#!/usr/bin/env python

from .src import *


def main():
    args, blast_options = HPCBlastArg()
    app = HPCBlast(args, blast_options)
    app.run()


if __name__ == "__main__":
    main()
