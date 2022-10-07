#!/usr/bin/env python

from .src import HPCBlast, HPCBlastArg


def main():
    HPCBlast(*HPCBlastArg()).run()


if __name__ == "__main__":
    main()
