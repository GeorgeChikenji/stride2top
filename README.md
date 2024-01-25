# STRIDE2TOP

STRIDE2TOP analyzes the topology of β-sheet-containing protein structures based on hydrogen-bond network analysis.

## Requirement
* STRIDE2TOP requires the output file of `stride -h` command or `stride` command in your `PATH`.


## Main Features
* Output the arrangement of β Strands in a Beta-Sheet.
* Output the arrangement in the [graphviz](http://www.graphviz.org) dot file style.
* Find and list the cycle paths in the graph.
* Find left-handed connections such as followings. (The first and the last β strands are Parallel.)
    * β - α - β
    * β - Loop - β
    * β - β (on other sheets)  - β


# Installation

## Build from source

To build from source, you need a compiler that supports C++14 features.
Go into the ./src directory and build with GNU make.
Specifying `-j N` option to `make` command is recommended, as it took some time to compile all the sources.
N is the number of jobs to run in parallel (and is usually the number of cores your CPU has).

```sh
cd src
make -j N
```

# Usage

Run with `--help` option to show the detailed help message.


# Licensing

- This software is distributed under the MIT License.
- This software contains a subset of [the Boost libraries][boost]. The Boost libraries are licensed as the Boost Software License 1.0.
- This software contains a copy of [Eigen][eigen] which is primarily MPL2 licensed and some files are
BSD or LGPL.
- This software contains a copy of [Google Test][googletest] which is BSD licensed.


[boost]: http://www.boost.org
[eigen]: http://eigen.tuxfamily.org
[googletest]: https://github.com/google/googletest
