# Triangle4CycleShuffle

This is a source code of the following paper:

Jacob Imola, Takao Murakami, Kamalika Chaudhuri, "Differentially Private Triangle and 4-Cycle Counting in the Shuffle Model," Proceedings of the 29th ACM Conference on Computer and Communications Security (CCS 2022), pp.1505-1518, 2022.

Full paper: https://arxiv.org/abs/2205.01429

# Directory Structure
- cpp/			&emsp;C/C++ codes.
- python/		&emsp;Python codes.
- LICENSE.txt		&emsp;MIT license.
- README.md		&emsp;This file.

# Usage

**(1) Install**

Install StatsLib (see cpp/README.md).

Install C/C++ codes as follows.
```
$ cd cpp/
$ make
$ cd ../
```

**(2) Download and preprocess Gplus**

Download the [Gplus dataset](https://snap.stanford.edu/data/ego-Gplus.html) and place the dataset in data/Gplus/.

Run the following commands.

```
$ cd python/
$ python3 Read_Gplus.py data/Gplus/gplus_combined.txt data/Gplus/edges.csv data/Gplus/deg.csv
$ cd ../
```

Then the edge file (edges.csv) and degree file (deg.csv) will be output in data/Gplus/.

**(3) Download and preprocess IMDB**

Download the [IMDB dataset](https://www.cise.ufl.edu/research/sparse/matrices/Pajek/IMDB.html) and place the dataset in data/IMDB/.

Run the following commands.

```
$ cd python/
$ python3 Read_IMDB.py data/IMDB/IMDB.mtx data/IMDB/edges.csv data/IMDB/deg.csv
$ cd ../
```

Then the edge file (edges.csv) and degree file (deg.csv) will be output in data/IMDB/.

**(4) Run subgraph counting algorithms**

Run the following commands ([Dataset] is "Gplus" or "IMDB").

```
$ cd cpp/
$ chmod +x run_SubgraphShuffle.sh
$ ./run_SubgraphShuffle.sh [Dataset (Gplus/IMDB)]
$ cd ../
```

Then experimental results for seven algorithms (WShuffle_triangle, WShuffle_triangle^* (w/o Lap), WShuffle_triangle^*, WLocal_triangle, ARR_triangle, WShuffle_square, WLocal_square) with n=20000 and (epsilon,delta)=(1,10^{-8}) will be output in data/[Dataset]/. For more details of parameter settings, see Usage of SubgraphShuffle.

# Execution Environment
We used CentOS 7.5 with gcc 4.8.5 and python 3.6.5.

# External Libraries used in our source code.
- [StatsLib](https://www.kthohr.com/statslib.html) is distributed under the [Apache License 2.0](https://github.com/kthohr/stats/blob/master/LICENSE).
