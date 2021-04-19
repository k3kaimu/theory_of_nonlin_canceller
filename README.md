# Theoretical Analysis Program for ``Theoretical Analysis of In-Band Full-Duplex Radios with Parallel Hammerstein Self-Interference Cancellers''

This is the source code of my paper "Theoretical Analysis of In-Band Full-Duplex Radios with Parallel Hammerstein Self-Interference Canceller" published on IEEE Transactions on Wireless Communications.

The author has confirmed that this program can compile and run on Ubuntu (on Windows Subsystem for Linux (WSL) 2).
I haven't checked it in other environments, but it seems to work fine.

Therefore, the following procedure is for Ubuntu, but you can run this program on Windows and macOS as well if you build the environment.


## How to compile and run on Ubuntu

### Step 0: Install build-essentials

```sh
$ sudo apt install build-essential
```

### Step 1: Install and activate a compiler for D Programming Language

```sh
$ curl -fsS https://dlang.org/install.sh | bash -s dmd-2.091.1
$ source ~/dlang/dmd-2.091.1/activate
```

### Step 2: Build this program

```sh
$ cd directory/of/this/program
$ dub build --build=release
```


### Step 3: Run

```sh
$ ./analysis
```


### Step 4: Plot results of theoretical curves

To plot the results, you needs python3(and its libraries: matplotlib, numpy, and scipy) on your environment.
To install those tools, please see other instructions on the Internet.

After install python3 and the libraries, run following commands: 

```sh
$ cd plot
$ sh ./run.sh
```


## Citation

Kazuki Komatsu, Yuichi Miyaji, Hideyuki Uehara, ``Theoretical Analysis of In-Band Full-Duplex Radios with Parallel Hammerstein Self-Interference Cancellers,'' IEEE Transactions on Wireless Communications, (Early Access).

