# How to build
```
mkdir build && cd build && cmake .. && make -j && cp ../in_5K.fluid . && cp ../in_15K.fluid .
```

# How to run
```
./fluidanimate 4 1 in_5K.fluid output_file
```

or 

```
./fluidanimate 4 1 in_15K.fluid output_file
```

# Gist
```Cell```s are what is shared among threads. We kept track addresses accessed by threads. Then we found which thread pairs accessed how many common addresses. Pinning will be done with the help of this.

# Questions

#define CACHELINE_SIZE 128  // aydin: how is this 128?
