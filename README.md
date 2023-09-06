# How to build
```
mkdir build && cd build && cmake .. && make -j && cp ../in_5K.fluid .
```

# How to run
```
./fluidanimate 4 1 in_5K.fluid output_file
```

# How to optimize cache coherence
gi: 0, sx: 0, ex: 8, sy: 0, ey: 21, sz: 0, ez: 8
gi: 1, sx: 0, ex: 8, sy: 0, ey: 21, sz: 8, ez: 15
gi: 2, sx: 8, ex: 15, sy: 0, ey: 21, sz: 0, ez: 8
gi: 3, sx: 8, ex: 15, sy: 0, ey: 21, sz: 8, ez: 15

# Tip
```Cell```s are what is shared among threads.

# Questions

#define CACHELINE_SIZE 128  // aydin: how is this 128?
