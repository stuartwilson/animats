# animats

Pseudo-physical models of the interaction between animats and their environments. 

Pre-requisites:

Build and install morphologica from:

https://github.com/ABRG-Models/morphologica

Make sure these packages are installed (Debian/Ubuntu example):

sudo apt install python python-numpy xterm

Now build:

```bash
cd animats
mkdir build
cd build
cmake ..
make
cd ..
```

To run:

Run using the following

```bash
./build/sim/process w0 logs/w0 1
```
where w0 is the sim name, logs/w0 will be the logged output and 1 is a seed
for the random number generator.
