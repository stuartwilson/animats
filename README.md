# NeoArealize

A model of reaction-diffusion pattern formation in neocortex.

Pre-requisites:

Build and install morphologica from:

https://github.com/ABRG-Models/morphologica

Make sure these packages are installed (Debian/Ubuntu example):

sudo apt install python python-numpy xterm

Now build:

```bash
cd NeoArealize
mkdir build
cd build
cmake ..
make
cd ..
```

To run:

```bash
python sim2.py
```

or, more recently, as I've removed the need for the sim.py master program, just

```bash
./build/sim/process w0 logs/w0 1
```
where w0 is the sim name, logs/w0 will be the logged output and 1 is a seed
for the random number generator.
