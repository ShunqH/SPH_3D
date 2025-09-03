# SPH 3D CPP

This is an SPH code for solving the 3D problem.

The code is fully rewritten from SPH_CPP (1D) that I uploaded before, with an improved data structure, neighbor particle searching, and more. 

## Running the Simulation

1. **Configure your setup:**
    ```bash
    python ./configure.py --setup=3Dsod
    ```

2. **Compile the code:**
    ```bash
    make
    ```

3. **Navigate to the `bin` directory:**
    ```bash
    cd ./bin
    ```

4. **Run the simulation:**
    ```bash
    ./xeno.sph -i ../input/3Dsod.in
    ```

## Analyzing the Output

You can use the Python scripts located in `./pyscript` to read and analyze the output files.

## Modifying the Input

To change the input parameters, edit:
```bash
./input/3Dsod.in
```

## Modifying the Setup

To change the setup, edit one of the files in:
```bash
./src/setup_xxx.cpp
```

You can also create your own setup file by following the examples.

Then modify the dictionary in
```bash
./configure.py
```
and rerun:
```bash
python ./configure.py --setup=xxx
```
Or manually edit:
```bash
config.mk
```

Finally, recompile the code:
```bash
make
```