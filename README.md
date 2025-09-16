# SPH 3D CPP

This is an SPH code for solving the 3D problem.

The code is fully rewritten from SPH_CPP (1D) that I uploaded before, with an improved data structure, neighbor particle searching, and more. 

## Running the Simulation

1. **Configure your setup:**
    ```bash
    python ./configure.py --setup=sod
    ```

    To enable OpenMP, add the `--openmp` flag:
    ```bash
    python ./configure.py --setup=sod --openmp
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
    ./xeno.sph -i ../input/sod.in
    ```

## Analyzing the Output

You can use the Python scripts located in `./pyscript` to read and analyze the output files.
    ```bash
    cd pyscript
    python data.py
    ```

## Modifying the Input

To change the input parameters, edit:
```bash
./input/sod.in
```

## Modifying the Setup

To change the setup, edit one of the files in:
```bash
./src/setup/setup_xxx.cpp
```

You can also create your own setup file by following the examples. (Check `./src/setup/setup_box.cpp`. )

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