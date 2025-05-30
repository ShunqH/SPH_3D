# SPH 1D CPP

This is an SPH code for solving the 3D Sod Shock Tube problem.

This code is fully rewrited from SPH_CPP(1D) I uploaded before, with an improved data structrue, neighbor particles searching, and more. 

## Running the Simulation

1. **Compile the code:**
    ```bash
    make
    ```

2. **Navigate to the `bin` directory:**
    ```bash
    cd ./bin
    ```

3. **Run the simulation:**
    ```bash
    sod.sph -i ../input/sod.in
    ```

## Analyzing the Output

You can use the Python script located in `./pyscript` to read the output files.

## Modifying the Setup

To change the setup, edit the file located at:
```bash
./input/sod.in


