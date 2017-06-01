
# Build and install

```bash
cargo build --release
cargo install
```

# Start default simulation

The default simulation is hard-coded in Posidonius, it can be changed by modifying the source code and re-building it again.

```bash
posidonius default target/default.bin target/default_history.bin
```

# Prepare initial snapshots

Using the generator, the user can design his/her own simulation with a python script. The generator will create a json file with the simulation description, which can be read later on by Posidonius to start it.

```bash
python generator/case3.py target/case3.json
python generator/case4.py target/case4.json
python generator/case7.py target/case7.json
python generator/example.py target/example.json
```

# Start simulation from initial snapshot

A simulation can be started using its json file (describing the simulation). The recovery and historic snapshot file names should be specified. The former will contain the information needed to resume interrupted simulations, while he later stores the evolution of the simulation over the years.

```bash
posidonius start target/case3.json target/case3.bin target/case3_history.bin
posidonius start target/case4.json target/case4.bin target/case4_history.bin
posidonius start target/case7.json target/case7.bin target/case7_history.bin
posidonius start target/example.json target/example.bin target/example_history.bin
```

# Resume interrupted simulations

Interrupted simulations can be restored using the recovery snapshot file. The historic snapshot filename has to be specified also to continue storing the history of the simulation.

```bash
posidonius resume target/case3.bin target/case3_history.bin
posidonius resume target/case4.bin target/case4_history.bin
posidonius resume target/case7.bin target/case7_history.bin
posidonius resume target/example.bin target/example_history.bin
```

# Analyse simulations

While a simulation is in progress or when it has ended, the historic snapshot file can be interpreted to generate a plain text file and a plot with the history of the simulation:

```bash
python analysis/process.py target/case3_history.bin
python analysis/process.py target/case4_history.bin
python analysis/process.py target/case7_history.bin
python analysis/process.py target/example_history.bin
```

To explore what possible resonances might be present in the system:

```
python analysis/process_timed_resonances.py target/case3_history.bin
python analysis/process_timed_resonances.py target/case4_history.bin
python analysis/process_timed_resonances.py target/case7_history.bin
python analysis/process_timed_resonances.py target/example_history.bin
```

Finally, to study a given resonance (e.g., 3:2) between planet one and two:

```
python analysis/process_single_resonance.py target/case3_history.bin 1 2 3 2
python analysis/process_single_resonance.py target/case4_history.bin 1 2 3 2
python analysis/process_single_resonance.py target/case7_history.bin 1 2 3 2
python analysis/process_single_resonance.py target/example_history.bin 1 2 3 2
```

