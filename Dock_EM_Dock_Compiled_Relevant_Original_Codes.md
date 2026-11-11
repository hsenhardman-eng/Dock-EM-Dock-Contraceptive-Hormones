# Compilation of Relevant Original Dock–EM–Dock Codes

This compilation combines the previously extracted original code excerpts from the five uploaded Colab-exported Python files. Code has not been rewritten. Excerpts retain their source line ranges. They are **not** a single executable script; original notebooks include alternatives, debugging, and dependencies on external molecular files. Inclusion does not verify that every branch was executed.

Ligands: ethinyl estradiol (EE), levonorgestrel (LNG), norelgestromin (NGMN), medroxyprogesterone acetate (MPA). Each receptor has 12 possible ordered pairs of distinct ligands. The excerpts cover code relevant to the workflow, but should not be taken as proof of complete source coverage for all 60 pairs.


---

# Progesterone receptor (PR)

# PR — original Dock–EM–Dock code excerpts
**Source:** `pr_energyminimization.py` (26,652 lines).
**Provenance:** All code excerpts below are copied verbatim from the uploaded Colab-exported Python source, with original line ranges. Sections are excerpts, not a standalone runnable script. Original complete file is included in the ZIP.
**Scope:** The original notebook may contain debugging, corrections, alternate parameterizations, or non-executable Colab shell commands. These excerpts do not imply all code paths were executed successfully.

## Original code excerpts
### Source lines 653–741
```python

for package in [
    "openmm",
    "openff-toolkit",
    "openmmforcefields",
    "openff-interchange"
]:
    try:
        print(package, importlib.metadata.version(package))
    except importlib.metadata.PackageNotFoundError:
        print(package, "NOT INSTALLED")

print("\n========== PROTEIN FORCE FIELD ==========")

protein = app.PDBFile(str(protein_path))

forcefield = app.ForceField(
    "amber14-all.xml"
)

print("Protein atoms:", protein.topology.getNumAtoms())
print("Protein residues:", protein.topology.getNumResidues())

try:
    system = forcefield.createSystem(
        protein.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None
    )

    print("AMBER14 parameterization: SUCCESS")
    print("System particles:", system.getNumParticles())

except Exception as error:
    print("AMBER14 parameterization: FAILED")
    print("Error:", str(error)[:2000])

print("\n========== EE TOOLKIT CHECK ==========")

try:
    from openff.toolkit import Molecule, ForceField

    ee = Molecule.from_file(
        str(ligand_path),
        file_format="SDF",
        allow_undefined_stereo=False
    )

    print("OpenFF ligand loading: SUCCESS")
    print("EE atoms:", ee.n_atoms)
    print("EE bonds:", ee.n_bonds)
    print("EE charge:", ee.total_charge)

except Exception as error:
    print("OpenFF ligand loading: FAILED")
    print("Error:", str(error)[:2000])

print("\nNo energy minimization performed.")

import sys
import subprocess
import importlib.metadata

print("========== INSTALLING FORCE-FIELD DEPENDENCIES ==========")

command = [
    sys.executable, "-m", "pip", "install",
    "-q",
    "openff-toolkit",
    "openmmforcefields",
    "openff-interchange"
]

result = subprocess.run(
    command,
    capture_output=True,
    text=True
)

print("Installation return code:", result.returncode)

if result.returncode != 0:
    print("\nInstallation error:")
    print(result.stderr[-3000:])

print("\n========== PACKAGE VERIFICATION ==========")

for package in [
    "openff-toolkit",
```

### Source lines 1227–1322
```python

# Convert the ligand to an OpenMM topology and positions.
ee_topology = ee.to_topology().to_openmm()
ee_positions = ee.conformers[0].m_as("nanometer")

print("PR atoms:", protein.topology.getNumAtoms())
print("EE atoms:", ee.n_atoms)

# Combine receptor and ligand topologies.
modeller = app.Modeller(
    protein.topology,
    protein.positions
)

modeller.add(
    ee_topology,
    ee_positions * unit.nanometer
)

print("Combined atoms:", modeller.topology.getNumAtoms())

# AMBER14 protein + OpenFF Sage 2.2.1 ligand.
forcefield = app.ForceField("amber14-all.xml")

generator = SMIRNOFFTemplateGenerator(
    molecules=[ee],
    forcefield="openff-2.2.1"
)

forcefield.registerTemplateGenerator(generator.generator)

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=None
)

print("\n========== SYSTEM VALIDATION ==========")
print("System particles:", system.getNumParticles())
print("System forces:", system.getNumForces())

assert system.getNumParticles() == modeller.topology.getNumAtoms()

# Save the combined system and starting coordinates.
system_path = work / "PR_EE_combined_system.xml"
structure_path = work / "PR_EE_preminimization.pdb"

system_path.write_text(XmlSerializer.serialize(system))

with open(structure_path, "w") as f:
    app.PDBFile.writeFile(
        modeller.topology,
        modeller.positions,
        f,
        keepIds=True
    )

# Evaluate initial potential energy.
integrator = VerletIntegrator(0.001 * unit.picoseconds)
platform = Platform.getPlatformByName("Reference")

context = Context(system, integrator, platform)
context.setPositions(modeller.positions)

state = context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(
    unit.kilojoule_per_mole
)

print("Initial potential energy (kJ/mol):", round(energy, 3))
print("Energy finite:", bool(np.isfinite(energy)))

assert np.isfinite(energy), "Initial potential energy is not finite."

print("\nSaved system:", system_path)
print("Saved structure:", structure_path)
print("SUCCESS: Combined PR–EE system assembled.")
print("No energy minimization performed.")

del context, integrator
''')

env = os.environ.copy()
env["MPLBACKEND"] = "Agg"
env["PATH"] = str(env_path / "bin") + os.pathsep + env.get("PATH", "")

result = subprocess.run(
    [str(env_path / "bin/python"), str(script_path)],
    capture_output=True,
    text=True,
    env=env,
    timeout=600
)

print(result.stdout)

```

### Source lines 1336–1477
```python

script_path.write_text(r'''
from pathlib import Path
import numpy as np

from openmm import (
    XmlSerializer, CustomExternalForce,
    LangevinMiddleIntegrator, LocalEnergyMinimizer,
    Platform, Context, unit
)
from openmm import app

work = Path("/content/PR_Dock_EM_Dock")

system = XmlSerializer.deserialize(
    (work / "PR_EE_combined_system.xml").read_text()
)

structure = app.PDBFile(
    str(work / "PR_EE_preminimization.pdb")
)

positions = structure.positions

print("========== PR–EE ENERGY MINIMIZATION ==========")
print("Total atoms:", system.getNumParticles())

# Harmonic restraints on PR protein heavy atoms.
restraint = CustomExternalForce(
    "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
)

restraint.addGlobalParameter(
    "k",
    1000.0 * unit.kilojoule_per_mole / unit.nanometer**2
)

for parameter in ["x0", "y0", "z0"]:
    restraint.addPerParticleParameter(parameter)

protein_heavy_atoms = 0

for atom in structure.topology.atoms():
    if atom.residue.name != "UNL" and atom.element.symbol != "H":
        xyz = positions[atom.index].value_in_unit(unit.nanometer)
        restraint.addParticle(atom.index, xyz)
        protein_heavy_atoms += 1

assert protein_heavy_atoms > 0, "No protein heavy atoms restrained."

system.addForce(restraint)

print("Restrained protein heavy atoms:", protein_heavy_atoms)
print("Restraint strength: 1000 kJ/mol/nm²")
print("Maximum iterations: 1000")
print("Nonbonded treatment: NoCutoff")
print("Solvation: No explicit solvent")

integrator = LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1 / unit.picosecond,
    0.002 * unit.picoseconds
)

platform = Platform.getPlatformByName("CPU")
context = Context(system, integrator, platform)
context.setPositions(positions)

initial_energy = context.getState(
    getEnergy=True
).getPotentialEnergy().value_in_unit(
    unit.kilojoule_per_mole
)

print("\nInitial restrained energy (kJ/mol):",
      round(initial_energy, 3))

LocalEnergyMinimizer.minimize(
    context,
    tolerance=10.0,
    maxIterations=1000
)

state = context.getState(
    getEnergy=True,
    getPositions=True
)

final_energy = state.getPotentialEnergy().value_in_unit(
    unit.kilojoule_per_mole
)

final_positions = state.getPositions()

print("Final restrained energy (kJ/mol):",
      round(final_energy, 3))

print("Energy change (kJ/mol):",
      round(final_energy - initial_energy, 3))

assert np.isfinite(final_energy), "Final energy is not finite."

# Save minimized structure.
output_path = work / "PR_EE_minimized_complex.pdb"

with open(output_path, "w") as f:
    app.PDBFile.writeFile(
        structure.topology,
        final_positions,
        f,
        keepIds=True
    )

print("\nSaved:", output_path)
print("SUCCESS: PR–EE minimization completed.")
print("Convergence and structural validation still required.")

del context, integrator
''')

env = os.environ.copy()
env["MPLBACKEND"] = "Agg"
env["PATH"] = str(env_path / "bin") + os.pathsep + env.get("PATH", "")

result = subprocess.run(
    [str(env_path / "bin/python"), str(script_path)],
    capture_output=True,
    text=True,
    env=env,
    timeout=900
)

print(result.stdout)

if result.returncode != 0:
    print("========== ERROR DETAILS ==========")
    print(result.stderr[-6000:])

print("Process return code:", result.returncode)

from pathlib import Path
import numpy as np
```

### Source lines 1772–1860
```python
print("Restrained PR heavy atoms:", len(protein_heavy))
print("Restrained EE atoms:", 0)
print("Restraint strength: 1000 kJ/mol/nm²")
print("Maximum iterations: 1000")
print("Solvent model: None (NoCutoff)")

integrator = VerletIntegrator(
    0.001 * unit.picoseconds
)

platform = Platform.getPlatformByName("CPU")

context = Context(system, integrator, platform)
context.setPositions(positions)

initial_energy = context.getState(
    getEnergy=True
).getPotentialEnergy().value_in_unit(
    unit.kilojoule_per_mole
)

print("\nInitial potential energy (kJ/mol):",
      round(initial_energy, 3))

LocalEnergyMinimizer.minimize(
    context,
    tolerance=10.0,
    maxIterations=1000
)

state = context.getState(
    getEnergy=True,
    getPositions=True,
    getForces=True
)

final_energy = state.getPotentialEnergy().value_in_unit(
    unit.kilojoule_per_mole
)

final_positions = state.getPositions()

forces = state.getForces(asNumpy=True).value_in_unit(
    unit.kilojoule_per_mole / unit.nanometer
)

force_rms = float(np.sqrt(np.mean(np.sum(forces**2, axis=1))))

print("Final potential energy (kJ/mol):",
      round(final_energy, 3))

print("Energy change (kJ/mol):",
      round(final_energy - initial_energy, 3))

print("Final RMS force (kJ/mol/nm):",
      round(force_rms, 4))

assert np.isfinite(final_energy)
assert np.isfinite(force_rms)

output_path = work / "PR_EE_minimized_corrected.pdb"

with open(output_path, "w") as f:
    app.PDBFile.writeFile(
        structure.topology,
        final_positions,
        f,
        keepIds=True
    )

print("\nSaved:", output_path)
print("SUCCESS: Corrected minimization completed.")
print("Convergence and structural validation still required.")

del context, integrator
''')

env = os.environ.copy()
env["MPLBACKEND"] = "Agg"
env["PATH"] = str(env_path / "bin") + os.pathsep + env.get("PATH", "")

result = subprocess.run(
    [str(env_path / "bin/python"), str(script_path)],
    capture_output=True,
    text=True,
    env=env,
    timeout=900
)

```

### Source lines 2896–2984
```python
    str(ee_sdf),
    removeHs=False
)[0]

assert ee is not None

xyz = np.asarray(
    ee.GetConformer().GetPositions()
)

heavy_indices = [
    atom.GetIdx()
    for atom in ee.GetAtoms()
    if atom.GetSymbol() != "H"
]

center = xyz[heavy_indices].mean(axis=0).tolist()
box_size = [30.0, 30.0, 30.0]

print("========== PR–EE VINA VALIDATION ==========")
print("Receptor:", receptor_path.name)
print("Box center (Å):", [round(x, 3) for x in center])
print("Box size (Å):", box_size)

v = Vina(
    sf_name="vina",
    cpu=2,
    seed=2026,
    verbosity=1
)

v.set_receptor(str(receptor_path))

print("\nReceptor loaded successfully.")

v.compute_vina_maps(
    center=center,
    box_size=box_size
)

print("\nSUCCESS: Vina affinity maps computed.")
print("The PR–EE receptor is readable by Vina.")
print("No second ligand docking performed.")
''')

env = os.environ.copy()
env["PATH"] = str(env_path / "bin") + os.pathsep + env.get("PATH", "")
env["MPLBACKEND"] = "Agg"
env["OPENMM_CPU_THREADS"] = "2"

result = subprocess.run(
    [str(env_path / "bin/python"), str(script_path)],
    capture_output=True,
    text=True,
    env=env,
    timeout=600
)

print(result.stdout)

if result.returncode != 0:
    print("========== ERROR DETAILS ==========")
    print(result.stderr[-6000:])

print("Process return code:", result.returncode)

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

print("========== LNG INPUT FILE SEARCH ==========")

folders = [
    Path("/content"),
    Path("/content/PR_Dock_EM_Dock")
]

candidates = []

for folder in folders:
    if not folder.exists():
        continue

    for path in folder.iterdir():
        if path.is_file() and path.suffix.lower() in (
            ".sdf", ".mol", ".mol2"
        ):
            if "lng" in path.name.lower() or "levonorgestrel" in path.name.lower():
                candidates.append(path)
```

### Source lines 3179–3276
```python
ligand_path = work / "PR_EE_to_LNG_ligand.pdbqt"
ee_sdf = work / "PR_EE_minimized_ligand.sdf"

output_path = work / "PR_EE_to_LNG_docked.pdbqt"

for path in (receptor_path, ligand_path, ee_sdf):
    assert path.exists(), f"Missing required file: {path}"

ee = Chem.SDMolSupplier(
    str(ee_sdf), removeHs=False
)[0]

assert ee is not None

xyz = np.asarray(ee.GetConformer().GetPositions())

heavy_indices = [
    atom.GetIdx()
    for atom in ee.GetAtoms()
    if atom.GetSymbol() != "H"
]

center = xyz[heavy_indices].mean(axis=0).tolist()

print("========== PR: EE → LNG SEQUENTIAL DOCKING ==========")
print("Receptor:", receptor_path.name)
print("Incoming ligand:", ligand_path.name)
print("Box center (Å):", [round(x, 3) for x in center])
print("Box size (Å): [30.0, 30.0, 30.0]")
print("Exhaustiveness: 16")
print("Requested poses: 10")
print("Random seed: 2026")

v = Vina(
    sf_name="vina",
    cpu=2,
    seed=2026,
    verbosity=1
)

v.set_receptor(str(receptor_path))
v.set_ligand_from_file(str(ligand_path))

v.compute_vina_maps(
    center=center,
    box_size=[30.0, 30.0, 30.0]
)

print("\n========== RUNNING DOCKING ==========")

v.dock(
    exhaustiveness=16,
    n_poses=10
)

energies = v.energies(n_poses=10)

print("\n========== DOCKING RESULTS ==========")

for i, row in enumerate(energies, start=1):
    print(
        f"Pose {i}: "
        f"Affinity = {row[0]:.3f} kcal/mol"
    )

v.write_poses(
    str(output_path),
    n_poses=10,
    overwrite=True
)

print("\nSaved:", output_path)
print("SUCCESS: PR EE → LNG docking completed.")
''')

env = os.environ.copy()
env["PATH"] = str(env_path / "bin") + os.pathsep + env.get("PATH", "")
env["MPLBACKEND"] = "Agg"

result = subprocess.run(
    [str(env_path / "bin/python"), str(script_path)],
    capture_output=True,
    text=True,
    env=env,
    timeout=1200
)

print(result.stdout)

if result.returncode != 0:
    print("========== ERROR DETAILS ==========")
    print(result.stderr[-6000:])

print("Process return code:", result.returncode)

from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
```

### Source lines 4254–4434
```python
heavy_atoms = [
    line
    for line in atom_lines
    if line[77:79].strip().upper() not in ("H", "HD")
]

print("========== MPA PDBQT VALIDATION ==========")
print("PDBQT atom records:", len(atom_lines))
print("PDBQT heavy atoms:", len(heavy_atoms))

assert len(heavy_atoms) == 28, (
    f"Expected 28 MPA heavy atoms, found {len(heavy_atoms)}."
)

# Verify that Vina accepts the prepared ligand
# together with the existing minimized PR–EE receptor.
validation_script = work / "validate_pr_ee_mpa_vina.py"

validation_script.write_text(r'''
from vina import Vina
from pathlib import Path

work = Path("/content/PR_Dock_EM_Dock")

v = Vina(
    sf_name="vina",
    cpu=2,
    seed=2026,
    verbosity=0
)

v.set_receptor(
    str(work / "PR_EE_rigid_complex.pdbqt")
)

v.set_ligand_from_file(
    str(work / "PR_EE_to_MPA_ligand.pdbqt")
)

print("SUCCESS: Vina loaded the PR–EE receptor and MPA ligand.")
print("No docking performed.")
''')

check = subprocess.run(
    [str(python_exe), str(validation_script)],
    capture_output=True,
    text=True,
    env=env,
    timeout=120
)

print(check.stdout)

if check.returncode != 0:
    print("========== VINA ERROR DETAILS ==========")
    print(check.stderr[-5000:])
    raise RuntimeError("Vina receptor/ligand validation failed.")

print("Saved:", mpa_pdbqt)
print("SUCCESS: MPA is prepared for sequential docking.")
print("No docking performed.")

import os
import subprocess
from pathlib import Path

work = Path("/content/PR_Dock_EM_Dock")
env_path = Path("/content/pr_openff_env")
script_path = work / "dock_pr_ee_to_mpa.py"

script_path.write_text(r'''
from pathlib import Path
import numpy as np
from rdkit import Chem
from vina import Vina

work = Path("/content/PR_Dock_EM_Dock")

receptor_path = work / "PR_EE_rigid_complex.pdbqt"
ligand_path = work / "PR_EE_to_MPA_ligand.pdbqt"
ee_path = work / "PR_EE_minimized_ligand.sdf"

output_path = work / "PR_EE_to_MPA_docked.pdbqt"

for path in (receptor_path, ligand_path, ee_path):
    assert path.exists(), f"Missing file: {path}"

# Recalculate the docking-box center from minimized EE.
ee = Chem.SDMolSupplier(
    str(ee_path),
    removeHs=False
)[0]

assert ee is not None

xyz = np.asarray(
    ee.GetConformer().GetPositions()
)

heavy_indices = [
    atom.GetIdx()
    for atom in ee.GetAtoms()
    if atom.GetSymbol() != "H"
]

center = xyz[heavy_indices].mean(axis=0).tolist()

print("========== PR: EE → MPA SEQUENTIAL DOCKING ==========")
print("Receptor:", receptor_path.name)
print("Incoming ligand:", ligand_path.name)
print("Box center (Å):", [round(x, 3) for x in center])
print("Box size (Å): [30.0, 30.0, 30.0]")
print("Exhaustiveness: 16")
print("Requested poses: 10")
print("Random seed: 2026")

v = Vina(
    sf_name="vina",
    cpu=2,
    seed=2026,
    verbosity=1
)

v.set_receptor(str(receptor_path))
v.set_ligand_from_file(str(ligand_path))

v.compute_vina_maps(
    center=center,
    box_size=[30.0, 30.0, 30.0]
)

print("\n========== RUNNING DOCKING ==========")

v.dock(
    exhaustiveness=16,
    n_poses=10
)

energies = v.energies(n_poses=10)

assert len(energies) > 0, "No docking poses returned."

print("\n========== DOCKING RESULTS ==========")

for i, row in enumerate(energies, start=1):
    print(
        f"Pose {i}: "
        f"Affinity = {row[0]:.3f} kcal/mol"
    )

v.write_poses(
    str(output_path),
    n_poses=len(energies),
    overwrite=True
)

assert output_path.exists()
assert output_path.stat().st_size > 0

print("\nSaved:", output_path)
print("Best affinity (kcal/mol):",
      round(float(energies[0][0]), 3))
print("SUCCESS: PR EE → MPA docking completed.")
''')

env = os.environ.copy()
env["PATH"] = str(env_path / "bin") + os.pathsep + env.get("PATH", "")
env["MPLBACKEND"] = "Agg"

result = subprocess.run(
    [str(env_path / "bin/python"), str(script_path)],
    capture_output=True,
    text=True,
    env=env,
    timeout=1200
)

print(result.stdout)

if result.returncode != 0:
    print("========== ERROR DETAILS ==========")
```

### Source lines 5210–5391
```python

heavy_atoms = [
    line
    for line in atom_lines
    if line[77:79].strip().upper() not in ("H", "HD")
]

print("========== NGMN PDBQT VALIDATION ==========")
print("PDBQT atom records:", len(atom_lines))
print("PDBQT heavy atoms:", len(heavy_atoms))

assert len(heavy_atoms) == 24, (
    f"Expected 24 NGMN heavy atoms, found {len(heavy_atoms)}."
)

# Check receptor and ligand compatibility with Vina.
validation_script = work / "validate_pr_ee_ngmn_vina.py"

validation_script.write_text(r'''
from pathlib import Path
from vina import Vina

work = Path("/content/PR_Dock_EM_Dock")

v = Vina(
    sf_name="vina",
    cpu=2,
    seed=2026,
    verbosity=0
)

v.set_receptor(
    str(work / "PR_EE_rigid_complex.pdbqt")
)

v.set_ligand_from_file(
    str(work / "PR_EE_to_NGMN_ligand.pdbqt")
)

print("SUCCESS: Vina loaded PR–EE and NGMN.")
print("No docking performed.")
''')

check = subprocess.run(
    [str(python_exe), str(validation_script)],
    capture_output=True,
    text=True,
    env=env,
    timeout=120
)

print(check.stdout)

if check.returncode != 0:
    print("========== VINA ERROR DETAILS ==========")
    print(check.stderr[-5000:])
    raise RuntimeError("Vina compatibility check failed.")

print("Saved:", ngmn_pdbqt)
print("SUCCESS: NGMN prepared for sequential docking.")
print("No docking performed.")

import os
import subprocess
from pathlib import Path

work = Path("/content/PR_Dock_EM_Dock")
env_path = Path("/content/pr_openff_env")
script_path = work / "dock_pr_ee_to_ngmn.py"

script_path.write_text(r'''
from pathlib import Path
import numpy as np
from rdkit import Chem
from vina import Vina

work = Path("/content/PR_Dock_EM_Dock")

receptor_path = work / "PR_EE_rigid_complex.pdbqt"
ligand_path = work / "PR_EE_to_NGMN_ligand.pdbqt"
ee_path = work / "PR_EE_minimized_ligand.sdf"

output_path = work / "PR_EE_to_NGMN_docked.pdbqt"

for path in (receptor_path, ligand_path, ee_path):
    assert path.exists(), f"Missing required file: {path}"

# Calculate the search-box center from minimized EE.
ee = Chem.SDMolSupplier(
    str(ee_path),
    removeHs=False
)[0]

assert ee is not None, "Could not read minimized EE."

xyz = np.asarray(
    ee.GetConformer().GetPositions()
)

heavy_indices = [
    atom.GetIdx()
    for atom in ee.GetAtoms()
    if atom.GetSymbol() != "H"
]

center = xyz[heavy_indices].mean(axis=0).tolist()
box_size = [30.0, 30.0, 30.0]

print("========== PR: EE → NGMN SEQUENTIAL DOCKING ==========")
print("Receptor:", receptor_path.name)
print("Incoming ligand:", ligand_path.name)
print("Box center (Å):", [round(x, 3) for x in center])
print("Box size (Å):", box_size)
print("Exhaustiveness: 16")
print("Requested poses: 10")
print("Random seed: 2026")

v = Vina(
    sf_name="vina",
    cpu=2,
    seed=2026,
    verbosity=1
)

v.set_receptor(str(receptor_path))
v.set_ligand_from_file(str(ligand_path))

v.compute_vina_maps(
    center=center,
    box_size=box_size
)

print("\n========== RUNNING DOCKING ==========")

v.dock(
    exhaustiveness=16,
    n_poses=10
)

energies = v.energies(n_poses=10)

assert len(energies) > 0, "No docking poses returned."

print("\n========== DOCKING RESULTS ==========")

for i, row in enumerate(energies, start=1):
    print(
        f"Pose {i}: "
        f"Affinity = {float(row[0]):.3f} kcal/mol"
    )

v.write_poses(
    str(output_path),
    n_poses=len(energies),
    overwrite=True
)

assert output_path.exists(), "Docking output missing."
assert output_path.stat().st_size > 0, "Docking output is empty."

print("\nSaved:", output_path)
print(
    "Best affinity (kcal/mol):",
    round(float(energies[0][0]), 3)
)
print("SUCCESS: PR EE → NGMN docking completed.")
''')

env = os.environ.copy()
env["PATH"] = (
    str(env_path / "bin")
    + os.pathsep
    + env.get("PATH", "")
)
env["MPLBACKEND"] = "Agg"

result = subprocess.run(
    [str(env_path / "bin/python"), str(script_path)],
    capture_output=True,
    text=True,
    env=env,
    timeout=1200
```

### Source lines 7897–8082
```python
modeller.add(
    ligand_topology,
    ligand_positions
)

protein_atoms = list(protein.topology.atoms())
combined_atoms = list(modeller.topology.atoms())

protein_atom_count = len(protein_atoms)
ligand_atom_count = ligand.n_atoms

assert len(combined_atoms) == (
    protein_atom_count + ligand_atom_count
)

print("\nProtein atoms:", protein_atom_count)
print("Ligand atoms:", ligand_atom_count)
print("Combined atoms:", len(combined_atoms))

# AMBER14 for PR, OpenFF Sage 2.2.1 for MPA.
forcefield = app.ForceField(
    "amber14-all.xml"
)

generator = SMIRNOFFTemplateGenerator(
    molecules=[ligand],
    forcefield="openff-2.2.1"
)

forcefield.registerTemplateGenerator(
    generator.generator
)

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=None,
    rigidWater=False
)

assert system.getNumParticles() == len(combined_atoms)

print("\nOpenMM particles:",
      system.getNumParticles())
print("OpenMM forces:",
      system.getNumForces())

# Restrain only PR heavy atoms.
restraint = mm.CustomExternalForce(
    "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
)

restraint.addGlobalParameter(
    "k",
    1000.0 * unit.kilojoule_per_mole
    / unit.nanometer**2
)

for name in ("x0", "y0", "z0"):
    restraint.addPerParticleParameter(name)

initial_xyz_nm = np.asarray(
    modeller.positions.value_in_unit(
        unit.nanometer
    )
)

restrained_indices = []

for atom in protein_atoms:
    if (
        atom.element is not None
        and atom.element.symbol != "H"
    ):
        idx = atom.index

        restraint.addParticle(
            idx,
            initial_xyz_nm[idx].tolist()
        )

        restrained_indices.append(idx)

assert len(restrained_indices) == 2036

system.addForce(restraint)

print("Restrained PR heavy atoms:",
      len(restrained_indices))
print("Restrained ligand atoms: 0")

integrator = mm.VerletIntegrator(
    0.001 * unit.picoseconds
)

platform = mm.Platform.getPlatformByName("CPU")

simulation = app.Simulation(
    modeller.topology,
    system,
    integrator,
    platform
)

simulation.context.setPositions(
    modeller.positions
)

initial_state = simulation.context.getState(
    getEnergy=True
)

initial_energy = initial_state.getPotentialEnergy().value_in_unit(
    unit.kilojoule_per_mole
)

print("\nInitial potential energy (kJ/mol):",
      round(initial_energy, 3))

print("Minimizing (maximum 1000 iterations)...")

simulation.minimizeEnergy(
    maxIterations=1000
)

final_state = simulation.context.getState(
    getEnergy=True,
    getPositions=True
)

final_energy = final_state.getPotentialEnergy().value_in_unit(
    unit.kilojoule_per_mole
)

final_positions = final_state.getPositions(
    asNumpy=True
)

final_xyz_nm = final_positions.value_in_unit(
    unit.nanometer
)

# Displacement measurements.
initial_xyz_a = initial_xyz_nm * 10.0
final_xyz_a = final_xyz_nm * 10.0

protein_shift = (
    final_xyz_a[:protein_atom_count]
    - initial_xyz_a[:protein_atom_count]
)

ligand_shift = (
    final_xyz_a[protein_atom_count:]
    - initial_xyz_a[protein_atom_count:]
)

protein_heavy_indices = np.asarray(
    restrained_indices,
    dtype=int
)

ligand_heavy_indices = np.asarray([
    atom.molecule_atom_index
    for atom in ligand.atoms
    if atom.symbol != "H"
], dtype=int)

protein_rmsd = np.sqrt(np.mean(np.sum(
    protein_shift[protein_heavy_indices]**2,
    axis=1
)))

ligand_rmsd = np.sqrt(np.mean(np.sum(
    ligand_shift[ligand_heavy_indices]**2,
    axis=1
)))

initial_ligand_centroid = initial_xyz_a[
    protein_atom_count + ligand_heavy_indices
].mean(axis=0)

final_ligand_centroid = final_xyz_a[
    protein_atom_count + ligand_heavy_indices
].mean(axis=0)

centroid_shift = np.linalg.norm(
```

### Source lines 8579–8667
```python
    removeHs=False
)[0]

assert mpa is not None

mpa_xyz = np.asarray(
    mpa.GetConformer().GetPositions()
)

heavy_indices = [
    atom.GetIdx()
    for atom in mpa.GetAtoms()
    if atom.GetSymbol() != "H"
]

center = mpa_xyz[heavy_indices].mean(axis=0)

print("\nFixed MPA centroid (Å):",
      np.round(center, 3).tolist())

print("Grid dimensions (Å): [30, 30, 30]")

print("\nInitializing Vina...")

vina = Vina(
    sf_name="vina",
    cpu=2,
    seed=42
)

vina.set_receptor(
    str(receptor_path)
)

vina.compute_vina_maps(
    center=center.tolist(),
    box_size=[30.0, 30.0, 30.0]
)

print("\nSUCCESS: Vina receptor grid generated.")
print("PR–MPA rigid receptor is readable.")
print("No docking performed.")
''')

env = os.environ.copy()
env["PATH"] = (
    str(env_path / "bin")
    + os.pathsep
    + env.get("PATH", "")
)
env["MPLBACKEND"] = "Agg"

result = subprocess.run(
    [str(env_path / "bin/python"), str(script_path)],
    capture_output=True,
    text=True,
    env=env,
    timeout=300
)

print(result.stdout)

if result.returncode != 0:
    print("========== ERROR DETAILS ==========")
    print(result.stderr[-6000:])
    raise RuntimeError(
        "PR–MPA Vina compatibility check failed."
    )

print("Process return code:", result.returncode)

import os
import subprocess
from pathlib import Path

work = Path("/content/PR_Dock_EM_Dock")
env_path = Path("/content/pr_openff_env")

ee_sdf = Path("/content/ee.sdf")
ee_pdbqt = work / "PR_MPA_to_EE_ligand.pdbqt"

assert ee_sdf.exists(), "Original EE SDF is missing."

env = os.environ.copy()
env["PATH"] = (
    str(env_path / "bin")
    + os.pathsep
    + env.get("PATH", "")
)
```

### Source lines 8817–8905
```python
    str(fixed_mpa_path),
    removeHs=False
)[0]

assert fixed_mpa is not None

xyz = np.asarray(
    fixed_mpa.GetConformer().GetPositions()
)

heavy_indices = [
    atom.GetIdx()
    for atom in fixed_mpa.GetAtoms()
    if atom.GetSymbol() != "H"
]

center = xyz[heavy_indices].mean(axis=0)

print("Receptor:", receptor_path.name)
print("Incoming ligand:", incoming_path.name)
print("Grid center (Å):",
      np.round(center, 3).tolist())
print("Grid size (Å): [30.0, 30.0, 30.0]")

v = Vina(
    sf_name="vina",
    cpu=2,
    seed=42
)

v.set_receptor(str(receptor_path))

v.set_ligand_from_file(
    str(incoming_path)
)

v.compute_vina_maps(
    center=center.tolist(),
    box_size=[30.0, 30.0, 30.0]
)

print("\nRunning Vina docking...")
print("Exhaustiveness: 32")
print("Requested poses: 20")

v.dock(
    exhaustiveness=32,
    n_poses=20
)

energies = np.asarray(
    v.energies(n_poses=20)
)

assert len(energies) > 0
assert np.isfinite(energies).all()

v.write_poses(
    str(output_path),
    n_poses=20,
    overwrite=True
)

assert output_path.exists()

print("\n========== DOCKING RESULTS ==========")

print("Returned poses:", len(energies))
print("Best Vina affinity (kcal/mol):",
      round(float(energies[0, 0]), 3))

print("\nPose | Affinity (kcal/mol)")
print("-" * 30)

for idx, row in enumerate(energies, start=1):
    print(
        f"{idx:>4} | {row[0]:>8.3f}"
    )

print("\nSaved:", output_path)

print("\nSUCCESS: MPA → EE docking completed.")
print("Fixed PR–MPA receptor was not minimized.")
''')

env = os.environ.copy()
env["PATH"] = (
    str(env_path / "bin")
    + os.pathsep
```

### Source lines 9996–10084
```python
    str(fixed_mpa_path),
    removeHs=False
)[0]

assert mpa is not None

positions = np.asarray(
    mpa.GetConformer().GetPositions()
)

heavy_indices = [
    atom.GetIdx()
    for atom in mpa.GetAtoms()
    if atom.GetSymbol() != "H"
]

center = positions[heavy_indices].mean(axis=0)

print("Rigid receptor:", receptor_path.name)
print("Incoming ligand:", incoming_path.name)
print("Grid center (Å):",
      np.round(center, 3).tolist())
print("Grid size (Å): [30.0, 30.0, 30.0]")

v = Vina(
    sf_name="vina",
    cpu=2,
    seed=42
)

v.set_receptor(str(receptor_path))
v.set_ligand_from_file(str(incoming_path))

v.compute_vina_maps(
    center=center.tolist(),
    box_size=[30.0, 30.0, 30.0]
)

print("\nRunning docking...")
print("Exhaustiveness: 32")
print("Requested poses: 20")

v.dock(
    exhaustiveness=32,
    n_poses=20
)

# Request up to 20 pose energies.
energies = np.asarray(
    v.energies(n_poses=20)
)

assert energies.ndim == 2
assert len(energies) >= 1
assert np.isfinite(energies).all()

v.write_poses(
    str(output_path),
    n_poses=20,
    overwrite=True
)

assert output_path.exists()

saved_lines = output_path.read_text(
    errors="replace"
).splitlines()

model_count = sum(
    line.startswith("MODEL")
    for line in saved_lines
)

if model_count == 0:
    model_count = 1 if any(
        line.startswith(("ATOM  ", "HETATM"))
        for line in saved_lines
    ) else 0

assert model_count >= 1

print("\n========== DOCKING RESULTS ==========")

print("Energy rows returned:", len(energies))
print("Saved PDBQT models:", model_count)

print(
    "Best Vina affinity (kcal/mol):",
    round(float(energies[0, 0]), 3)
```

### Source lines 11195–11283
```python
    str(fixed_mpa_path),
    removeHs=False
)[0]

assert mpa is not None

positions = np.asarray(
    mpa.GetConformer().GetPositions()
)

heavy_indices = [
    atom.GetIdx()
    for atom in mpa.GetAtoms()
    if atom.GetSymbol() != "H"
]

center = positions[heavy_indices].mean(axis=0)

print("Rigid receptor:", receptor_path.name)
print("Incoming ligand:", incoming_path.name)
print("Grid center (Å):",
      np.round(center, 3).tolist())
print("Grid size (Å): [30.0, 30.0, 30.0]")

v = Vina(
    sf_name="vina",
    cpu=2,
    seed=42
)

v.set_receptor(str(receptor_path))
v.set_ligand_from_file(str(incoming_path))

v.compute_vina_maps(
    center=center.tolist(),
    box_size=[30.0, 30.0, 30.0]
)

print("\nRunning docking...")
print("Exhaustiveness: 32")
print("Requested poses: 20")

v.dock(
    exhaustiveness=32,
    n_poses=20
)

energies = np.asarray(
    v.energies(n_poses=20)
)

assert energies.ndim == 2
assert len(energies) >= 1
assert np.isfinite(energies).all()

v.write_poses(
    str(output_path),
    n_poses=20,
    overwrite=True
)

assert output_path.exists()

saved_lines = output_path.read_text(
    errors="replace"
).splitlines()

model_count = sum(
    line.startswith("MODEL")
    for line in saved_lines
)

if model_count == 0:
    model_count = 1 if any(
        line.startswith(("ATOM  ", "HETATM"))
        for line in saved_lines
    ) else 0

assert model_count >= 1

print("\n========== DOCKING RESULTS ==========")

print("Energy rows returned:", len(energies))
print("Saved PDBQT models:", model_count)

print(
    "Best Vina affinity (kcal/mol):",
    round(float(energies[0, 0]), 3)
)
```

### Source lines 15280–15368
```python
incoming = work / "PR_LNG_to_MPA_incoming.pdbqt"
fixed_sdf = work / "PR_LNG_minimized_ligand.sdf"
output = work / "PR_LNG_to_MPA_docked.pdbqt"
log = work / "PR_LNG_to_MPA_docking_log.txt"

fixed = Chem.SDMolSupplier(str(fixed_sdf), removeHs=False)[0]
assert fixed is not None

heavy = [
    a.GetIdx() for a in fixed.GetAtoms()
    if a.GetAtomicNum() > 1
]

xyz = np.asarray(fixed.GetConformer().GetPositions())
center = xyz[heavy].mean(axis=0)

print("========== PR LNG→MPA DOCKING ==========")
print("Grid center (Å):", np.round(center, 4).tolist())
print("Grid size (Å): [30, 30, 30]")
print("Exhaustiveness: 32")
print("Seed: 42")
print("CPU: 2")
print("Requested poses: 20")

v = Vina(sf_name="vina", cpu=2, seed=42, verbosity=1)

v.set_receptor(rigid_pdbqt_filename=str(receptor))
v.set_ligand_from_file(str(incoming))

v.compute_vina_maps(
    center=center.tolist(),
    box_size=[30.0, 30.0, 30.0]
)

v.dock(exhaustiveness=32, n_poses=20)

energies = v.energies(n_poses=20)
assert len(energies) > 0

v.write_poses(
    pdbqt_filename=str(output),
    n_poses=20,
    overwrite=True
)

saved = sum(
    line.startswith("MODEL ")
    for line in output.read_text().splitlines()
)

if saved == 0 and "ATOM" in output.read_text():
    saved = 1

assert saved > 0

print("\n========== DOCKING RESULTS ==========")
print("Returned energy rows:", len(energies))
print("Saved poses:", saved)
print("Best Vina score (kcal/mol):", round(float(energies[0][0]), 4))

for i, row in enumerate(energies, 1):
    print(f"Pose {i:>2}: {float(row[0]):.4f} kcal/mol")

log.write_text(
    "PR LNG→MPA sequential docking\n"
    f"Grid center (Å): {center.tolist()}\n"
    "Grid size (Å): [30, 30, 30]\n"
    "Exhaustiveness: 32\n"
    "Seed: 42\n"
    "CPU: 2\n"
    "Requested poses: 20\n"
    f"Returned energy rows: {len(energies)}\n"
    f"Saved poses: {saved}\n"
    f"Best score (kcal/mol): {float(energies[0][0]):.6f}\n"
)

print("\nSaved docking poses:", output)
print("Saved docking log:", log)
print("SUCCESS: PR LNG→MPA docking completed.")
print("Fixed PR–LNG receptor unchanged.")
''')

result = subprocess.run(
    [str(env_path / "bin/python"), str(script)],
    capture_output=True,
    text=True,
    env=env,
    timeout=1800
)
```

### Source lines 15894–15982
```python
)[0]

assert fixed is not None

heavy_indices = [
    atom.GetIdx()
    for atom in fixed.GetAtoms()
    if atom.GetAtomicNum() > 1
]

coords = np.asarray(
    fixed.GetConformer().GetPositions()
)

center = coords[heavy_indices].mean(axis=0)

print("========== PR LNG→NGMN DOCKING ==========")
print("Grid center (Å):", np.round(center, 4).tolist())
print("Grid size (Å): [30, 30, 30]")
print("Exhaustiveness: 32")
print("Seed: 42")
print("CPU: 2")
print("Requested poses: 20")

v = Vina(
    sf_name="vina",
    cpu=2,
    seed=42,
    verbosity=1
)

v.set_receptor(
    rigid_pdbqt_filename=str(receptor)
)

v.set_ligand_from_file(str(incoming))

v.compute_vina_maps(
    center=center.tolist(),
    box_size=[30.0, 30.0, 30.0]
)

v.dock(
    exhaustiveness=32,
    n_poses=20
)

energies = v.energies(n_poses=20)
assert len(energies) > 0, "No docking poses returned."

v.write_poses(
    pdbqt_filename=str(output),
    n_poses=20,
    overwrite=True
)

saved = sum(
    line.startswith("MODEL ")
    for line in output.read_text().splitlines()
)

if saved == 0 and "ATOM" in output.read_text():
    saved = 1

assert saved > 0, "No poses saved."

print("\n========== DOCKING RESULTS ==========")
print("Returned energy rows:", len(energies))
print("Saved poses:", saved)
print(
    "Best local Vina score (kcal/mol):",
    round(float(energies[0][0]), 4)
)

for i, row in enumerate(energies, 1):
    print(
        f"Pose {i:>2}: "
        f"{float(row[0]):.4f} kcal/mol"
    )

log.write_text(
    "PR LNG→NGMN sequential docking\n"
    f"Grid center (Å): {center.tolist()}\n"
    "Grid size (Å): [30, 30, 30]\n"
    "Exhaustiveness: 32\n"
    "Seed: 42\n"
    "CPU: 2\n"
    "Requested poses: 20\n"
    f"Returned energy rows: {len(energies)}\n"
```

### Source lines 16156–16244
```python
fixed = Chem.SDMolSupplier(
    str(work / "PR_LNG_minimized_ligand.sdf"),
    removeHs=False
)[0]

assert fixed is not None

heavy = [
    a.GetIdx()
    for a in fixed.GetAtoms()
    if a.GetAtomicNum() > 1
]

xyz = np.asarray(fixed.GetConformer().GetPositions())
center = xyz[heavy].mean(axis=0)

print("========== PR LNG→NGMN DOCKING ==========")
print("Grid center (Å):", np.round(center, 4).tolist())
print("Grid size (Å): [30, 30, 30]")
print("Exhaustiveness: 32")
print("Seed: 42")
print("CPU: 2")
print("Requested poses: 20")

v = Vina(
    sf_name="vina",
    cpu=2,
    seed=42,
    verbosity=1
)

v.set_receptor(
    rigid_pdbqt_filename=str(
        work / "PR_LNG_rigid_complex.pdbqt"
    )
)

v.set_ligand_from_file(
    str(work / "PR_LNG_to_NGMN_incoming.pdbqt")
)

v.compute_vina_maps(
    center=center.tolist(),
    box_size=[30.0, 30.0, 30.0]
)

v.dock(exhaustiveness=32, n_poses=20)

energies = v.energies(n_poses=20)
assert len(energies) > 0

output = work / "PR_LNG_to_NGMN_docked.pdbqt"

v.write_poses(
    pdbqt_filename=str(output),
    n_poses=20,
    overwrite=True
)

saved = sum(
    line.startswith("MODEL ")
    for line in output.read_text().splitlines()
)

if saved == 0 and "ATOM" in output.read_text():
    saved = 1

assert saved > 0

print("\n========== RESULTS ==========")
print("Returned energy rows:", len(energies))
print("Saved poses:", saved)
print(
    "Best local Vina score (kcal/mol):",
    round(float(energies[0][0]), 4)
)

for i, row in enumerate(energies, 1):
    print(
        f"Pose {i:>2}: "
        f"{float(row[0]):.4f} kcal/mol"
    )

log = work / "PR_LNG_to_NGMN_docking_log.txt"

log.write_text(
    "PR LNG→NGMN sequential docking\n"
    f"Grid center (Å): {center.tolist()}\n"
    "Grid size (Å): [30, 30, 30]\n"
```

### Source lines 17629–17928
```python

hydrogen_count = sum(
    atom.GetAtomicNum() == 1
    for atom in saved.GetAtoms()
)

print("Reference heavy atoms:", len(mapped_xyz))
print("Reconstructed heavy atoms:", len(heavy_indices))
print("Reconstructed hydrogens:", hydrogen_count)
print("Total ligand atoms:", saved.GetNumAtoms())

print("\n========== VALIDATION ==========")
print("Expected CIP stereocenters:", expected_cip)
print("Observed CIP stereocenters:", observed_cip)
print("Stereocenters matched:", len(expected_cip), "/", len(expected_cip))
print(
    "Maximum heavy-atom shift during H addition (Å):",
    round(max_heavy_shift, 8)
)
print(
    "Maximum heavy-atom shift after SDF save (Å):",
    round(save_roundtrip_shift, 8)
)
print(
    "Minimum nonbonded heavy-atom distance (Å):",
    round(minimum_nonbonded_distance, 4)
)

print("\nSaved:", output_path)
print("SUCCESS: PR–NGMN ligand reconstruction completed.")
print("Original docked heavy-atom coordinates preserved.")
print("No energy minimization performed.")

from pathlib import Path
import os
import subprocess

work = Path("/content/PR_Dock_EM_Dock")
venv = Path("/content/pr_openff_env")

script_path = work / "minimize_pr_ngmn_complex.py"

script_path.write_text(r'''
from pathlib import Path
import numpy as np
from rdkit import Chem

from openmm import (
    CustomExternalForce,
    LangevinMiddleIntegrator,
    Platform,
    unit
)
from openmm.app import (
    PDBFile,
    Modeller,
    ForceField,
    NoCutoff,
    Simulation
)

from openff.toolkit import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

work = Path("/content/PR_Dock_EM_Dock")

protein_path = work / "PR_NGMN_repaired_protein.pdb"
ligand_path = work / "PR_NGMN_reconstructed_ligand.sdf"

protein_output = work / "PR_NGMN_minimized_protein.pdb"
ligand_output = work / "PR_NGMN_minimized_ligand.sdf"
complex_output = work / "PR_NGMN_minimized_complex.pdb"

for path in (protein_path, ligand_path):
    assert path.exists(), f"Missing input: {path}"

print("========== PR NGMN RESTRAINED MINIMIZATION ==========")

protein = PDBFile(str(protein_path))

rdmol = Chem.SDMolSupplier(
    str(ligand_path),
    removeHs=False
)[0]

assert rdmol is not None
assert rdmol.GetNumAtoms() == 53

offmol = Molecule.from_rdkit(
    rdmol,
    allow_undefined_stereo=False,
    hydrogens_are_explicit=True
)

# AM1-BCC partial charges, OpenFF Sage 2.2.1 parameters.
offmol.assign_partial_charges(
    partial_charge_method="am1bcc"
)

generator = SMIRNOFFTemplateGenerator(
    molecules=[offmol],
    forcefield="openff-2.2.1"
)

forcefield = ForceField("amber14-all.xml")
forcefield.registerTemplateGenerator(
    generator.generator
)

modeller = Modeller(
    protein.topology,
    protein.positions
)

ligand_topology = offmol.to_topology().to_openmm()
ligand_positions = (
    np.asarray(rdmol.GetConformer().GetPositions())
    * unit.angstrom
)

modeller.add(
    ligand_topology,
    ligand_positions
)

protein_atom_count = protein.topology.getNumAtoms()
ligand_atom_count = rdmol.GetNumAtoms()

print("Protein atoms:", protein_atom_count)
print("Ligand atoms:", ligand_atom_count)
print("Combined atoms:", modeller.topology.getNumAtoms())

assert protein_atom_count == 4138
assert ligand_atom_count == 53

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=NoCutoff,
    constraints=None,
    rigidWater=False
)

# Harmonic positional restraints on protein heavy atoms.
# 1000 kJ/mol/nm²; protein hydrogens and NGMN are free.
restraint_k = 1000.0

restraints = CustomExternalForce(
    "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
)

restraints.addGlobalParameter(
    "k",
    restraint_k * unit.kilojoule_per_mole / unit.nanometer**2
)

for parameter in ("x0", "y0", "z0"):
    restraints.addPerParticleParameter(parameter)

initial_positions = modeller.positions.value_in_unit(
    unit.nanometer
)

protein_heavy_indices = []

for atom in modeller.topology.atoms():
    if atom.index >= protein_atom_count:
        continue

    if atom.element.symbol == "H":
        continue

    protein_heavy_indices.append(atom.index)

    restraints.addParticle(
        atom.index,
        [
            float(value)
            for value in initial_positions[atom.index]
        ]
    )

system.addForce(restraints)

integrator = LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1 / unit.picosecond,
    0.002 * unit.picoseconds
)

platform = Platform.getPlatformByName("Reference")

simulation = Simulation(
    modeller.topology,
    system,
    integrator,
    platform
)

simulation.context.setPositions(modeller.positions)

initial_energy = simulation.context.getState(
    getEnergy=True
).getPotentialEnergy().value_in_unit(
    unit.kilojoule_per_mole
)

print("Restrained PR heavy atoms:", len(protein_heavy_indices))
print("Initial potential energy (kJ/mol):", round(initial_energy, 4))

simulation.minimizeEnergy(
    tolerance=10 * unit.kilojoule_per_mole / unit.nanometer,
    maxIterations=2000
)

state = simulation.context.getState(
    getPositions=True,
    getEnergy=True
)

final_energy = state.getPotentialEnergy().value_in_unit(
    unit.kilojoule_per_mole
)

final_positions = state.getPositions(
    asNumpy=True
).value_in_unit(unit.angstrom)

initial_xyz = np.asarray(
    modeller.positions.value_in_unit(unit.angstrom)
)

pr_indices = np.asarray(
    protein_heavy_indices,
    dtype=int
)

ligand_indices = np.arange(
    protein_atom_count,
    protein_atom_count + ligand_atom_count
)

ligand_heavy_local = np.asarray([
    atom.GetIdx()
    for atom in rdmol.GetAtoms()
    if atom.GetAtomicNum() > 1
])

ligand_heavy_global = (
    protein_atom_count + ligand_heavy_local
)

def rmsd(a, b):
    return float(
        np.sqrt(
            np.mean(
                np.sum((a - b)**2, axis=1)
            )
        )
    )

protein_rmsd = rmsd(
    initial_xyz[pr_indices],
    final_positions[pr_indices]
)

ligand_rmsd = rmsd(
    initial_xyz[ligand_heavy_global],
    final_positions[ligand_heavy_global]
)

initial_centroid = initial_xyz[
    ligand_heavy_global
].mean(axis=0)

final_centroid = final_positions[
    ligand_heavy_global
].mean(axis=0)

centroid_shift = float(
    np.linalg.norm(final_centroid - initial_centroid)
)

# Save minimized protein.
protein_only = Modeller(
    modeller.topology,
    state.getPositions()
)

protein_only.delete([
    residue
    for residue in protein_only.topology.residues()
    if residue.name not in {
        residue.name
        for residue in protein.topology.residues()
    }
])

# Write protein directly from its original topology,
# using the corresponding minimized coordinates.
with open(protein_output, "w") as handle:
```

### Source lines 18430–18729
```python
    removeHs=False,
    sanitize=True
)[0]

assert rdmol is not None
assert rdmol.GetNumAtoms() == 53

oxime = rdmol.GetBondBetweenAtoms(2, 22)
assert oxime is not None
assert oxime.GetStereo() == Chem.BondStereo.STEREOE

offmol = Molecule.from_rdkit(
    rdmol,
    allow_undefined_stereo=False,
    hydrogens_are_explicit=True
)

print("Strict OpenFF stereochemistry validation: PASSED")

# Match the previous PR minimization protocol.
offmol.assign_partial_charges(
    partial_charge_method="am1bcc"
)

generator = SMIRNOFFTemplateGenerator(
    molecules=[offmol],
    forcefield="openff-2.2.1"
)

forcefield = ForceField("amber14-all.xml")
forcefield.registerTemplateGenerator(generator.generator)

modeller = Modeller(
    protein.topology,
    protein.positions
)

ligand_topology = offmol.to_topology().to_openmm()

ligand_xyz = np.asarray(
    rdmol.GetConformer().GetPositions(),
    dtype=float
)

modeller.add(
    ligand_topology,
    ligand_xyz * unit.angstrom
)

protein_atom_count = protein.topology.getNumAtoms()
ligand_atom_count = rdmol.GetNumAtoms()

assert protein_atom_count == 4138
assert ligand_atom_count == 53
assert modeller.topology.getNumAtoms() == 4191

print("Protein atoms:", protein_atom_count)
print("Ligand atoms:", ligand_atom_count)
print("Combined atoms:", modeller.topology.getNumAtoms())

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=NoCutoff,
    constraints=None,
    rigidWater=False
)

# Restrain PR heavy atoms only.
restraint_k = 1000.0

restraints = CustomExternalForce(
    "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
)

restraints.addGlobalParameter(
    "k",
    restraint_k
    * unit.kilojoule_per_mole
    / unit.nanometer**2
)

for parameter in ("x0", "y0", "z0"):
    restraints.addPerParticleParameter(parameter)

initial_nm = np.asarray(
    modeller.positions.value_in_unit(unit.nanometer)
)

protein_heavy_indices = []

for atom in modeller.topology.atoms():
    if atom.index >= protein_atom_count:
        continue
    if atom.element.symbol == "H":
        continue

    protein_heavy_indices.append(atom.index)

    restraints.addParticle(
        atom.index,
        initial_nm[atom.index].tolist()
    )

system.addForce(restraints)

integrator = LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1 / unit.picosecond,
    0.002 * unit.picoseconds
)

platform = Platform.getPlatformByName("Reference")

simulation = Simulation(
    modeller.topology,
    system,
    integrator,
    platform
)

simulation.context.setPositions(modeller.positions)

initial_energy = (
    simulation.context.getState(getEnergy=True)
    .getPotentialEnergy()
    .value_in_unit(unit.kilojoule_per_mole)
)

print("Restrained PR heavy atoms:", len(protein_heavy_indices))
print(
    "Initial potential energy (kJ/mol):",
    round(initial_energy, 4)
)

simulation.minimizeEnergy(
    tolerance=10
    * unit.kilojoule_per_mole
    / unit.nanometer,
    maxIterations=2000
)

state = simulation.context.getState(
    getPositions=True,
    getEnergy=True
)

final_energy = (
    state.getPotentialEnergy()
    .value_in_unit(unit.kilojoule_per_mole)
)

final_positions = np.asarray(
    state.getPositions(asNumpy=True).value_in_unit(
        unit.angstrom
    )
)

initial_positions = np.asarray(
    modeller.positions.value_in_unit(unit.angstrom)
)

protein_heavy_indices = np.asarray(
    protein_heavy_indices,
    dtype=int
)

ligand_heavy_local = np.asarray([
    atom.GetIdx()
    for atom in rdmol.GetAtoms()
    if atom.GetAtomicNum() > 1
], dtype=int)

ligand_heavy_global = (
    protein_atom_count + ligand_heavy_local
)

def rmsd(a, b):
    return float(
        np.sqrt(np.mean(np.sum((a - b)**2, axis=1)))
    )

protein_rmsd = rmsd(
    initial_positions[protein_heavy_indices],
    final_positions[protein_heavy_indices]
)

ligand_rmsd = rmsd(
    initial_positions[ligand_heavy_global],
    final_positions[ligand_heavy_global]
)

initial_centroid = initial_positions[
    ligand_heavy_global
].mean(axis=0)

final_centroid = final_positions[
    ligand_heavy_global
].mean(axis=0)

centroid_shift = float(
    np.linalg.norm(final_centroid - initial_centroid)
)

# Evaluate PR–NGMN heavy-atom proximity.
initial_min_distance = float(
    cdist(
        initial_positions[protein_heavy_indices],
        initial_positions[ligand_heavy_global]
    ).min()
)

final_min_distance = float(
    cdist(
        final_positions[protein_heavy_indices],
        final_positions[ligand_heavy_global]
    ).min()
)

assert np.isfinite(final_energy)
assert np.isfinite(final_positions).all()
assert final_energy <= initial_energy + 1e-3

# Validate the result before writing output files.
assert protein_rmsd < 2.0, (
    f"Excessive PR displacement: {protein_rmsd:.4f} Å"
)
assert ligand_rmsd < 5.0, (
    f"Excessive NGMN displacement: {ligand_rmsd:.4f} Å"
)
assert final_min_distance > 1.5, (
    f"Severe PR–NGMN heavy-atom overlap: "
    f"{final_min_distance:.4f} Å"
)

# Save minimized protein.
with open(protein_output, "w") as handle:
    PDBFile.writeFile(
        protein.topology,
        final_positions[:protein_atom_count] * unit.angstrom,
        handle,
        keepIds=True
    )

# Save minimized NGMN with validated connectivity.
minimized_ligand = Chem.Mol(rdmol)
ligand_conf = minimized_ligand.GetConformer()

for local_idx in range(ligand_atom_count):
    position = final_positions[
        protein_atom_count + local_idx
    ]
    ligand_conf.SetAtomPosition(
        local_idx,
        tuple(float(x) for x in position)
    )

writer = Chem.SDWriter(str(ligand_output))
writer.write(minimized_ligand)
writer.close()

# Save the full minimized complex.
with open(complex_output, "w") as handle:
    PDBFile.writeFile(
        modeller.topology,
        state.getPositions(),
        handle,
        keepIds=True
    )

print("\n========== MINIMIZATION RESULTS ==========")
print(
    "Final potential energy (kJ/mol):",
    round(final_energy, 4)
)
print(
    "Energy change (kJ/mol):",
    round(final_energy - initial_energy, 4)
)
print(
    "PR heavy-atom RMSD (Å):",
    round(protein_rmsd, 4)
)
print(
    "NGMN heavy-atom RMSD (Å):",
    round(ligand_rmsd, 4)
)
print(
    "NGMN centroid shift (Å):",
    round(centroid_shift, 4)
)
print(
    "Minimized NGMN centroid (Å):",
    np.round(final_centroid, 4)
)
print(
    "Initial minimum PR–NGMN distance (Å):",
    round(initial_min_distance, 4)
)
print(
    "Final minimum PR–NGMN distance (Å):",
```

### Source lines 20074–20162
```python
    "minimum_interligand_heavy_distance_A": minimum_interligand_distance,
    "interligand_heavy_atom_pairs_under_4A": pairs_under_4A,
    "interligand_heavy_atom_pairs_under_5A": pairs_under_5A,
    "contact_cutoff_A": contact_cutoff,
    "fixed_ligand_PR_contact_residues": sorted_residues(fixed_contacts),
    "incoming_ligand_PR_contact_residues": sorted_residues(lng_contacts),
    "shared_PR_contact_residues": sorted_residues(shared_contacts),
    "union_PR_contact_residue_count": len(union_contacts),
    "PR_contact_jaccard": jaccard,
    "minimum_fixed_ligand_PR_distance_A": minimum_ngmn_pr_distance,
    "minimum_incoming_ligand_PR_distance_A": minimum_lng_pr_distance
}

result_path.write_text(
    json.dumps(results, indent=2) + "\n"
)

print("Fixed NGMN heavy atoms:", len(fixed_xyz))
print("Incoming LNG heavy atoms:", len(lng_xyz))
print("Best local Vina score (kcal/mol):", best_score)

print("\n========== SPATIAL GEOMETRY ==========")
print("NGMN centroid (Å):", np.round(fixed_centroid, 4))
print("LNG centroid (Å):", np.round(lng_centroid, 4))
print("Centroid distance (Å):", round(centroid_distance, 4))
print(
    "Minimum interligand heavy distance (Å):",
    round(minimum_interligand_distance, 4)
)
print("Interligand heavy-atom pairs <4 Å:", pairs_under_4A)
print("Interligand heavy-atom pairs <5 Å:", pairs_under_5A)

print("\n========== PR CONTACT ANALYSIS ==========")
print("NGMN PR contact residues:", len(fixed_contacts))
print(sorted_residues(fixed_contacts))
print("LNG PR contact residues:", len(lng_contacts))
print(sorted_residues(lng_contacts))
print("Shared PR contact residues:", len(shared_contacts))
print(sorted_residues(shared_contacts))
print("PR contact Jaccard:", round(jaccard, 4))
print(
    "Minimum NGMN–PR distance (Å):",
    round(minimum_ngmn_pr_distance, 4)
)
print(
    "Minimum LNG–PR distance (Å):",
    round(minimum_lng_pr_distance, 4)
)

print("\nSaved:", result_path)
print("SUCCESS: PR NGMN→LNG geometry analysis completed.")
print("SASA analysis not yet performed.")

from pathlib import Path
import os
import subprocess

work = Path("/content/PR_Dock_EM_Dock")
venv = Path("/content/pr_openff_env")
script_path = work / "analyze_pr_ngmn_to_lng_sasa.py"

script_path.write_text(r'''
from pathlib import Path
import json
import numpy as np
import freesasa

work = Path("/content/PR_Dock_EM_Dock")

protein_path = work / "PR_NGMN_minimized_protein.pdb"
docked_path = work / "PR_NGMN_to_LNG_docked.pdbqt"
results_path = work / "PR_NGMN_to_LNG_results.json"

for path in (protein_path, docked_path, results_path):
    assert path.exists(), f"Missing input: {path}"

results = json.loads(results_path.read_text())

print("========== PR NGMN→LNG SASA ANALYSIS ==========")

bondi = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    "F": 1.47,
    "CL": 1.75,
    "BR": 1.85,
```

### Source lines 20819–20907
```python
    "minimum_interligand_heavy_distance_A": minimum_interligand_distance,
    "interligand_heavy_atom_pairs_under_4A": pairs_under_4A,
    "interligand_heavy_atom_pairs_under_5A": pairs_under_5A,
    "contact_cutoff_A": contact_cutoff,
    "fixed_ligand_PR_contact_residues": sorted_residues(fixed_contacts),
    "incoming_ligand_PR_contact_residues": sorted_residues(mpa_contacts),
    "shared_PR_contact_residues": sorted_residues(shared_contacts),
    "union_PR_contact_residue_count": len(union_contacts),
    "PR_contact_jaccard": jaccard,
    "minimum_fixed_ligand_PR_distance_A": minimum_ngmn_pr_distance,
    "minimum_incoming_ligand_PR_distance_A": minimum_mpa_pr_distance
}

result_path.write_text(
    json.dumps(results, indent=2) + "\n"
)

print("Fixed NGMN heavy atoms:", len(fixed_xyz))
print("Incoming MPA heavy atoms:", len(mpa_xyz))
print("Best local Vina score (kcal/mol):", best_score)

print("\n========== SPATIAL GEOMETRY ==========")
print("NGMN centroid (Å):", np.round(fixed_centroid, 4))
print("MPA centroid (Å):", np.round(mpa_centroid, 4))
print("Centroid distance (Å):", round(centroid_distance, 4))
print(
    "Minimum interligand heavy distance (Å):",
    round(minimum_interligand_distance, 4)
)
print("Interligand heavy-atom pairs <4 Å:", pairs_under_4A)
print("Interligand heavy-atom pairs <5 Å:", pairs_under_5A)

print("\n========== PR CONTACT ANALYSIS ==========")
print("NGMN PR contact residues:", len(fixed_contacts))
print(sorted_residues(fixed_contacts))
print("MPA PR contact residues:", len(mpa_contacts))
print(sorted_residues(mpa_contacts))
print("Shared PR contact residues:", len(shared_contacts))
print(sorted_residues(shared_contacts))
print("PR contact Jaccard:", round(jaccard, 4))
print(
    "Minimum NGMN–PR distance (Å):",
    round(minimum_ngmn_pr_distance, 4)
)
print(
    "Minimum MPA–PR distance (Å):",
    round(minimum_mpa_pr_distance, 4)
)

print("\nSaved:", result_path)
print("SUCCESS: PR NGMN→MPA geometry analysis completed.")
print("SASA analysis not yet performed.")

from pathlib import Path
import os
import subprocess

work = Path("/content/PR_Dock_EM_Dock")
venv = Path("/content/pr_openff_env")
script_path = work / "analyze_pr_ngmn_to_mpa_sasa.py"

script_path.write_text(r'''
from pathlib import Path
import json
import freesasa

work = Path("/content/PR_Dock_EM_Dock")

protein_path = work / "PR_NGMN_minimized_protein.pdb"
docked_path = work / "PR_NGMN_to_MPA_docked.pdbqt"
results_path = work / "PR_NGMN_to_MPA_results.json"

for path in (protein_path, docked_path, results_path):
    assert path.exists(), f"Missing input: {path}"

results = json.loads(results_path.read_text())

print("========== PR NGMN→MPA SASA ANALYSIS ==========")

bondi = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    "F": 1.47,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98
```

### Source lines 21221–21385
```python
    score = float(score)
    centroid = float(centroid)
    minimum_distance = float(minimum_distance)
    jaccard = float(jaccard)
    shielded_fraction = float(shielded_fraction)

    assert all(np.isfinite(value) for value in [
        score,
        centroid,
        minimum_distance,
        jaccard,
        shielded_fraction
    ]), f"Non-finite measurement in {filename}"

    assert centroid >= 0
    assert minimum_distance >= 0
    assert 0 <= jaccard <= 1
    assert 0 <= shielded_fraction <= 1

    rows.append({
        "Sequence": f"{first}→{second}",
        "First ligand": first,
        "Incoming ligand": second,
        "Best local Vina score (kcal/mol)": score,
        "Centroid distance (Å)": centroid,
        "Minimum interligand distance (Å)": minimum_distance,
        "Shared PR residues (n)": len(shared),
        "Shared PR residue identities": ", ".join(shared),
        "PR contact Jaccard": jaccard,
        "PR-shielded fraction (%)": shielded_fraction * 100
    })

    print(f"PASS: {first}→{second}")

df = pd.DataFrame(rows)

assert len(df) == 12
assert df["Sequence"].nunique() == 12

# Save machine-readable comparison.
csv_path = work / "PR_Dock_EM_Dock_12_sequence_comparison.csv"
xlsx_path = work / "PR_Dock_EM_Dock_12_sequence_comparison.xlsx"
json_path = work / "PR_Dock_EM_Dock_12_sequence_summary.json"

df.to_csv(csv_path, index=False)

with pd.ExcelWriter(
    xlsx_path,
    engine="openpyxl"
) as writer:
    df.to_excel(
        writer,
        sheet_name="PR 12 Sequences",
        index=False
    )

summary = {
    "receptor": "PR",
    "number_of_sequences": len(df),
    "sequences": rows,
    "interpretive_limitations": [
        "Local Vina scores are not experimental binding free energies.",
        "Original rigid-receptor scores and Dock–EM–Dock scores "
        "were generated using different protocols.",
        "PR shielding excludes the fixed first ligand.",
        "Contact overlap does not establish cooperativity.",
        "All spatial measurements describe best-ranked docked poses."
    ]
}

json_path.write_text(
    json.dumps(summary, indent=2) + "\n"
)

print("\n========== COMPLETE PR COMPARISON ==========")

display_df = df.copy()

for column in [
    "Best local Vina score (kcal/mol)",
    "Centroid distance (Å)",
    "Minimum interligand distance (Å)",
    "PR contact Jaccard",
    "PR-shielded fraction (%)"
]:
    display_df[column] = display_df[column].round(4)

print(display_df.to_string(index=False))

print("\n========== SUMMARY STATISTICS ==========")
print(
    "Local Vina score range (kcal/mol):",
    round(df["Best local Vina score (kcal/mol)"].min(), 4),
    "to",
    round(df["Best local Vina score (kcal/mol)"].max(), 4)
)
print(
    "Centroid distance range (Å):",
    round(df["Centroid distance (Å)"].min(), 4),
    "to",
    round(df["Centroid distance (Å)"].max(), 4)
)
print(
    "Minimum interligand distance range (Å):",
    round(df["Minimum interligand distance (Å)"].min(), 4),
    "to",
    round(df["Minimum interligand distance (Å)"].max(), 4)
)
print(
    "PR-shielded fraction range (%):",
    round(df["PR-shielded fraction (%)"].min(), 4),
    "to",
    round(df["PR-shielded fraction (%)"].max(), 4)
)

print("\n========== OUTPUT FILES ==========")
print("CSV:", csv_path)
print("Excel:", xlsx_path)
print("JSON:", json_path)
print("SUCCESS: All 12 PR sequences consolidated.")

from pathlib import Path
import json
import re

work = Path("/content/PR_Dock_EM_Dock")
ligands = ["EE", "LNG", "NGMN", "MPA"]

sequences = [
    (first, second)
    for first in ligands
    for second in ligands
    if first != second
]

keywords = (
    "score", "affinity", "energy",
    "centroid", "distance", "jaccard",
    "contact", "shared", "sasa",
    "shield", "fraction"
)

def summarize(value, max_length=180):
    if isinstance(value, (dict, list)):
        rendered = json.dumps(value, ensure_ascii=False)
    else:
        rendered = repr(value)

    if len(rendered) > max_length:
        return rendered[:max_length] + "..."

    return rendered

print("========== PR JSON SCHEMA DIAGNOSTIC ==========")

missing_files = []
schema_groups = {}

for first, second in sequences:
    name = f"PR_{first}_to_{second}_results.json"
    path = work / name

    print("\n" + "=" * 70)
    print(f"SEQUENCE: {first}→{second}")
    print("FILE:", path)
```

### Source lines 21590–21791
```python
        ])

        if isolated is not None and bound is not None:
            isolated = float(isolated)
            bound = float(bound)
            if isolated > 0:
                return (isolated - bound) / isolated

    return None

rows = []
errors = []

print("========== PR 12-SEQUENCE CONSOLIDATION ==========")

for a, b in sequences:
    seq = f"{a}→{b}"
    path = work / f"PR_{a}_to_{b}_results.json"

    try:
        d = json.loads(path.read_text())

        score = numeric(score_for(d), "Vina score", seq)
        centroid = numeric(
            centroid_for(d), "centroid distance", seq
        )
        minimum = numeric(
            minimum_distance_for(d, a, b),
            "minimum interligand distance", seq
        )
        jaccard = numeric(
            jaccard_for(d), "PR contact Jaccard", seq
        )
        fraction = numeric(
            shielding_for(d, a, b),
            "PR-shielded SASA fraction", seq
        )

        shared = shared_for(d)

        if not isinstance(shared, list):
            raise ValueError(
                f"{seq}: shared residue identities unavailable"
            )

        if not 0 <= jaccard <= 1:
            raise ValueError(f"{seq}: invalid Jaccard")

        if not 0 <= fraction <= 1:
            raise ValueError(f"{seq}: invalid SASA fraction")

        if minimum < 0 or centroid < 0:
            raise ValueError(f"{seq}: invalid distance")

        rows.append({
            "Sequence": seq,
            "Fixed ligand": a,
            "Incoming ligand": b,
            "Best local Vina score (kcal/mol)": score,
            "Centroid distance (Å)": centroid,
            "Minimum interligand distance (Å)": minimum,
            "Shared PR residues (n)": len(shared),
            "Shared PR residue identities": ", ".join(shared),
            "PR contact Jaccard": jaccard,
            "PR-shielded fraction (%)": fraction * 100
        })

        print(f"PASS: {seq}")

    except Exception as exc:
        errors.append(f"{seq}: {exc}")
        print(f"NEEDS REVIEW: {seq}: {exc}")

print("\n========== VALIDATION SUMMARY ==========")
print("Successfully normalized:", len(rows), "/ 12")
print("Issues:", len(errors))

if errors:
    print("\nIssues requiring correction:")
    for error in errors:
        print(" -", error)
    print(
        "\nNo consolidated files written. "
        "Original JSON files remain unchanged."
    )

else:
    df = pd.DataFrame(rows)

    assert len(df) == 12
    assert df["Sequence"].nunique() == 12

    csv_path = work / "PR_Dock_EM_Dock_12_sequence_comparison.csv"
    xlsx_path = work / "PR_Dock_EM_Dock_12_sequence_comparison.xlsx"
    json_path = work / "PR_Dock_EM_Dock_12_sequence_summary.json"

    df.to_csv(csv_path, index=False)

    with pd.ExcelWriter(
        xlsx_path, engine="openpyxl"
    ) as writer:
        df.to_excel(
            writer,
            sheet_name="PR 12 Sequences",
            index=False
        )

    json_path.write_text(
        json.dumps({
            "receptor": "PR",
            "sequence_count": 12,
            "results": rows,
            "notes": [
                "Original result JSON files were not modified.",
                "All scores are local Vina docking scores.",
                "SASA shielding excludes the fixed ligand.",
                "Shared receptor contacts do not prove cooperativity.",
                "Original rigid-receptor and Dock–EM–Dock "
                "scores are not directly comparable as "
                "binding-affinity changes."
            ]
        }, indent=2) + "\n"
    )

    print("\n========== COMPLETE PR COMPARISON ==========")

    with pd.option_context(
        "display.max_rows", 20,
        "display.max_columns", None,
        "display.width", 220
    ):
        print(df.round(4).to_string(index=False))

    print("\n========== SUMMARY STATISTICS ==========")

    for column in [
        "Best local Vina score (kcal/mol)",
        "Centroid distance (Å)",
        "Minimum interligand distance (Å)",
        "PR contact Jaccard",
        "PR-shielded fraction (%)"
    ]:
        print(
            column,
            "| min:", round(df[column].min(), 4),
            "| max:", round(df[column].max(), 4)
        )

    print("\n========== OUTPUT FILES ==========")
    print("CSV:", csv_path)
    print("Excel:", xlsx_path)
    print("JSON:", json_path)
    print("SUCCESS: All 12 PR sequences consolidated.")

from pathlib import Path
import json
import math
import re
import pandas as pd

work = Path("/content/PR_Dock_EM_Dock")
target = work / "PR_LNG_to_MPA_results.json"

assert target.exists(), f"Missing: {target}"
data = json.loads(target.read_text())

print("========== LNG→MPA SCORE VERIFICATION ==========")

# Inspect possible saved docking artifacts.
candidates = sorted(
    p for p in work.iterdir()
    if p.is_file()
    and "LNG" in p.name.upper()
    and "MPA" in p.name.upper()
    and p.suffix.lower() in (".pdbqt", ".txt", ".log")
)

verified_scores = []

for path in candidates:
    print("Checking:", path.name)

    content = path.read_text(errors="replace")

    # Vina PDBQT pose annotation.
    for match in re.finditer(
        r"REMARK\s+VINA\s+RESULT:\s*"
        r"(-?\d+(?:\.\d+)?)",
        content
    ):
        verified_scores.append(
            (float(match.group(1)), path.name, "PDBQT")
        )
        break

    # Explicitly labeled best-score log entries.
    for match in re.finditer(
        r"(?:best\s+(?:local\s+)?(?:vina\s+)?"
        r"(?:score|affinity))"
        r"[^\n:]*:\s*(-?\d+(?:\.\d+)?)",
        content,
        flags=re.IGNORECASE
```

### Source lines 22129–22385
```python
    try:
        assert path.exists(), f"Missing file: {path}"

        data = json.loads(path.read_text())

        score = required_number(
            choose(data, [
                ("best_local_vina_score_kcal_mol",),
                ("docking_score_kcal_mol",),
                ("docking", "best_affinity_kcal_mol")
            ]),
            "Vina score",
            sequence
        )

        centroid = required_number(
            choose(data, [
                ("centroid_distance_A",),
                ("interligand_centroid_distance_A",),
                ("best_pose_spatial_analysis",
                 "centroid_distance_angstrom"),
                ("spatial_analysis",
                 "centroid_distance_angstrom")
            ]),
            "centroid distance",
            sequence
        )

        minimum = required_number(
            choose(data, [
                ("minimum_interligand_heavy_distance_A",),
                ("minimum_interligand_distance_A",),
                ("best_pose_spatial_analysis",
                 f"minimum_{a}_{b}_distance_angstrom"),
                ("spatial_analysis",
                 "minimum_interligand_distance_angstrom")
            ]),
            "minimum interligand distance",
            sequence
        )

        jaccard = required_number(
            choose(data, [
                ("PR_contact_jaccard",),
                ("PR_contact_jaccard_4_A",),
                ("contact_jaccard_4A",),
                ("best_pose_spatial_analysis",
                 "contact_residue_jaccard_similarity"),
                ("spatial_analysis",
                 "contact_jaccard_similarity")
            ]),
            "PR contact Jaccard",
            sequence
        )

        shielding = required_number(
            extract_sasa(data, b),
            "PR-shielded fraction",
            sequence
        )

        shared = extract_shared(data)

        if not isinstance(shared, list):
            raise ValueError(
                "Shared PR residue identities unavailable"
            )

        if not (0 <= shielding <= 1):
            raise ValueError("Invalid shielding fraction")

        if not (0 <= jaccard <= 1):
            raise ValueError("Invalid Jaccard index")

        if centroid < 0 or minimum < 0:
            raise ValueError("Invalid distance")

        rows.append({
            "Sequence": sequence,
            "Fixed ligand": a,
            "Incoming ligand": b,
            "Original rigid score (kcal/mol)": original_scores[(a, b)],
            "Dock-EM-Dock local Vina score (kcal/mol)": score,
            "Centroid distance (Å)": centroid,
            "Minimum interligand distance (Å)": minimum,
            "Shared PR residues (n)": len(shared),
            "Shared PR residue identities": ", ".join(shared),
            "PR contact Jaccard": jaccard,
            "PR-shielded fraction (%)": shielding * 100
        })

        print(f"PASS: {sequence}")

    except Exception as exc:
        errors.append(f"{sequence}: {exc}")
        print(f"NEEDS REVIEW: {sequence}: {exc}")

print("\n========== VALIDATION ==========")
print("Normalized:", len(rows), "/ 12")
print("Errors:", len(errors))

if errors:
    for error in errors:
        print(" -", error)

    print(
        "\nNo consolidated files written. "
        "Original JSON files remain unchanged."
    )

else:
    df = pd.DataFrame(rows)

    assert len(df) == 12
    assert df["Sequence"].nunique() == 12

    csv_path = work / "PR_Dock_EM_Dock_12_sequence_comparison.csv"
    xlsx_path = work / "PR_Dock_EM_Dock_12_sequence_comparison.xlsx"
    json_path = work / "PR_Dock_EM_Dock_12_sequence_summary.json"

    df.to_csv(csv_path, index=False)

    with pd.ExcelWriter(
        xlsx_path,
        engine="openpyxl"
    ) as writer:
        df.to_excel(
            writer,
            sheet_name="PR 12 Sequences",
            index=False
        )

    json_path.write_text(
        json.dumps({
            "receptor": "PR",
            "sequence_count": 12,
            "sequences": rows,
            "interpretation_limits": [
                "The two score columns use different protocols.",
                "Score differences are not changes in binding free energy.",
                "SASA shielding excludes the fixed ligand.",
                "PR contact overlap does not prove cooperativity.",
                "Spatial measurements use best-ranked docked poses."
            ]
        }, indent=2) + "\n"
    )

    print("\n========== COMPLETE PR COMPARISON ==========")

    with pd.option_context(
        "display.max_rows", 20,
        "display.max_columns", None,
        "display.width", 240
    ):
        print(df.round(4).to_string(index=False))

    print("\n========== SUMMARY STATISTICS ==========")

    for column in [
        "Dock-EM-Dock local Vina score (kcal/mol)",
        "Centroid distance (Å)",
        "Minimum interligand distance (Å)",
        "PR contact Jaccard",
        "PR-shielded fraction (%)"
    ]:
        print(
            column,
            "| minimum:", round(df[column].min(), 4),
            "| maximum:", round(df[column].max(), 4)
        )

    print("\n========== SAVED FILES ==========")
    print("CSV:", csv_path)
    print("Excel:", xlsx_path)
    print("JSON:", json_path)
    print("SUCCESS: All 12 PR sequences consolidated.")

from pathlib import Path
import json
import re

work = Path("/content/PR_Dock_EM_Dock")
assert work.exists(), f"Directory not found: {work}"

ligands = ["EE", "LNG", "NGMN", "MPA"]

keywords = [
    "forcefield", "force_field", "amber", "openff",
    "sage", "solvent", "implicit", "explicit",
    "vacuum", "nonbonded", "cutoff", "pme",
    "restraint", "minimiz", "tolerance",
    "convergence", "iteration", "platform",
    "potential_energy", "rmsd", "centroid_shift",
    "am1bcc", "charge"
]

def relevant(text):
    lower = text.lower()
    return any(word in lower for word in keywords)

def flatten(obj, prefix=""):
    if isinstance(obj, dict):
        for key, value in obj.items():
            path = f"{prefix}.{key}" if prefix else key
            yield from flatten(value, path)

    elif isinstance(obj, list):
        if len(obj) <= 8:
            for i, value in enumerate(obj):
                yield from flatten(value, f"{prefix}[{i}]")
        else:
            yield prefix, f"[list of {len(obj)} items]"

    else:
        yield prefix, obj

print("=" * 75)
print("PR DOCK–EM–DOCK: MINIMIZATION PROTOCOL AUDIT")
print("=" * 75)

for ligand in ligands:
    print("\n" + "=" * 75)
    print(f"FIRST LIGAND: {ligand}")
    print("=" * 75)

    # Files beginning with the first-ligand prefix.
    matches = sorted(
        p for p in work.iterdir()
        if p.is_file()
        and p.name.upper().startswith(f"PR_{ligand}_")
    )

    print("\n[1] RELEVANT FILE INVENTORY")

    for p in matches:
        if any(word in p.name.lower() for word in [
            "minim", "protocol", "setup", "prepare",
            "repair", "energy", "stereo", "log"
        ]):
            print(
                f"  {p.name} "
                f"({p.stat().st_size:,} bytes)"
            )

    print("\n[2] MINIMIZED STRUCTURE FILES")

    for suffix in [
        "minimized_protein.pdb",
        "minimized_ligand.sdf",
        "minimized_complex.pdb"
    ]:
        path = work / f"PR_{ligand}_{suffix}"
        print(
            f"  {suffix}: "
            f"{'FOUND' if path.exists() else 'NOT FOUND'}"
        )

```

### Source lines 25563–25729
```python
            "Sequence": f"{first}→{second}",
            "Protein_alignment_RMSD_A": round(sup.rms, 4),
            "First_to_progesterone_centroid_A": round(
                first_centroid, 4
            ),
            "First_to_progesterone_minimum_A": round(
                first_minimum, 4
            ),
            "Second_to_progesterone_centroid_A": round(
                second_centroid, 4
            ),
            "Second_to_progesterone_minimum_A": round(
                second_minimum, 4
            )
        })

df = pd.DataFrame(rows)

print("=" * 76)
print("PR EXPERIMENTAL ORTHOSTERIC-SITE SPATIAL COMPARISON")
print("=" * 76)
print(df.to_string(index=False))

output = work / "PR_12_sequence_experimental_site_distances.csv"
df.to_csv(output, index=False)

print("\nSaved:", output)
print("No docking or minimization performed.")

from pathlib import Path
import pandas as pd

work = Path("/content/PR_Dock_EM_Dock")

spatial = pd.read_csv(
    work / "PR_12_sequence_binding_site_mapping_corrected.csv"
)

orthosteric = pd.read_csv(
    work / "PR_12_sequence_orthosteric_pocket_comparison.csv"
)

reference = pd.read_csv(
    work / "PR_12_sequence_experimental_site_distances.csv"
)

df = spatial.merge(
    orthosteric,
    on="Sequence",
    validate="one_to_one"
).merge(
    reference,
    on="Sequence",
    validate="one_to_one"
)

if len(df) != 12:
    raise RuntimeError(
        f"Expected 12 sequences; found {len(df)}"
    )

df["First_ligand_site"] = (
    "Experimental orthosteric pocket"
)

df["Second_ligand_site"] = (
    "Spatially distinct predicted secondary region"
)

df["Binding_pattern"] = (
    "Orthosteric-first / secondary-site second"
)

df["Allostery_status"] = (
    "Putative only; functional allostery not established"
)

columns = [
    "Sequence",
    "First_ligand_site",
    "Second_ligand_site",
    "Binding_pattern",
    "First_to_progesterone_centroid_A",
    "Second_to_progesterone_centroid_A",
    "Minimum_interligand_distance_A",
    "Contact_Jaccard",
    "First_orthosteric_contacts",
    "Second_orthosteric_contacts",
    "Allostery_status"
]

final = df[columns].copy()

print("=" * 80)
print("FINAL PR DOCK–EM–DOCK BINDING-SITE CLASSIFICATION")
print("=" * 80)

print(final.to_string(index=False))

csv_path = work / "PR_12_sequence_final_binding_site_classification.csv"
xlsx_path = work / "PR_12_sequence_final_binding_site_classification.xlsx"

final.to_csv(csv_path, index=False)
final.to_excel(xlsx_path, index=False)

print("\nSaved:")
print(csv_path)
print(xlsx_path)

print("\nSequences classified:", len(final))
print("No docking or minimization performed.")

from google.colab import drive
from pathlib import Path
import shutil
import zipfile
from datetime import datetime

# Connect Google Drive
drive.mount("/content/drive")

# Project folder containing all PR calculations
source = Path("/content/PR_Dock_EM_Dock")

# Permanent Google Drive backup location
backup_folder = Path(
    "/content/drive/MyDrive/Dock_EM_Dock_Research/Progesterone_Receptor"
)
backup_folder.mkdir(parents=True, exist_ok=True)

if not source.is_dir():
    raise FileNotFoundError(
        f"PR project folder not found: {source}"
    )

# Timestamp prevents overwriting earlier backups
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
zip_path = backup_folder / f"PR_COMPLETE_BACKUP_{timestamp}.zip"

print("=" * 70)
print("BACKING UP COMPLETE PR DOCK–EM–DOCK PROJECT")
print("=" * 70)

files = [p for p in source.rglob("*") if p.is_file()]

print("Source:", source)
print("Files to back up:", len(files))
print("Destination:", zip_path)

# Save every file while preserving folder structure
with zipfile.ZipFile(
    zip_path, "w", compression=zipfile.ZIP_DEFLATED,
    compresslevel=6
) as archive:
    for file in files:
        archive.write(
            file,
            arcname=file.relative_to(source.parent)
        )

# Verify archive integrity
with zipfile.ZipFile(zip_path, "r") as archive:
    bad_file = archive.testzip()
    archived_files = archive.namelist()

if bad_file is not None:
    raise RuntimeError(
```


---

# Estrogen receptor alpha (ERalpha)

# ERalpha — original Dock–EM–Dock code excerpts
**Source:** `er_ee_energyminimization.py` (13,387 lines).
**Provenance:** All code excerpts below are copied verbatim from the uploaded Colab-exported Python source, with original line ranges. Sections are excerpts, not a standalone runnable script. Original complete file is included in the ZIP.
**Scope:** The original notebook may contain debugging, corrections, alternate parameterizations, or non-executable Colab shell commands. These excerpts do not imply all code paths were executed successfully.

## Original code excerpts
### Source lines 684–958
```python
            fixer_h.positions,
            f,
            keepIds=True
        )

    print("Saved:", output_file)

except Exception as e:
    print("Hydrogen addition encountered an error.")
    print("Error type:", type(e).__name__)
    print("Details:", str(e))

from openmm.app import PDBFile, ForceField, NoCutoff, HBonds

# Load hydrogenated estrogen receptor
protein = PDBFile("/content/er_protein_h.pdb")

# Load AMBER14 force field
forcefield = ForceField("amber14-all.xml")

print("=== FORCE FIELD COMPATIBILITY CHECK ===")
print("Total atoms:", protein.topology.getNumAtoms())

try:
    system = forcefield.createSystem(
        protein.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds
    )

    print("\nSUCCESS!")
    print("AMBER14 force field is compatible.")
    print("Number of particles:", system.getNumParticles())
    print("Protein is ready for energy minimization.")

except Exception as e:
    print("\nForce-field compatibility failed.")
    print("Error type:", type(e).__name__)
    print("Details:", str(e))

from openmm.app import PDBFile, ForceField, NoCutoff, HBonds

# Load hydrogenated estrogen receptor
protein = PDBFile("/content/er_protein_h.pdb")

# Load AMBER14 force field
forcefield = ForceField("amber14-all.xml")

print("=== FORCE FIELD COMPATIBILITY CHECK ===")
print("Total atoms:", protein.topology.getNumAtoms())

try:
    system = forcefield.createSystem(
        protein.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds
    )

    print("\nSUCCESS!")
    print("AMBER14 force field is compatible.")
    print("Number of particles:", system.getNumParticles())
    print("Protein is ready for energy minimization.")

except Exception as e:
    print("\nForce-field compatibility failed.")
    print("Error type:", type(e).__name__)
    print("Details:", str(e))

import openmm as mm
from openmm import unit
from openmm.app import (
    PDBFile, ForceField, Simulation,
    NoCutoff, HBonds
)

# Load hydrogenated protein
protein = PDBFile("/content/er_protein_h.pdb")

# Set up AMBER14
forcefield = ForceField("amber14-all.xml")

system = forcefield.createSystem(
    protein.topology,
    nonbondedMethod=NoCutoff,
    constraints=HBonds
)

# Restrain heavy atoms to their initial positions
restraint = mm.CustomExternalForce(
    "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
)

restraint.addGlobalParameter(
    "k", 1000.0
)

for name in ["x0", "y0", "z0"]:
    restraint.addPerParticleParameter(name)

positions_nm = protein.positions.value_in_unit(
    unit.nanometer
)

for atom in protein.topology.atoms():
    if atom.element.symbol != "H":
        pos = positions_nm[atom.index]
        restraint.addParticle(
            atom.index,
            [pos.x, pos.y, pos.z]
        )

system.addForce(restraint)

# Initialize simulation
integrator = mm.VerletIntegrator(
    0.001 * unit.picoseconds
)

simulation = Simulation(
    protein.topology,
    system,
    integrator
)

simulation.context.setPositions(protein.positions)

# Calculate initial potential energy
initial = simulation.context.getState(
    getEnergy=True
).getPotentialEnergy()

print("Initial potential energy:", initial)

# Run restrained minimization
print("Starting restrained energy minimization...")

simulation.minimizeEnergy(
    maxIterations=1000
)

# Calculate final potential energy
state = simulation.context.getState(
    getEnergy=True,
    getPositions=True
)

final = state.getPotentialEnergy()

print("Final potential energy:", final)

# Save minimized structure
output_file = "/content/er_protein_minimized.pdb"

with open(output_file, "w") as f:
    PDBFile.writeFile(
        simulation.topology,
        state.getPositions(),
        f,
        keepIds=True
    )

print("Saved:", output_file)
print("Restrained minimization completed!")

import numpy as np
from openmm.app import PDBFile
from openmm import unit

original = PDBFile("/content/er_protein_h.pdb")
minimized = PDBFile("/content/er_protein_minimized.pdb")

# Verify that the atom order is unchanged
original_atoms = list(original.topology.atoms())
minimized_atoms = list(minimized.topology.atoms())

assert len(original_atoms) == len(minimized_atoms)

for a, b in zip(original_atoms, minimized_atoms):
    assert (
        a.name == b.name
        and a.residue.name == b.residue.name
        and a.residue.id == b.residue.id
    )

# Select backbone atoms
backbone = [
    atom.index
    for atom in original_atoms
    if atom.name in ("N", "CA", "C", "O")
]

# Coordinates in Angstroms
pos1 = np.array(
    original.positions.value_in_unit(unit.angstrom)
)[backbone]

pos2 = np.array(
    minimized.positions.value_in_unit(unit.angstrom)
)[backbone]

# RMSD without alignment: measures movement
# relative to the original coordinate frame
difference = pos2 - pos1

rmsd = np.sqrt(
    np.mean(np.sum(difference**2, axis=1))
)

print("=== BACKBONE DISPLACEMENT CHECK ===")
print("Backbone atoms:", len(backbone))
print(f"Backbone RMSD (unaligned): {rmsd:.4f} Å")
print("Verification complete!")

import os
from rdkit import Chem

print("=== AVAILABLE DOCKING FILES ===")

for filename in sorted(os.listdir("/content")):
    if filename.startswith("estrogenreceptor_ee_out_1"):
        path = os.path.join("/content", filename)
        print(f"{filename} | {os.path.getsize(path)} bytes")

print("\n=== EXISTING PREPARED FILES ===")

for filename in ["ee.sdf", "ee_prepared.sdf",
                 "er_protein_h.pdb", "er_protein_minimized.pdb"]:
    path = os.path.join("/content", filename)
    print(filename, "FOUND" if os.path.exists(path) else "NOT FOUND")

print("\n=== DOCKED COMPLEX CONTENT ===")

for filename in sorted(os.listdir("/content")):
    if filename.startswith("estrogenreceptor_ee_out_1") and filename.lower().endswith(".pdb"):
        path = os.path.join("/content", filename)
        with open(path) as f:
            lines = f.readlines()
        protein_atoms = sum(line.startswith("ATOM  ") for line in lines)
        ligand_atoms = sum(line.startswith("HETATM") for line in lines)
        print(filename)
        print("Protein ATOM records:", protein_atoms)
        print("HETATM records:", ligand_atoms)
        print("Residue names:", sorted(set(
            line[17:20].strip() for line in lines
            if line.startswith(("ATOM  ", "HETATM"))
        ))[-15:])

from rdkit import Chem
from collections import Counter
import os

complex_file = "/content/estrogenreceptor_ee_out_1._9-9.complex.pdb"
ligand_file = "/content/estrogenreceptor_ee_out_1._9-9.mol2"

print("=== DOCKED EE VERIFICATION ===")

with open(complex_file) as f:
    hetatm_lines = [
        line for line in f
        if line.startswith("HETATM")
    ]

print("HETATM count:", len(hetatm_lines))
print("Residue identities:",
      Counter(line[17:20].strip() for line in hetatm_lines))

print("\n=== DOCKED LIGAND FILE ===")
print("MOL2 exists:", os.path.exists(ligand_file))

if os.path.exists(ligand_file):
    with open(ligand_file) as f:
        mol2_text = f.read()

    print("MOL2 sections:",
          [line for line in mol2_text.splitlines()
```

### Source lines 1337–1425
```python

from openmm import app
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Molecule

# Load the prepared protein
protein = app.PDBFile("/content/er_protein_h.pdb")

# Load docked EE with its original coordinates
ee = Molecule.from_file(
    "/content/ee_docked_restored.sdf",
    allow_undefined_stereo=False
)

print("Protein atoms:", protein.topology.getNumAtoms())
print("EE atoms:", ee.n_atoms)

# Set up force field
forcefield = app.ForceField(
    "amber14-all.xml",
    "implicit/obc2.xml"
)

# Register ligand parameters
smirnoff = SMIRNOFFTemplateGenerator(
    molecules=[ee]
)

forcefield.registerTemplateGenerator(
    smirnoff.generator
)

print("Protein and EE force-field setup completed!")

# Test a more specific OpenFF import
from openff.toolkit.topology import Molecule

print("OpenFF Molecule imported successfully!")

import importlib.metadata as md

packages = [
    "torch",
    "torchvision",
    "torchaudio",
    "openff-toolkit",
    "openmm",
    "openmmforcefields"
]

for package in packages:
    try:
        print(f"{package}: {md.version(package)}")
    except md.PackageNotFoundError:
        print(f"{package}: NOT INSTALLED")

# Commented out IPython magic to ensure Python compatibility.

# %pip uninstall -y torchvision

import sys
import importlib.util
import subprocess

print("Python executable:", sys.executable)

for package in ["torch", "torchvision", "openmm", "openff"]:
    try:
        spec = importlib.util.find_spec(package)
        print(f"{package}: {'FOUND' if spec else 'NOT FOUND'}")
    except Exception as e:
        print(f"{package}: ERROR - {e}")

print("\nChecking pip environment:")
subprocess.run([
    sys.executable, "-m", "pip",
    "show", "torchvision"
])

import importlib.util
import sys
import subprocess

spec = importlib.util.find_spec("torchvision")

print("TorchVision location:", spec.origin)
print("Python:", sys.executable)

print("\nConda installation:")
```

### Source lines 1469–1768
```python
    allow_undefined_stereo=False
)

print("EE ligand loaded successfully!")
print("EE atoms:", ee.n_atoms)
print("EE conformers:", ee.n_conformers)

print("\nSTEP 1 COMPLETE!")

from openmm import app
from openmm import unit

# Load prepared ER protein
protein = app.PDBFile("/content/er_protein_h.pdb")

# Select protein force field
forcefield = app.ForceField(
    "amber14-all.xml",
    "amber14/tip3p.xml"
)

print("Checking ER protein force-field compatibility...")

try:
    protein_system = forcefield.createSystem(
        protein.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None
    )

    print("SUCCESS: ER protein parameterized!")
    print("System particles:",
          protein_system.getNumParticles())

except Exception as e:
    print("PROTEIN PARAMETERIZATION FAILED")
    print(type(e).__name__, str(e))

from openff.toolkit import Molecule, ForceField

# Load restored docked EE
ee = Molecule.from_file(
    "/content/ee_docked_restored.sdf",
    file_format="SDF",
    allow_undefined_stereo=False
)

print("EE loaded:", ee.n_atoms, "atoms")
print("EE conformers:", ee.n_conformers)

# Load small-molecule force field
ligand_forcefield = ForceField(
    "openff_unconstrained-2.2.1.offxml"
)

print("SUCCESS: OpenFF ligand force field loaded!")

# Test ligand parameterization
try:
    ee_system = ligand_forcefield.create_openmm_system(
        ee.to_topology()
    )

    print("SUCCESS: EE parameterized!")
    print("EE system particles:",
          ee_system.getNumParticles())

except Exception as e:
    print("EE PARAMETERIZATION FAILED")
    print(type(e).__name__, str(e))

from openmm import app
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Molecule
import numpy as np

# Load prepared receptor and docked ligand
protein = app.PDBFile("/content/er_protein_h.pdb")
ee = Molecule.from_file(
    "/content/ee_docked_restored.sdf",
    allow_undefined_stereo=False
)

# Build combined topology
modeller = app.Modeller(
    protein.topology,
    protein.positions
)

# Add EE at its docked coordinates
ee_topology = ee.to_topology().to_openmm()
ee_positions = ee.conformers[0].to(
    "nanometer"
).magnitude

from openmm import unit

modeller.add(
    ee_topology,
    ee_positions * unit.nanometer
)

print("Combined topology created!")
print("Total atoms:", modeller.topology.getNumAtoms())

# Protein force field
forcefield = app.ForceField(
    "amber14-all.xml",
    "amber14/tip3p.xml"
)

# Ligand force-field generator
smirnoff = SMIRNOFFTemplateGenerator(
    molecules=[ee],
    forcefield="openff-2.2.1"
)

forcefield.registerTemplateGenerator(
    smirnoff.generator
)

# Parameterize combined complex
complex_system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=None
)

print("SUCCESS: ER–EE complex parameterized!")
print("System particles:",
      complex_system.getNumParticles())

import openmm as mm
from openmm import app, unit
import numpy as np

print("=== ER–EE COMPLEX ENERGY MINIMIZATION ===")

# Restrain protein heavy atoms to preserve the
# overall receptor structure while allowing relaxation
restraint = mm.CustomExternalForce(
    "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
)
restraint.addGlobalParameter("k", 1000.0)
for parameter in ["x0", "y0", "z0"]:
    restraint.addPerParticleParameter(parameter)

positions_nm = modeller.positions.value_in_unit(unit.nanometer)

for atom in modeller.topology.atoms():
    if atom.index < 3819 and atom.element.symbol != "H":
        pos = positions_nm[atom.index]
        restraint.addParticle(atom.index, [pos.x, pos.y, pos.z])

complex_system.addForce(restraint)

# Initialize minimization
integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)

simulation = app.Simulation(
    modeller.topology,
    complex_system,
    integrator,
    mm.Platform.getPlatformByName("CPU")
)

simulation.context.setPositions(modeller.positions)

initial_energy = simulation.context.getState(
    getEnergy=True
).getPotentialEnergy()

print("Initial energy:", initial_energy)
print("Running energy minimization...")

simulation.minimizeEnergy(
    tolerance=10 * unit.kilojoule_per_mole / unit.nanometer,
    maxIterations=1000
)

state = simulation.context.getState(
    getEnergy=True,
    getPositions=True
)

final_energy = state.getPotentialEnergy()

print("Final energy:", final_energy)
print("Energy minimization finished!")

# Save minimized protein–ligand complex
output_file = "/content/ER_EE_minimized_complex.pdb"

with open(output_file, "w") as f:
    app.PDBFile.writeFile(
        modeller.topology,
        state.getPositions(),
        f,
        keepIds=True
    )

print("Saved:", output_file)

import pandas as pd
from google.colab import files

# Save energy results
results = pd.DataFrame({
    "Complex": ["ER-EE"],
    "Initial_Energy_kJ_mol": [-8721.704129206226],
    "Final_Energy_kJ_mol": [-21652.868916890242],
    "Energy_Change_kJ_mol": [-12931.164787684016],
    "Max_Iterations": [1000]
})

csv_path = "/content/ER_EE_minimization_results.csv"
results.to_csv(csv_path, index=False)

# Download both files
files.download("/content/ER_EE_minimized_complex.pdb")
files.download(csv_path)

print("ER–EE results saved and downloaded!")

import numpy as np
from scipy.spatial.distance import cdist
from openmm import app, unit

# Original coordinates
original = np.array(
    modeller.positions.value_in_unit(unit.angstrom)
)

# Minimized coordinates
minimized_pdb = app.PDBFile(
    "/content/ER_EE_minimized_complex.pdb"
)

minimized = np.array(
    minimized_pdb.positions.value_in_unit(unit.angstrom)
)

# Protein has 3819 atoms; EE has 46
protein_count = 3819

# Ligand coordinates before and after
ee_before = original[protein_count:]
ee_after = minimized[protein_count:]

# Protein coordinates after minimization
protein_after = minimized[:protein_count]

# Measure ligand movement
movement = np.linalg.norm(
    ee_after - ee_before, axis=1
)

# Measure minimum protein-ligand distance
distances = cdist(ee_after, protein_after)

print("=== ER–EE VALIDATION ===")
print("Mean EE atom displacement:",
      round(movement.mean(), 3), "Å")
print("Maximum EE atom displacement:",
      round(movement.max(), 3), "Å")
print("Minimum ER–EE atom distance:",
      round(distances.min(), 3), "Å")
print("Validation calculations complete!")

import numpy as np
from scipy.spatial.distance import cdist
from openmm import app, unit

pdb = app.PDBFile(
    "/content/ER_EE_minimized_complex.pdb"
)

coords = np.array(
    pdb.positions.value_in_unit(unit.angstrom)
)

atoms = list(pdb.topology.atoms())

protein_count = 3819

# Select only non-hydrogen atoms
protein_heavy = [
    a.index for a in atoms[:protein_count]
    if a.element.symbol != "H"
]

ligand_heavy = [
    a.index for a in atoms[protein_count:]
    if a.element.symbol != "H"
]

# Heavy-atom distances
dist = cdist(
    coords[ligand_heavy],
    coords[protein_heavy]
```

### Source lines 1769–1843
```python
)

minimum = dist.min()

# Count close heavy-atom pairs
clashes_2A = np.sum(dist < 2.0)
contacts_3A = np.sum(dist < 3.0)
contacts_4A = np.sum(dist < 4.0)

print("=== ER–EE HEAVY-ATOM VALIDATION ===")
print("Minimum heavy-atom distance:",
      round(minimum, 3), "Å")
print("Heavy-atom pairs below 2.0 Å:",
      clashes_2A)
print("Heavy-atom pairs below 3.0 Å:",
      contacts_3A)
print("Heavy-atom pairs below 4.0 Å:",
      contacts_4A)

import numpy as np
from scipy.spatial.distance import cdist
from openmm import app, unit

pdb = app.PDBFile(
    "/content/ER_EE_minimized_complex.pdb"
)

atoms = list(pdb.topology.atoms())
coords = np.array(
    pdb.positions.value_in_unit(unit.angstrom)
)

protein_count = 3819

protein_heavy = [
    a for a in atoms[:protein_count]
    if a.element.symbol != "H"
]

ligand_heavy = [
    a for a in atoms[protein_count:]
    if a.element.symbol != "H"
]

distances = cdist(
    coords[[a.index for a in ligand_heavy]],
    coords[[a.index for a in protein_heavy]]
)

print("=== CLOSE ER–EE CONTACTS (< 3.0 Å) ===")

for i, j in np.argwhere(distances < 3.0):
    ligand_atom = ligand_heavy[i]
    protein_atom = protein_heavy[j]
    residue = protein_atom.residue

    print(
        f"EE {ligand_atom.name} ({ligand_atom.element.symbol}) "
        f"— ER {residue.name}{residue.id} "
        f"{protein_atom.name} ({protein_atom.element.symbol}): "
        f"{distances[i, j]:.3f} Å"
    )

import numpy as np
from openmm import app, unit

pdb = app.PDBFile("/content/ER_EE_minimized_complex.pdb")
atoms = list(pdb.topology.atoms())
coords = np.array(
    pdb.positions.value_in_unit(unit.angstrom)
)

protein_count = 3819

# Find the two close-contact pairs
```

### Source lines 2502–2590
```python
)

print("Original protein heavy atoms:",
      protein_heavy_original)
print("Original EE heavy atoms:",
      ee_heavy_original)
print("Expected total heavy atoms:",
      expected_heavy)
print("PDBQT heavy-atom records:",
      heavy_total)

if heavy_total == expected_heavy:
    print("PASS: Heavy-atom counts match!")
else:
    print("WARNING: Heavy-atom counts differ.")
    print("Do not dock until this is resolved.")

import subprocess
import os

print("=== DOCKING NGMN INTO MINIMIZED ER–EE COMPLEX ===")

command = [
    "vina",
    "--receptor", "/content/ER_EE_receptor.pdbqt",
    "--ligand", "/content/NGMN_ligand.pdbqt",

    "--center_x", "9.258",
    "--center_y", "46.190",
    "--center_z", "131.409",

    "--size_x", "21.093",
    "--size_y", "25.252",
    "--size_z", "21.259",

    "--exhaustiveness", "16",
    "--num_modes", "10",

    "--out", "/content/NGMN_docked_ER_EE.pdbqt"
]

result = subprocess.run(
    command,
    capture_output=True,
    text=True
)

print(result.stdout)

if result.stderr:
    print("Additional messages:")
    print(result.stderr)

print("Return code:", result.returncode)

output_file = "/content/NGMN_docked_ER_EE.pdbqt"

if result.returncode == 0 and os.path.exists(output_file):
    print("Docking completed.")
    print("Saved:", output_file)
else:
    print("Docking did not complete successfully.")
    print("Please share the output before proceeding.")

import numpy as np
from scipy.spatial.distance import cdist

print("=== NGMN DOCKING POSE VALIDATION ===")

def read_atoms(path, first_model_only=False):
    atoms = []
    with open(path) as f:
        for line in f:
            if first_model_only and line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM  ", "HETATM")):
                atom_type = line.split()[-1]
                if atom_type not in {"H", "HD"}:
                    atoms.append([
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54])
                    ])
    return np.array(atoms)

ngmn = read_atoms(
    "/content/NGMN_docked_ER_EE.pdbqt",
    first_model_only=True
)
```

### Source lines 3544–3632
```python
print("PDBQT heavy atoms:", len(heavy_lines))
print("Output file:", output_path)
print("File size:", output_path.stat().st_size, "bytes")

assert len(heavy_lines) == mpa.GetNumHeavyAtoms(), (
    "Heavy-atom count mismatch"
)

print("\nSUCCESS: MPA is ready for docking!")

import subprocess
from pathlib import Path

print("=== DOCKING MPA INTO MINIMIZED ER–EE COMPLEX ===")

receptor = Path("/content/ER_EE_receptor.pdbqt")
ligand = Path("/content/MPA_ligand.pdbqt")
output = Path("/content/MPA_docked_ER_EE.pdbqt")

assert receptor.exists(), "ER–EE receptor file missing"
assert ligand.exists(), "MPA ligand file missing"

command = [
    "vina",
    "--receptor", str(receptor),
    "--ligand", str(ligand),

    "--center_x", "9.258",
    "--center_y", "46.190",
    "--center_z", "131.409",

    "--size_x", "21.093",
    "--size_y", "25.252",
    "--size_z", "21.259",

    "--exhaustiveness", "16",
    "--num_modes", "10",

    "--out", str(output)
]

result = subprocess.run(
    command,
    capture_output=True,
    text=True
)

print(result.stdout)

if result.stderr:
    print("\nAdditional messages:")
    print(result.stderr)

print("\nReturn code:", result.returncode)

if result.returncode == 0 and output.exists():
    print("\nSUCCESS: MPA docking completed!")
    print("Saved:", output)
    print("File size:", output.stat().st_size, "bytes")
else:
    print("\nDocking failed. Review the messages above.")

import numpy as np
from scipy.spatial.distance import cdist

print("=== VALIDATING MPA DOCKING POSES ===")

def read_heavy_coords(path):
    coords = []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                if line.split()[-1] not in {"H", "HD"}:
                    coords.append([
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54])
                    ])
    return np.array(coords)

protein_xyz = read_heavy_coords(
    "/content/ER_protein.pdbqt"
)
ee_xyz = read_heavy_coords(
    "/content/EE_minimized.pdbqt"
)

# Read saved docking models
models = []
```

### Source lines 3765–3853
```python
from meeko import MoleculePreparation, PDBQTWriterLegacy

mol = Chem.SDMolSupplier(ligand_file, removeHs=False)[0]
assert mol is not None, "Unable to read LNG SDF."
assert mol.GetNumConformers() > 0, "LNG SDF has no 3D coordinates."

mol = Chem.AddHs(mol, addCoords=True)
Chem.MolToMolFile(mol, "/content/LNG_prepared.sdf")

prep = MoleculePreparation()
setup = prep.prepare(mol)[0]
pdbqt, ok, err = PDBQTWriterLegacy.write_string(setup)
assert ok, err

with open("/content/LNG_ligand.pdbqt", "w") as f:
    f.write(pdbqt)

# Locate Vina executable
vina = shutil.which("vina")
assert vina, "Vina executable not found. Run this in your previous docking runtime."

output = "/content/LNG_docked_ER_EE.pdbqt"
command = [
    vina,
    "--receptor", receptor,
    "--ligand", "/content/LNG_ligand.pdbqt",
    "--center_x", "9.258",
    "--center_y", "46.190",
    "--center_z", "131.409",
    "--size_x", "21.093",
    "--size_y", "25.252",
    "--size_z", "21.259",
    "--exhaustiveness", "16",
    "--num_modes", "10",
    "--out", output
]

result = subprocess.run(command, capture_output=True, text=True)
print(result.stdout)
if result.returncode != 0:
    print(result.stderr)
    raise RuntimeError("LNG docking failed.")

print("Saved LNG docking poses:", output)

# Extract best docking score
import re
scores = re.findall(
    r"^\s*REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)",
    open(output).read(),
    flags=re.MULTILINE
)
if scores:
    print("Best saved LNG docking score:", scores[0], "kcal/mol")
else:
    print("Check the Vina score table above.")

import os, glob
from collections import Counter

print("=== ER–NGMN COMPLEX CHECK ===")

pdb_files = sorted(
    glob.glob("/content/*.pdb"),
    key=os.path.getmtime,
    reverse=True
)

for path in pdb_files:
    residues = Counter()
    atoms = 0

    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                atoms += 1
                key = (line[17:20].strip(), line[21:22].strip())
                residues[key] += 1

    print(f"\nFile: {os.path.basename(path)}")
    print(f"Total atoms: {atoms}")
    print("Nonstandard residues:")
    print({
        str(k): v for k, v in residues.items()
        if k[0] not in {
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU",
            "GLY","HIS","ILE","LEU","LYS","MET","PHE",
            "PRO","SER","THR","TRP","TYR","VAL","HOH"
        }
```

### Source lines 4006–4305
```python
    if trial:
        moving = Rotation.random(random_state=rng).apply(
            moving - moving.mean(axis=0)
        ) + dock_xyz.mean(axis=0)

    for iteration in range(30):
        cost = cdist(moving, dock_xyz)
        for i, e1 in enumerate(ref_elements):
            for j, e2 in enumerate(dock_elements):
                if e1 != e2:
                    cost[i, j] = 1e6

        rows, cols = linear_sum_assignment(cost)
        A = ref_xyz[rows]
        B = dock_xyz[cols]

        Ac = A.mean(axis=0)
        Bc = B.mean(axis=0)
        U, S, Vt = np.linalg.svd((A-Ac).T @ (B-Bc))
        correction = np.eye(3)
        correction[-1,-1] = np.linalg.det(U @ Vt)
        R = U @ correction @ Vt
        moving = (ref_xyz - Ac) @ R + Bc

    rmsd = np.sqrt(np.mean(
        np.sum((moving[rows] - dock_xyz[cols])**2, axis=1)
    ))

    if best is None or rmsd < best[0]:
        best = (rmsd, rows.copy(), cols.copy())

rmsd, rows, cols = best
print(f"Best heavy-atom matching RMSD: {rmsd:.3f} Å")

if rmsd > 1.0:
    raise ValueError(
        "Atom mapping is uncertain. Do not use this pose for minimization."
    )

# Assign actual docked coordinates to matched heavy atoms
conf = mol.GetConformer()
for i, j in zip(rows, cols):
    atom_idx = heavy_indices[i]
    conf.SetAtomPosition(atom_idx, Point3D(*dock_xyz[j]))

# Regenerate hydrogen positions
mol = Chem.AddHs(Chem.RemoveHs(mol), addCoords=True)
Chem.MolToMolFile(mol, output_path)

print("Saved:", output_path)
print("Atoms:", mol.GetNumAtoms())
print("Heavy atoms:", mol.GetNumHeavyAtoms())
print("Restoration complete.")

import openmm as mm
from openmm import app, unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Molecule
import numpy as np

print("=== ER–NGMN ENERGY MINIMIZATION ===")

# Same prepared receptor used for ER–EE
protein = app.PDBFile("/content/er_protein_h.pdb")

# NGMN restored at its CB-Dock3 docking coordinates
ngmn = Molecule.from_file(
    "/content/NGMN_CB_Dock3_restored.sdf",
    allow_undefined_stereo=False
)

# Build combined complex
modeller = app.Modeller(protein.topology, protein.positions)

ngmn_topology = ngmn.to_topology().to_openmm()
ngmn_positions = (
    ngmn.conformers[0].to("nanometer").magnitude
    * unit.nanometer
)

modeller.add(ngmn_topology, ngmn_positions)

# Same protein and ligand force fields as ER–EE
forcefield = app.ForceField(
    "amber14-all.xml",
    "amber14/tip3p.xml"
)

smirnoff = SMIRNOFFTemplateGenerator(
    molecules=[ngmn],
    forcefield="openff-2.2.1"
)
forcefield.registerTemplateGenerator(smirnoff.generator)

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=None
)

# Same positional restraints as ER–EE
restraint = mm.CustomExternalForce(
    "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
)
restraint.addGlobalParameter("k", 1000.0)

for name in ["x0", "y0", "z0"]:
    restraint.addPerParticleParameter(name)

positions_nm = modeller.positions.value_in_unit(unit.nanometer)
protein_atom_count = protein.topology.getNumAtoms()

for atom in modeller.topology.atoms():
    if atom.index < protein_atom_count and atom.element.symbol != "H":
        p = positions_nm[atom.index]
        restraint.addParticle(atom.index, [p.x, p.y, p.z])

system.addForce(restraint)

integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)

simulation = app.Simulation(
    modeller.topology,
    system,
    integrator,
    mm.Platform.getPlatformByName("CPU")
)

simulation.context.setPositions(modeller.positions)

initial = simulation.context.getState(
    getEnergy=True
).getPotentialEnergy()

print("Initial energy:", initial)
print("Minimizing...")

simulation.minimizeEnergy(
    tolerance=10 * unit.kilojoule_per_mole / unit.nanometer,
    maxIterations=1000
)

state = simulation.context.getState(
    getEnergy=True,
    getPositions=True
)

final = state.getPotentialEnergy()

print("Final energy:", final)

output = "/content/ER_NGMN_minimized_complex.pdb"

with open(output, "w") as f:
    app.PDBFile.writeFile(
        modeller.topology,
        state.getPositions(),
        f,
        keepIds=True
    )

print("Saved:", output)
print("Energy minimization complete.")

ngmn = Molecule.from_file(
    "/content/NGMN_CB_Dock3_restored.sdf",
    allow_undefined_stereo=False
)

import openmm as mm
from openmm import app, unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Molecule
import numpy as np

print("=== ER–NGMN ENERGY MINIMIZATION ===")

# Same prepared receptor used for ER–EE
protein = app.PDBFile("/content/er_protein_h.pdb")

# NGMN restored at its CB-Dock3 docking coordinates
ngmn = Molecule.from_file(
    "/content/NGMN_CB_Dock3_restored.sdf",
    allow_undefined_stereo=True
)

# Build combined complex
modeller = app.Modeller(protein.topology, protein.positions)

ngmn_topology = ngmn.to_topology().to_openmm()
ngmn_positions = (
    ngmn.conformers[0].to("nanometer").magnitude
    * unit.nanometer
)

modeller.add(ngmn_topology, ngmn_positions)

# Same protein and ligand force fields as ER–EE
forcefield = app.ForceField(
    "amber14-all.xml",
    "amber14/tip3p.xml"
)

smirnoff = SMIRNOFFTemplateGenerator(
    molecules=[ngmn],
    forcefield="openff-2.2.1"
)
forcefield.registerTemplateGenerator(smirnoff.generator)

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=None
)

# Same positional restraints as ER–EE
restraint = mm.CustomExternalForce(
    "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
)
restraint.addGlobalParameter("k", 1000.0)

for name in ["x0", "y0", "z0"]:
    restraint.addPerParticleParameter(name)

positions_nm = modeller.positions.value_in_unit(unit.nanometer)
protein_atom_count = protein.topology.getNumAtoms()

for atom in modeller.topology.atoms():
    if atom.index < protein_atom_count and atom.element.symbol != "H":
        p = positions_nm[atom.index]
        restraint.addParticle(atom.index, [p.x, p.y, p.z])

system.addForce(restraint)

integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)

simulation = app.Simulation(
    modeller.topology,
    system,
    integrator,
    mm.Platform.getPlatformByName("CPU")
)

simulation.context.setPositions(modeller.positions)

initial = simulation.context.getState(
    getEnergy=True
).getPotentialEnergy()

print("Initial energy:", initial)
print("Minimizing...")

simulation.minimizeEnergy(
    tolerance=10 * unit.kilojoule_per_mole / unit.nanometer,
    maxIterations=1000
)

state = simulation.context.getState(
    getEnergy=True,
    getPositions=True
)

final = state.getPotentialEnergy()

print("Final energy:", final)

output = "/content/ER_NGMN_minimized_complex.pdb"

with open(output, "w") as f:
    app.PDBFile.writeFile(
        modeller.topology,
        state.getPositions(),
        f,
        keepIds=True
    )

print("Saved:", output)
print("Energy minimization complete.")

import os, sys, shutil, subprocess, re
import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D
from scipy.spatial.distance import cdist

print("=== ER–NGMN + EE SEQUENTIAL DOCKING ===")

# Minimized ER–NGMN complex
complex_path = "/content/ER_NGMN_minimized_complex.pdb"
receptor_path = "/content/ER_NGMN_receptor.pdbqt"
ligand_path = "/content/EE_ligand_B.pdbqt"
output_path = "/content/EE_docked_ER_NGMN.pdbqt"

assert os.path.isfile(complex_path), "Minimized ER–NGMN complex missing."

# Extract minimized protein and bound NGMN
protein_lines = []
ngmn_lines = []

for line in open(complex_path):
```

### Source lines 4306–4355
```python
    if not line.startswith(("ATOM  ", "HETATM")):
        continue
    if line[17:20].strip() in {"NGM", "UNL", "UNK"}:
        ngmn_lines.append(line)
    else:
        protein_lines.append(line)

assert len(protein_lines) > 1000, "Protein extraction failed."
assert len(ngmn_lines) >= 24, "Bound NGMN extraction failed."

with open("/content/ER_NGMN_protein_minimized.pdb", "w") as f:
    f.writelines(protein_lines)
    f.write("END\n")

# Convert minimized protein to PDBQT
obabel = shutil.which("obabel")
assert obabel, "Open Babel is not installed."

subprocess.run([
    obabel,
    "/content/ER_NGMN_protein_minimized.pdb",
    "-O", "/content/ER_NGMN_protein.pdbqt",
    "-xr"
], check=True, capture_output=True, text=True)

# Transfer minimized NGMN coordinates into its original
# chemically defined structure, preserving atom order.
ngmn = Chem.SDMolSupplier(
    "/content/NGMN_CB_Dock3_restored.sdf",
    removeHs=False
)[0]
assert ngmn is not None

heavy_indices = [
    a.GetIdx() for a in ngmn.GetAtoms()
    if a.GetSymbol() != "H"
]

ngmn_heavy = [
    line for line in ngmn_lines
    if line[76:78].strip().upper() != "H"
]

assert len(heavy_indices) == len(ngmn_heavy) == 24

# OpenMM retains the original ligand atom order
conf = ngmn.GetConformer()

for idx, line in zip(heavy_indices, ngmn_heavy):
    xyz = [
```

### Source lines 4399–4628
```python
pdbqt, ok, error = PDBQTWriterLegacy.write_string(setup)
assert ok, error

with open(ligand_path, "w") as f:
    f.write(pdbqt)

# Center docking box on minimized bound NGMN
ngmn_xyz = np.array([
    [
        float(line[30:38]),
        float(line[38:46]),
        float(line[46:54])
    ]
    for line in ngmn_heavy
])

center = ngmn_xyz.mean(axis=0)
size = np.ptp(ngmn_xyz, axis=0) + 18.0

vina = shutil.which("vina")
assert vina, "AutoDock Vina executable missing."

command = [
    vina,
    "--receptor", receptor_path,
    "--ligand", ligand_path,
    "--center_x", str(center[0]),
    "--center_y", str(center[1]),
    "--center_z", str(center[2]),
    "--size_x", str(size[0]),
    "--size_y", str(size[1]),
    "--size_z", str(size[2]),
    "--exhaustiveness", "16",
    "--num_modes", "10",
    "--out", output_path
]

result = subprocess.run(
    command,
    capture_output=True,
    text=True
)

print(result.stdout)

if result.returncode != 0:
    print(result.stderr)
    raise RuntimeError("EE docking failed.")

scores = re.findall(
    r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)",
    open(output_path).read()
)

print("\n=== RESULT ===")
print("Complex: ER–NGMN + EE")
print("Ligand A retained: NGMN")
print("Ligand B: EE")
print("Best saved docking score:",
      scores[0] if scores else "Check Vina table",
      "kcal/mol")
print("Docked EE coordinates saved:", output_path)

import os, re, shutil, subprocess, sys
import numpy as np
from rdkit import Chem

print("=== ER–NGMN + MPA SEQUENTIAL DOCKING ===")

# Reuse the minimized ER–NGMN receptor from the EE run
receptor = "/content/ER_NGMN_receptor.pdbqt"
mpa_sdf = "/content/mpa_prepared.sdf"
ligand = "/content/MPA_ligand_ER_NGMN.pdbqt"
output = "/content/MPA_docked_ER_NGMN.pdbqt"

assert os.path.isfile(receptor), "ER–NGMN receptor missing."
assert os.path.isfile(mpa_sdf), "Prepared MPA SDF missing."

# Prepare MPA as flexible Ligand B
import six
sys.modules["rdkit.six"] = six

from meeko import MoleculePreparation, PDBQTWriterLegacy

mpa = Chem.SDMolSupplier(mpa_sdf, removeHs=False)[0]
assert mpa is not None, "Could not read MPA."
assert mpa.GetNumConformers() > 0, "MPA lacks 3D coordinates."

mpa = Chem.AddHs(mpa, addCoords=True)
setup = MoleculePreparation().prepare(mpa)[0]
pdbqt, ok, error = PDBQTWriterLegacy.write_string(setup)
assert ok, error

with open(ligand, "w") as f:
    f.write(pdbqt)

# Use the same docking box as the ER–NGMN + EE run
# Center and size were reported in the previous Vina output
center = (9.33096, 46.27, 131.31)
size = (26.232, 25.99, 23.391)

vina = shutil.which("vina")
assert vina, "AutoDock Vina executable missing."

command = [
    vina,
    "--receptor", receptor,
    "--ligand", ligand,
    "--center_x", str(center[0]),
    "--center_y", str(center[1]),
    "--center_z", str(center[2]),
    "--size_x", str(size[0]),
    "--size_y", str(size[1]),
    "--size_z", str(size[2]),
    "--exhaustiveness", "16",
    "--num_modes", "10",
    "--out", output
]

result = subprocess.run(
    command,
    capture_output=True,
    text=True
)

print(result.stdout)

if result.returncode != 0:
    print(result.stderr)
    raise RuntimeError("MPA docking failed.")

scores = re.findall(
    r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)",
    open(output).read()
)

print("\n=== RESULT ===")
print("Complex: ER–NGMN + MPA")
print("Ligand A retained: NGMN")
print("Ligand B: MPA")
print(
    "Best saved docking score:",
    scores[0] if scores else "Check Vina table",
    "kcal/mol"
)
print("Docked MPA coordinates saved:", output)

import os, re, shutil, subprocess

print("=== ER–NGMN + LNG SEQUENTIAL DOCKING ===")

# Reuse the minimized ER–NGMN receptor and prepared LNG
receptor = "/content/ER_NGMN_receptor.pdbqt"
ligand = "/content/LNG_ligand.pdbqt"
output = "/content/LNG_docked_ER_NGMN.pdbqt"

assert os.path.isfile(receptor), "Minimized ER–NGMN receptor missing."
assert os.path.isfile(ligand), "Prepared LNG ligand missing."

vina = shutil.which("vina")
assert vina, "AutoDock Vina executable missing."

# Same docking box and settings as ER–NGMN + EE and MPA
command = [
    vina,
    "--receptor", receptor,
    "--ligand", ligand,
    "--center_x", "9.33096",
    "--center_y", "46.27",
    "--center_z", "131.31",
    "--size_x", "26.232",
    "--size_y", "25.99",
    "--size_z", "23.391",
    "--exhaustiveness", "16",
    "--num_modes", "10",
    "--out", output
]

result = subprocess.run(
    command,
    capture_output=True,
    text=True
)

print(result.stdout)

if result.returncode != 0:
    print(result.stderr)
    raise RuntimeError("LNG docking failed.")

scores = re.findall(
    r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)",
    open(output).read()
)

print("\n=== RESULT ===")
print("Complex: ER–NGMN + LNG")
print("Ligand A retained: NGMN")
print("Ligand B: LNG")
print(
    "Best saved docking score:",
    scores[0] if scores else "Check Vina table",
    "kcal/mol"
)
print("Docked LNG coordinates saved:", output)

import glob
from rdkit import Chem
from collections import Counter

print("=== ER–MPA PREPARATION CHECK ===")

# Find the uploaded CB-Dock3 files
pdb_files = glob.glob("/content/*mpa*complex.pdb")
mol2_files = glob.glob("/content/*mpa*.mol2")

assert pdb_files, "ER–MPA complex PDB not found."
assert mol2_files, "MPA MOL2 file not found."

complex_path = pdb_files[0]
mol2_path = mol2_files[0]

print("Complex:", complex_path)
print("Ligand:", mol2_path)

# Inspect docked MPA in the complex
ligand_lines = [
    line for line in open(complex_path)
    if line.startswith(("ATOM  ", "HETATM"))
    and line[17:20].strip() == "UNL"
```

### Source lines 4749–5048
```python
    pos = dock_conf.GetAtomPosition(dock_idx)
    ref_conf.SetAtomPosition(
        ref_idx,
        Point3D(pos.x, pos.y, pos.z)
    )

# Restore hydrogens at the docked geometry.
restored = Chem.AddHs(ref_heavy, addCoords=True)

output = "/content/MPA_CB_Dock3_restored.sdf"
Chem.MolToMolFile(restored, output)

print("Saved:", output)
print("Total atoms:", restored.GetNumAtoms())
print("Heavy atoms:", restored.GetNumHeavyAtoms())
print("Restoration complete.")

import openmm as mm
from openmm import app, unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Molecule
import numpy as np
import csv

print("=== ER–MPA ENERGY MINIMIZATION ===")

# Same prepared estrogen receptor used previously
protein = app.PDBFile("/content/er_protein_h.pdb")

# MPA restored at its CB-Dock3 docking coordinates
mpa = Molecule.from_file(
    "/content/MPA_CB_Dock3_restored.sdf",
    allow_undefined_stereo=False
)

# Combine ER and MPA
modeller = app.Modeller(protein.topology, protein.positions)

mpa_positions = (
    mpa.conformers[0].to("nanometer").magnitude
    * unit.nanometer
)

modeller.add(
    mpa.to_topology().to_openmm(),
    mpa_positions
)

# Same force fields as previous complexes
forcefield = app.ForceField(
    "amber14-all.xml",
    "amber14/tip3p.xml"
)

smirnoff = SMIRNOFFTemplateGenerator(
    molecules=[mpa],
    forcefield="openff-2.2.1"
)
forcefield.registerTemplateGenerator(smirnoff.generator)

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=None
)

# Same protein heavy-atom restraints
restraint = mm.CustomExternalForce(
    "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
)
restraint.addGlobalParameter("k", 1000.0)

for name in ["x0", "y0", "z0"]:
    restraint.addPerParticleParameter(name)

positions_nm = modeller.positions.value_in_unit(unit.nanometer)
protein_atom_count = protein.topology.getNumAtoms()

for atom in modeller.topology.atoms():
    if atom.index < protein_atom_count and atom.element.symbol != "H":
        p = positions_nm[atom.index]
        restraint.addParticle(
            atom.index,
            [p.x, p.y, p.z]
        )

system.addForce(restraint)

# Energy minimization
integrator = mm.VerletIntegrator(
    0.001 * unit.picoseconds
)

simulation = app.Simulation(
    modeller.topology,
    system,
    integrator,
    mm.Platform.getPlatformByName("CPU")
)

simulation.context.setPositions(modeller.positions)

initial = simulation.context.getState(
    getEnergy=True
).getPotentialEnergy()

print("Initial energy:", initial)
print("Minimizing...")

simulation.minimizeEnergy(
    tolerance=10 * unit.kilojoule_per_mole / unit.nanometer,
    maxIterations=1000
)

state = simulation.context.getState(
    getEnergy=True,
    getPositions=True
)

final = state.getPotentialEnergy()

print("Final energy:", final)

# Save minimized ER–MPA complex
output = "/content/ER_MPA_minimized_complex.pdb"

with open(output, "w") as f:
    app.PDBFile.writeFile(
        modeller.topology,
        state.getPositions(),
        f,
        keepIds=True
    )

print("Saved:", output)
print("Energy minimization complete.")

import os, re, shutil, subprocess, sys
import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D

print("=== ER–MPA + EE SEQUENTIAL DOCKING ===")

complex_path = "/content/ER_MPA_minimized_complex.pdb"
receptor = "/content/ER_MPA_receptor.pdbqt"
output = "/content/EE_docked_ER_MPA.pdbqt"

assert os.path.isfile(complex_path)
assert os.path.isfile("/content/MPA_CB_Dock3_restored.sdf")

# Separate minimized protein and bound MPA
protein_lines = []
mpa_lines = []

for line in open(complex_path):
    if not line.startswith(("ATOM  ", "HETATM")):
        continue
    if line[17:20].strip() in {"UNK", "UNL", "MPA"}:
        mpa_lines.append(line)
    else:
        protein_lines.append(line)

assert len(protein_lines) == 3819, "Unexpected protein atom count"
assert len(mpa_lines) == 62, "Unexpected MPA atom count"

with open("/content/ER_MPA_protein.pdb", "w") as f:
    f.writelines(protein_lines)
    f.write("END\n")

# Reconstruct minimized MPA coordinates
mpa = Chem.SDMolSupplier(
    "/content/MPA_CB_Dock3_restored.sdf",
    removeHs=False
)[0]

heavy = [a.GetIdx() for a in mpa.GetAtoms()
         if a.GetSymbol() != "H"]

mpa_heavy_lines = [
    line for line in mpa_lines
    if line[76:78].strip().upper() != "H"
]

assert len(heavy) == len(mpa_heavy_lines) == 28

conf = mpa.GetConformer()

for idx, line in zip(heavy, mpa_heavy_lines):
    xyz = [float(line[i:i+8]) for i in (30, 38, 46)]
    conf.SetAtomPosition(idx, Point3D(*xyz))

mpa = Chem.AddHs(Chem.RemoveHs(mpa), addCoords=True)
Chem.MolToMolFile(
    mpa, "/content/MPA_minimized_for_docking.sdf"
)

# Convert minimized protein and bound MPA to rigid PDBQT
obabel = shutil.which("obabel")
assert obabel, "Open Babel missing"

subprocess.run([
    obabel, "/content/ER_MPA_protein.pdb",
    "-O", "/content/ER_MPA_protein.pdbqt", "-xr"
], check=True, capture_output=True)

subprocess.run([
    obabel, "/content/MPA_minimized_for_docking.sdf",
    "-O", "/content/MPA_minimized.pdbqt", "-xr"
], check=True, capture_output=True)

with open(receptor, "w") as f:
    for path in [
        "/content/ER_MPA_protein.pdbqt",
        "/content/MPA_minimized.pdbqt"
    ]:
        for line in open(path):
            if line.startswith(("ATOM  ", "HETATM")):
                f.write(line)

# Reuse prepared EE ligand
ee_ligand = "/content/EE_ligand_B.pdbqt"
assert os.path.isfile(ee_ligand), "EE ligand PDBQT missing"

# Center box on minimized MPA
xyz = np.array([
    [float(line[i:i+8]) for i in (30, 38, 46)]
    for line in mpa_heavy_lines
])

center = xyz.mean(axis=0)
size = np.ptp(xyz, axis=0) + 18.0

vina = shutil.which("vina")
assert vina, "AutoDock Vina missing"

command = [
    vina,
    "--receptor", receptor,
    "--ligand", ee_ligand,
    "--center_x", str(center[0]),
    "--center_y", str(center[1]),
    "--center_z", str(center[2]),
    "--size_x", str(size[0]),
    "--size_y", str(size[1]),
    "--size_z", str(size[2]),
    "--exhaustiveness", "16",
    "--num_modes", "10",
    "--out", output
]

result = subprocess.run(
    command, capture_output=True, text=True
)

print(result.stdout)

if result.returncode != 0:
    print(result.stderr)
    raise RuntimeError("Docking failed")

scores = re.findall(
    r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)",
    open(output).read()
)

print("\n=== RESULT ===")
print("Combination: ER–MPA + EE")
print("Ligand A retained: MPA")
print("Ligand B: EE")
print("Best docking score:",
      scores[0] if scores else "Check Vina output",
      "kcal/mol")
print("Saved pose:", output)

import subprocess, shutil, os

print("=== ER–MPA + NGMN SEQUENTIAL DOCKING ===")

receptor = "/content/ER_MPA_receptor.pdbqt"
ligand = "/content/NGMN_ligand_B.pdbqt"
output = "/content/NGMN_docked_ER_MPA.pdbqt"

assert os.path.isfile(receptor), "ER–MPA receptor missing"
assert os.path.isfile(ligand), "NGMN ligand PDBQT missing"

vina = shutil.which("vina")
assert vina, "AutoDock Vina not found"

command = [
    vina,
    "--receptor", receptor,
    "--ligand", ligand,
    "--center_x", "25.7129",
    "--center_y", "62.2403",
    "--center_z", "130.214",
    "--size_x", "23.836",
    "--size_y", "26.487",
    "--size_z", "27.652",
    "--exhaustiveness", "16",
```

### Source lines 5049–5229
```python
    "--num_modes", "10",
    "--out", output
]

result = subprocess.run(
    command,
    capture_output=True,
    text=True
)

print(result.stdout)

if result.returncode != 0:
    print(result.stderr)
    raise RuntimeError("Docking failed")

print("\nSaved pose:", output)

import glob
import os

print("=== FINDING NGMN FILES ===")

for path in sorted(glob.glob("/content/*")):
    if "ngmn" in os.path.basename(path).lower():
        print(os.path.basename(path))

import os
import subprocess
import shutil

print("=== ER–MPA + NGMN SEQUENTIAL DOCKING ===")

receptor = "/content/ER_MPA_receptor.pdbqt"
ligand = "/content/NGMN_ligand.pdbqt"
output = "/content/NGMN_docked_ER_MPA.pdbqt"

assert os.path.isfile(receptor), "ER–MPA receptor missing"
assert os.path.isfile(ligand), "NGMN ligand missing"

vina = shutil.which("vina")
assert vina, "AutoDock Vina not found"

command = [
    vina,
    "--receptor", receptor,
    "--ligand", ligand,
    "--center_x", "25.7129",
    "--center_y", "62.2403",
    "--center_z", "130.214",
    "--size_x", "23.836",
    "--size_y", "26.487",
    "--size_z", "27.652",
    "--exhaustiveness", "16",
    "--num_modes", "10",
    "--out", output
]

result = subprocess.run(
    command,
    capture_output=True,
    text=True
)

print(result.stdout)

if result.returncode != 0:
    print(result.stderr)
    raise RuntimeError("Docking failed")

print("\n=== RESULT ===")
print("Ligand A retained: MPA")
print("Ligand B docked: NGMN")
print("Saved pose:", output)

import os
import glob
import shutil
import subprocess

print("=== ER–MPA + LNG SEQUENTIAL DOCKING ===")

receptor = "/content/ER_MPA_receptor.pdbqt"
output = "/content/LNG_docked_ER_MPA.pdbqt"

# Locate an existing prepared LNG ligand
candidates = [
    "/content/LNG_ligand_B.pdbqt",
    "/content/LNG_ligand.pdbqt",
    "/content/LNG_ligand_ER_NGMN.pdbqt"
]

ligand = next(
    (p for p in candidates if os.path.isfile(p)),
    None
)

if ligand is None:
    print("Available LNG files:")
    for p in glob.glob("/content/*"):
        if "lng" in os.path.basename(p).lower():
            print(os.path.basename(p))
    raise FileNotFoundError(
        "Prepared LNG ligand not found. "
        "Please send the file list above."
    )

assert os.path.isfile(receptor), "ER–MPA receptor missing"

vina = shutil.which("vina")
assert vina, "AutoDock Vina not found"

print("Using ligand:", ligand)

command = [
    vina,
    "--receptor", receptor,
    "--ligand", ligand,
    "--center_x", "25.7129",
    "--center_y", "62.2403",
    "--center_z", "130.214",
    "--size_x", "23.836",
    "--size_y", "26.487",
    "--size_z", "27.652",
    "--exhaustiveness", "16",
    "--num_modes", "10",
    "--out", output
]

result = subprocess.run(
    command,
    capture_output=True,
    text=True
)

print(result.stdout)

if result.returncode != 0:
    print(result.stderr)
    raise RuntimeError("Docking failed")

print("\n=== RESULT ===")
print("Ligand A retained: MPA")
print("Ligand B docked: LNG")
print("Saved pose:", output)

import glob
import os
from collections import Counter
from rdkit import Chem

print("=== ER–LNG PREPARATION CHECK ===")

# Find CB-Dock3 files
pdb_files = glob.glob("/content/*lng*complex.pdb")
mol2_files = glob.glob("/content/*lng*.mol2")

assert pdb_files, "ER–LNG complex PDB not found."
assert mol2_files, "LNG MOL2 file not found."

complex_path = pdb_files[0]
mol2_path = mol2_files[0]

print("Complex:", complex_path)
print("Ligand:", mol2_path)

# Inspect ligand atoms in the complex
ligand_lines = [
    line for line in open(complex_path)
    if line.startswith(("ATOM  ", "HETATM"))
    and line[17:20].strip() == "UNL"
]

elements = Counter(
    line[76:78].strip().upper()
    for line in ligand_lines
)

print("\nDocked ligand atoms:", len(ligand_lines))
print("Elements:", dict(elements))

```

### Source lines 5381–5557
```python
# Transfer docked heavy-atom coordinates
for dock_idx, ref_idx in enumerate(best_match):
    pos = dock_conf.GetAtomPosition(dock_idx)
    ref_conf.SetAtomPosition(
        ref_idx,
        Point3D(pos.x, pos.y, pos.z)
    )

# Restore explicit hydrogens
restored = Chem.AddHs(ref_heavy, addCoords=True)

output = "/content/LNG_CB_Dock3_restored.sdf"
Chem.MolToMolFile(restored, output)

print("Saved:", output)
print("Total atoms:", restored.GetNumAtoms())
print("Heavy atoms:", restored.GetNumHeavyAtoms())
print("Restoration complete.")

import openmm as mm
from openmm import app, unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Molecule

print("=== ER–LNG ENERGY MINIMIZATION ===")

# Same prepared ER used for the other complexes
protein = app.PDBFile("/content/er_protein_h.pdb")

# Restored LNG at its CB-Dock3 Pocket 1 coordinates
lng = Molecule.from_file(
    "/content/LNG_CB_Dock3_restored.sdf",
    allow_undefined_stereo=False
)

# Combine ER and LNG
modeller = app.Modeller(protein.topology, protein.positions)

lng_positions = (
    lng.conformers[0].to("nanometer").magnitude
    * unit.nanometer
)

modeller.add(
    lng.to_topology().to_openmm(),
    lng_positions
)

# Same protein and ligand force fields
forcefield = app.ForceField(
    "amber14-all.xml",
    "amber14/tip3p.xml"
)

smirnoff = SMIRNOFFTemplateGenerator(
    molecules=[lng],
    forcefield="openff-2.2.1"
)
forcefield.registerTemplateGenerator(smirnoff.generator)

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=None
)

# Restrain protein heavy atoms
restraint = mm.CustomExternalForce(
    "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
)
restraint.addGlobalParameter("k", 1000.0)

for name in ["x0", "y0", "z0"]:
    restraint.addPerParticleParameter(name)

positions_nm = modeller.positions.value_in_unit(unit.nanometer)
protein_atom_count = protein.topology.getNumAtoms()

for atom in modeller.topology.atoms():
    if (
        atom.index < protein_atom_count
        and atom.element.symbol != "H"
    ):
        p = positions_nm[atom.index]
        restraint.addParticle(
            atom.index,
            [p.x, p.y, p.z]
        )

system.addForce(restraint)

# Energy minimization
integrator = mm.VerletIntegrator(
    0.001 * unit.picoseconds
)

simulation = app.Simulation(
    modeller.topology,
    system,
    integrator,
    mm.Platform.getPlatformByName("CPU")
)

simulation.context.setPositions(modeller.positions)

initial = simulation.context.getState(
    getEnergy=True
).getPotentialEnergy()

print("Initial energy:", initial)
print("Minimizing...")

simulation.minimizeEnergy(
    tolerance=10 * unit.kilojoule_per_mole / unit.nanometer,
    maxIterations=1000
)

state = simulation.context.getState(
    getEnergy=True,
    getPositions=True
)

final = state.getPotentialEnergy()

print("Final energy:", final)

# Save minimized ER–LNG complex
output = "/content/ER_LNG_minimized_complex.pdb"

with open(output, "w") as f:
    app.PDBFile.writeFile(
        modeller.topology,
        state.getPositions(),
        f,
        keepIds=True
    )

print("Saved:", output)
print("Energy minimization complete.")

import os
import subprocess
import shutil
import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D

print("=== PREPARING MINIMIZED ER–LNG RECEPTOR ===")

complex_path = "/content/ER_LNG_minimized_complex.pdb"
ligand_sdf = "/content/LNG_CB_Dock3_restored.sdf"

assert os.path.isfile(complex_path)
assert os.path.isfile(ligand_sdf)

protein_lines = []
lng_lines = []

for line in open(complex_path):
    if not line.startswith(("ATOM  ", "HETATM")):
        continue
    if line[17:20].strip() in {"UNK", "UNL", "LNG"}:
        lng_lines.append(line)
    else:
        protein_lines.append(line)

print("Protein atoms:", len(protein_lines))
print("LNG atoms:", len(lng_lines))

assert len(protein_lines) == 3819, "Unexpected protein atom count"
assert len(lng_lines) == 51, "Unexpected LNG atom count"

protein_pdb = "/content/ER_LNG_protein.pdb"

with open(protein_pdb, "w") as f:
    f.writelines(protein_lines)
    f.write("END\n")
```

### Source lines 5604–5802
```python
        for line in open(path):
            if line.startswith(("ATOM  ", "HETATM")):
                f.write(line)

# Determine docking box from minimized LNG heavy atoms
heavy_xyz = np.array([
    list(conf.GetAtomPosition(atom.GetIdx()))
    for atom in lng.GetAtoms()
    if atom.GetSymbol() != "H"
])

center = heavy_xyz.mean(axis=0)
size = np.ptp(heavy_xyz, axis=0) + 18.0

print("\n=== DOCKING BOX ===")
print("Center:", np.round(center, 4))
print("Size:", np.round(size, 4))
print("\nSaved receptor:", receptor)
print("ER–LNG preparation complete.")

import os
import shutil
import subprocess

print("=== ER–LNG + EE SEQUENTIAL DOCKING ===")

receptor = "/content/ER_LNG_receptor.pdbqt"
ligand = "/content/EE_ligand_B.pdbqt"
output = "/content/EE_docked_ER_LNG.pdbqt"

assert os.path.isfile(receptor), "ER–LNG receptor missing"
assert os.path.isfile(ligand), "EE ligand missing"

vina = shutil.which("vina")
assert vina, "AutoDock Vina not found"

command = [
    vina,
    "--receptor", receptor,
    "--ligand", ligand,
    "--center_x", "9.4395",
    "--center_y", "46.227",
    "--center_z", "131.0945",
    "--size_x", "23.377",
    "--size_y", "27.651",
    "--size_z", "23.314",
    "--exhaustiveness", "16",
    "--num_modes", "10",
    "--out", output
]

result = subprocess.run(
    command,
    capture_output=True,
    text=True
)

print(result.stdout)

if result.returncode != 0:
    print(result.stderr)
    raise RuntimeError("Docking failed")

print("\n=== RESULT ===")
print("Ligand A retained: LNG")
print("Ligand B docked: EE")
print("Saved pose:", output)

import os
import shutil
import subprocess

print("=== ER–LNG + NGMN SEQUENTIAL DOCKING ===")

receptor = "/content/ER_LNG_receptor.pdbqt"
ligand = "/content/NGMN_ligand.pdbqt"
output = "/content/NGMN_docked_ER_LNG.pdbqt"

assert os.path.isfile(receptor), "ER–LNG receptor missing"
assert os.path.isfile(ligand), "NGMN ligand missing"

vina = shutil.which("vina")
assert vina, "AutoDock Vina not found"

command = [
    vina,
    "--receptor", receptor,
    "--ligand", ligand,
    "--center_x", "9.4395",
    "--center_y", "46.227",
    "--center_z", "131.0945",
    "--size_x", "23.377",
    "--size_y", "27.651",
    "--size_z", "23.314",
    "--exhaustiveness", "16",
    "--num_modes", "10",
    "--out", output
]

result = subprocess.run(
    command,
    capture_output=True,
    text=True
)

print(result.stdout)

if result.returncode != 0:
    print(result.stderr)
    raise RuntimeError("Docking failed")

print("\n=== RESULT ===")
print("Ligand A retained: LNG")
print("Ligand B docked: NGMN")
print("Saved pose:", output)

import os
import shutil
import subprocess

print("=== ER–LNG + MPA SEQUENTIAL DOCKING ===")

receptor = "/content/ER_LNG_receptor.pdbqt"
ligand = "/content/MPA_ligand_ER_NGMN.pdbqt"
output = "/content/MPA_docked_ER_LNG.pdbqt"

assert os.path.isfile(receptor), "ER–LNG receptor missing"
assert os.path.isfile(ligand), "MPA ligand missing"

vina = shutil.which("vina")
assert vina, "AutoDock Vina not found"

command = [
    vina,
    "--receptor", receptor,
    "--ligand", ligand,
    "--center_x", "9.4395",
    "--center_y", "46.227",
    "--center_z", "131.0945",
    "--size_x", "23.377",
    "--size_y", "27.651",
    "--size_z", "23.314",
    "--exhaustiveness", "16",
    "--num_modes", "10",
    "--out", output
]

result = subprocess.run(
    command,
    capture_output=True,
    text=True
)

print(result.stdout)

if result.returncode != 0:
    print(result.stderr)
    raise RuntimeError("Docking failed")

print("\n=== RESULT ===")
print("Ligand A retained: LNG")
print("Ligand B docked: MPA")
print("Saved pose:", output)

import os

ligands = ["EE", "NGMN", "MPA", "LNG"]

print("=== BINDING POSITION FILE CHECK ===\n")

for ligand_a in ligands:
    complex_path = f"/content/ER_{ligand_a}_minimized_complex.pdb"

    print(f"ER–{ligand_a}")
    print("  Minimized complex:",
          "FOUND" if os.path.isfile(complex_path) else "MISSING")

    for ligand_b in ligands:
        if ligand_a == ligand_b:
            continue

        pose_path = (
            f"/content/{ligand_b}_docked_ER_{ligand_a}.pdbqt"
        )

        print(
            f"  {ligand_a} + {ligand_b}:",
            "FOUND" if os.path.isfile(pose_path) else "MISSING"
        )

    print()

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

ligands = ["EE", "NGMN", "MPA", "LNG"]
results = []

```

### Source lines 5853–5955
```python

for ligand_a in ligands:
    a_path = f"/content/ER_{ligand_a}_minimized_complex.pdb"
    xyz_a = get_ligand_a_coordinates(a_path)

    for ligand_b in ligands:
        if ligand_a == ligand_b:
            continue

        b_path = f"/content/{ligand_b}_docked_ER_{ligand_a}.pdbqt"
        xyz_b = get_ligand_b_coordinates(b_path)

        assert xyz_a.shape[1] == xyz_b.shape[1] == 3
        assert len(xyz_a) > 0 and len(xyz_b) > 0

        distances = cdist(xyz_a, xyz_b)

        centroid_distance = np.linalg.norm(
            xyz_a.mean(axis=0) - xyz_b.mean(axis=0)
        )

        results.append({
            "Ligand A": ligand_a,
            "Ligand B": ligand_b,
            "Centroid distance (Å)": round(centroid_distance, 3),
            "Minimum heavy-atom distance (Å)": round(
                distances.min(), 3
            ),
            "A heavy atoms": len(xyz_a),
            "B heavy atoms": len(xyz_b)
        })

df = pd.DataFrame(results)

print("=== ALL 12 LIGAND POSITION COMPARISONS ===")
print(df.to_string(index=False))

output = "/content/ligand_binding_distances.csv"
df.to_csv(output, index=False)

print("\nSaved:", output)

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

ligands = ["EE", "NGMN", "MPA", "LNG"]
cutoff = 4.0
results = []

def nearby_residues(complex_path, ligand_xyz):
    protein_xyz = []
    residue_ids = []

    for line in open(complex_path):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue

        # Protein residues are standard ATOM records
        if not line.startswith("ATOM  "):
            continue

        element = line[76:78].strip().upper()
        if element in {"H", "D"}:
            continue

        xyz = [
            float(line[i:i+8])
            for i in (30, 38, 46)
        ]

        residue = (
            line[21].strip(),
            line[22:26].strip(),
            line[17:20].strip()
        )

        protein_xyz.append(xyz)
        residue_ids.append(residue)

    distances = cdist(
        np.array(protein_xyz),
        ligand_xyz
    )

    close_indices = np.where(
        distances.min(axis=1) <= cutoff
    )[0]

    return {
        residue_ids[i]
        for i in close_indices
    }

def format_residues(residues):
    return ", ".join(
        f"{name}{number}"
        for chain, number, name in sorted(
            residues,
            key=lambda x: (
                x[0],
                int(x[1]) if x[1].lstrip("-").isdigit() else 0
            )
```

### Source lines 5973–6061
```python
        shared = residues_a & residues_b

        results.append({
            "Ligand A": ligand_a,
            "Ligand B": ligand_b,
            "A residues": format_residues(residues_a),
            "B residues": format_residues(residues_b),
            "Shared residues": format_residues(shared),
            "Shared count": len(shared)
        })

df_sites = pd.DataFrame(results)

print("=== BINDING-SITE RESIDUE COMPARISON ===\n")

for _, row in df_sites.iterrows():
    print(f"{row['Ligand A']} + {row['Ligand B']}")
    print("Ligand A residues:", row["A residues"])
    print("Ligand B residues:", row["B residues"])
    print("Shared residues:", row["Shared residues"] or "None")
    print("Shared count:", row["Shared count"])
    print()

output = "/content/binding_site_residue_comparison.csv"
df_sites.to_csv(output, index=False)

print("Saved:", output)

import sys, subprocess, importlib.util

if importlib.util.find_spec("py3Dmol") is None:
    subprocess.run(
        [sys.executable, "-m", "pip", "-q", "install", "py3Dmol"],
        check=True
    )

import py3Dmol

complex_path = "/content/ER_MPA_minimized_complex.pdb"
pose_path = "/content/EE_docked_ER_MPA.pdbqt"

# Separate protein and retained MPA using known atom count
atom_lines = [
    line for line in open(complex_path)
    if line.startswith(("ATOM  ", "HETATM"))
]

protein_pdb = "".join(atom_lines[:3819]) + "END\n"
mpa_pdb = "".join(atom_lines[3819:]) + "END\n"

# Extract first Vina docking mode
pose_lines = []
for line in open(pose_path):
    if line.startswith("ENDMDL"):
        break
    if line.startswith(("ATOM  ", "HETATM")):
        pose_lines.append(line)

ee_pdb = "".join(pose_lines) + "END\n"

view = py3Dmol.view(width=850, height=600)

# Estrogen receptor: transparent cartoon
view.addModel(protein_pdb, "pdb")
view.setStyle(
    {"model": 0},
    {"cartoon": {"color": "lightgray", "opacity": 0.65}}
)

# Retained MPA: blue
view.addModel(mpa_pdb, "pdb")
view.setStyle(
    {"model": 1},
    {"stick": {"colorscheme": "blueCarbon", "radius": 0.18}}
)

# Docked EE: orange
view.addModel(ee_pdb, "pdb")
view.setStyle(
    {"model": 2},
    {"stick": {"colorscheme": "orangeCarbon", "radius": 0.18}}
)

view.zoomTo({"model": 1})
view.zoom(0.8)
view.setBackgroundColor("white")
view.show()

print("Blue = MPA (Ligand A, retained)")
```

### Source lines 6097–6396
```python

        center_a = xyz_a.mean(axis=0)
        center_b = xyz_b.mean(axis=0)

        results.append({
            "Ligand A": ligand_a,
            "Ligand B": ligand_b,
            "A-to-reference center (Å)": round(
                np.linalg.norm(center_a - reference_center), 3
            ),
            "B-to-reference center (Å)": round(
                np.linalg.norm(center_b - reference_center), 3
            ),
            "B-to-reference minimum (Å)": round(
                cdist(xyz_b, orthosteric_reference).min(), 3
            )
        })

df_reference = pd.DataFrame(results)

print("=== BINDING POSITIONS RELATIVE TO EE REFERENCE ===")
print(df_reference.to_string(index=False))

output = "/content/orthosteric_reference_comparison.csv"
df_reference.to_csv(output, index=False)

print("\nSaved:", output)

import numpy as np

ligands = ["EE", "NGMN", "MPA", "LNG"]

def protein_ca_coordinates(path):
    coords = {}

    for line in open(path):
        if not line.startswith("ATOM  "):
            continue

        if line[12:16].strip() != "CA":
            continue

        key = (
            line[21].strip(),
            line[22:26].strip(),
            line[26].strip()
        )

        coords[key] = np.array([
            float(line[i:i+8])
            for i in (30, 38, 46)
        ])

    return coords

reference = protein_ca_coordinates(
    "/content/ER_EE_minimized_complex.pdb"
)

print("=== RECEPTOR ALIGNMENT CHECK ===")

for ligand in ligands:
    path = f"/content/ER_{ligand}_minimized_complex.pdb"
    current = protein_ca_coordinates(path)

    common = sorted(set(reference) & set(current))

    ref_xyz = np.array([reference[k] for k in common])
    cur_xyz = np.array([current[k] for k in common])

    # RMSD in the existing coordinate frame
    direct_rmsd = np.sqrt(
        np.mean(np.sum((ref_xyz - cur_xyz)**2, axis=1))
    )

    # Kabsch alignment to identify rigid-body displacement
    p = cur_xyz - cur_xyz.mean(axis=0)
    q = ref_xyz - ref_xyz.mean(axis=0)

    u, s, vt = np.linalg.svd(p.T @ q)

    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)

    rotation = u @ correction @ vt
    aligned = p @ rotation + ref_xyz.mean(axis=0)

    aligned_rmsd = np.sqrt(
        np.mean(np.sum((aligned - ref_xyz)**2, axis=1))
    )

    print(
        f"ER–{ligand}: "
        f"{len(common)} matched Cα atoms | "
        f"Direct RMSD = {direct_rmsd:.3f} Å | "
        f"Aligned RMSD = {aligned_rmsd:.3f} Å"
    )

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

ligands = ["EE", "NGMN", "MPA", "LNG"]

reference = get_ligand_a_coordinates(
    "/content/ER_EE_minimized_complex.pdb"
)

results = []

def reference_proximity(xyz):
    distances = cdist(xyz, reference)
    minimum = distances.min(axis=1)

    return {
        "Closest distance (Å)": round(float(minimum.min()), 3),
        "Atoms within 4 Å": int(np.sum(minimum <= 4.0)),
        "Heavy atoms": len(xyz),
        "Percent near reference": round(
            100 * np.mean(minimum <= 4.0), 1
        )
    }

for ligand_a in ligands:
    complex_path = f"/content/ER_{ligand_a}_minimized_complex.pdb"
    xyz_a = get_ligand_a_coordinates(complex_path)

    for ligand_b in ligands:
        if ligand_a == ligand_b:
            continue

        xyz_b = get_ligand_b_coordinates(
            f"/content/{ligand_b}_docked_ER_{ligand_a}.pdbqt"
        )

        a_metrics = reference_proximity(xyz_a)
        b_metrics = reference_proximity(xyz_b)

        results.append({
            "Ligand A": ligand_a,
            "Ligand B": ligand_b,
            "A near reference (%)": a_metrics["Percent near reference"],
            "B near reference (%)": b_metrics["Percent near reference"],
            "B closest distance (Å)": b_metrics["Closest distance (Å)"]
        })

df_overlap = pd.DataFrame(results)

print("=== EE REFERENCE POCKET PROXIMITY ===")
print(df_overlap.to_string(index=False))

df_overlap.to_csv(
    "/content/EE_reference_pocket_proximity.csv",
    index=False
)

print("\nSaved: /content/EE_reference_pocket_proximity.csv")

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

ligands = ["EE", "NGMN", "MPA", "LNG"]
results = []

for ligand_a in ligands:
    xyz_a = get_ligand_a_coordinates(
        f"/content/ER_{ligand_a}_minimized_complex.pdb"
    )

    for ligand_b in ligands:
        if ligand_a == ligand_b:
            continue

        xyz_b = get_ligand_b_coordinates(
            f"/content/{ligand_b}_docked_ER_{ligand_a}.pdbqt"
        )

        distances = cdist(xyz_a, xyz_b)

        results.append({
            "Ligand A": ligand_a,
            "Ligand B": ligand_b,
            "Minimum distance (Å)": round(float(distances.min()), 3),
            "Atom pairs < 2.0 Å": int(np.sum(distances < 2.0)),
            "Atom pairs < 2.5 Å": int(np.sum(distances < 2.5)),
            "Atom pairs < 3.0 Å": int(np.sum(distances < 3.0))
        })

df_clashes = pd.DataFrame(results)

print("=== LIGAND–LIGAND STERIC COMPATIBILITY ===")
print(df_clashes.to_string(index=False))

output = "/content/ligand_steric_compatibility.csv"
df_clashes.to_csv(output, index=False)

print("\nSaved:", output)

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

ligands = ["EE", "NGMN", "MPA", "LNG"]
results = []

for a in ligands:
    path = f"/content/ER_{a}_minimized_complex.pdb"

    # Protein heavy atoms only
    protein_xyz = []

    for line in open(path):
        if not line.startswith("ATOM  "):
            continue

        element = line[76:78].strip().upper()

        if element in {"H", "D"}:
            continue

        protein_xyz.append([
            float(line[i:i+8])
            for i in (30, 38, 46)
        ])

    protein_xyz = np.array(protein_xyz)

    for b in ligands:
        if a == b:
            continue

        xyz_b = get_ligand_b_coordinates(
            f"/content/{b}_docked_ER_{a}.pdbqt"
        )

        distances = cdist(xyz_b, protein_xyz)

        results.append({
            "Ligand A": a,
            "Ligand B": b,
            "Closest B–receptor distance (Å)": round(
                float(distances.min()), 3
            ),
            "Pairs < 2.0 Å": int(np.sum(distances < 2.0)),
            "Pairs < 2.5 Å": int(np.sum(distances < 2.5))
        })

df = pd.DataFrame(results)

print("=== LIGAND B–RECEPTOR STERIC SCREEN ===")
print(df.to_string(index=False))

output = "/content/ligand_receptor_steric_screen.csv"
df.to_csv(output, index=False)

print("\nSaved:", output)

import os
import numpy as np
from scipy.spatial.distance import cdist

print("=== ER–EE + NGMN MODE 2 STERIC CHECK ===")

receptor_path = "/content/ER_EE_minimized_complex.pdb"
pose_path = "/content/NGMN_docked_ER_EE.pdbqt"

assert os.path.isfile(receptor_path)
assert os.path.isfile(pose_path)

protein_xyz = []

for line in open(receptor_path):
    if not line.startswith("ATOM  "):
        continue

    element = line[76:78].strip().upper()
    if element in {"H", "D"}:
        continue

    protein_xyz.append([
        float(line[i:i+8])
        for i in (30, 38, 46)
    ])

protein_xyz = np.array(protein_xyz)

# Extract Vina Mode 2
mode2_lines = []
current_mode = None

for line in open(pose_path):
    if line.startswith("MODEL"):
        current_mode = int(line.split()[1])
        continue

    if line.startswith("ENDMDL"):
        if current_mode == 2:
            break
        current_mode = None
```

### Source lines 6397–6415
```python
        continue

    if current_mode == 2 and line.startswith(("ATOM  ", "HETATM")):
        atom_type = line[77:79].strip().upper()
        if atom_type not in {"H", "HD", "D"}:
            mode2_lines.append(line)

assert mode2_lines, "Mode 2 not found in the docking file"

mode2_xyz = np.array([
    [float(line[i:i+8]) for i in (30, 38, 46)]
    for line in mode2_lines
])

distances = cdist(mode2_xyz, protein_xyz)

print("Mode 2 heavy atoms:", len(mode2_xyz))
print("Closest NGMN–receptor distance:",
      round(float(distances.min()), 3), "Å")
```

### Source lines 6499–6793
```python

for a in ligands:
    complex_path = f"/content/ER_{a}_minimized_complex.pdb"

    xyz_a = get_ligand_a_coordinates(complex_path)
    protein_xyz, protein_residues = get_protein_atoms(complex_path)

    for b in ligands:
        if a == b:
            continue

        mode = 2 if (a == "EE" and b == "NGMN") else 1
        pose_path = f"/content/{b}_docked_ER_{a}.pdbqt"
        xyz_b = get_pose(pose_path, mode)

        ab = cdist(xyz_a, xyz_b)
        centroid = np.linalg.norm(
            xyz_a.mean(axis=0) - xyz_b.mean(axis=0)
        )

        distance_rows.append({
            "Ligand A": a,
            "Ligand B": b,
            "Selected mode": mode,
            "Centroid distance (Å)": round(float(centroid), 3),
            "Minimum A-B distance (Å)": round(float(ab.min()), 3)
        })

        dist_a = cdist(xyz_a, protein_xyz)
        dist_b = cdist(xyz_b, protein_xyz)

        contacts_a = set(protein_residues[
            np.any(dist_a <= 4.0, axis=0)
        ])
        contacts_b = set(protein_residues[
            np.any(dist_b <= 4.0, axis=0)
        ])

        shared = contacts_a & contacts_b

        contact_rows.append({
            "Ligand A": a,
            "Ligand B": b,
            "Selected mode": mode,
            "A contact residues": ", ".join(sorted(contacts_a)),
            "B contact residues": ", ".join(sorted(contacts_b)),
            "Shared residues": ", ".join(sorted(shared)),
            "Number shared": len(shared)
        })

df_dist = pd.DataFrame(distance_rows)
df_contacts = pd.DataFrame(contact_rows)

print("=== UPDATED SPATIAL DISTANCES ===")
print(df_dist.to_string(index=False))

print("\n=== RECEPTOR CONTACT COMPARISON ===")
print(df_contacts[
    ["Ligand A", "Ligand B", "Selected mode",
     "Number shared", "Shared residues"]
].to_string(index=False))

df_dist.to_csv(
    "/content/final_ligand_spatial_distances.csv",
    index=False
)

df_contacts.to_csv(
    "/content/final_receptor_contact_comparison.csv",
    index=False
)

print("\nSaved both updated CSV files.")

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.cluster.hierarchy import linkage, fcluster

ligands = ["EE", "LNG", "NGMN", "MPA"]
placements = []

for a in ligands:
    complex_path = f"/content/ER_{a}_minimized_complex.pdb"
    xyz_a = get_ligand_a_coordinates(complex_path)

    placements.append({
        "Combination": f"{a} as A",
        "Ligand": a,
        "Role": "A",
        "Coordinates": xyz_a.mean(axis=0)
    })

    for b in ligands:
        if a == b:
            continue

        mode = 2 if (a == "EE" and b == "NGMN") else 1

        xyz_b = get_pose(
            f"/content/{b}_docked_ER_{a}.pdbqt",
            mode
        )

        placements.append({
            "Combination": f"{a} → {b}",
            "Ligand": b,
            "Role": "B",
            "Coordinates": xyz_b.mean(axis=0)
        })

coords = np.array([p["Coordinates"] for p in placements])

# Cluster ligand centroids at a provisional 6 Å threshold.
# This is spatial grouping, not validated pocket assignment.
clusters = fcluster(
    linkage(coords, method="complete"),
    t=6.0,
    criterion="distance"
)

rows = []

for p, cluster in zip(placements, clusters):
    x, y, z = p["Coordinates"]

    rows.append({
        "Combination": p["Combination"],
        "Ligand": p["Ligand"],
        "Role": p["Role"],
        "Centroid X": round(float(x), 3),
        "Centroid Y": round(float(y), 3),
        "Centroid Z": round(float(z), 3),
        "Spatial group": int(cluster)
    })

df_sites = pd.DataFrame(rows)

print("=== LIGAND POSITION GROUPS ===")
print(df_sites.to_string(index=False))

print("\n=== GROUP SUMMARY ===")
print(
    df_sites.groupby("Spatial group")
    .agg(
        Count=("Ligand", "size"),
        Ligands=("Ligand", lambda x: ", ".join(sorted(set(x)))),
        Roles=("Role", lambda x: ", ".join(sorted(set(x))))
    )
    .to_string()
)

df_sites.to_csv(
    "/content/ligand_spatial_groups.csv",
    index=False
)

print("\nSaved: /content/ligand_spatial_groups.csv")

import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist

ligands = ["EE", "LNG", "NGMN", "MPA"]
records = []

for a in ligands:
    path = f"/content/ER_{a}_minimized_complex.pdb"
    protein_xyz, residues = get_protein_atoms(path)

    for role, b in [("A", a)] + [
        ("B", x) for x in ligands if x != a
    ]:
        if role == "A":
            xyz = get_ligand_a_coordinates(path)
            combination = f"{a} as A"
        else:
            mode = 2 if (a == "EE" and b == "NGMN") else 1
            xyz = get_pose(
                f"/content/{b}_docked_ER_{a}.pdbqt",
                mode
            )
            combination = f"{a} → {b}"

        group = int(df_sites.loc[
            (df_sites["Combination"] == combination) &
            (df_sites["Role"] == role),
            "Spatial group"
        ].iloc[0])

        d = cdist(xyz, protein_xyz)
        contacts = set(residues[np.any(d <= 4.0, axis=0)])

        records.append({
            "Group": group,
            "Combination": combination,
            "Role": role,
            "Ligand": b,
            "Contacts": contacts
        })

print("=== RECEPTOR CONTACT FINGERPRINTS BY SPATIAL GROUP ===")

for group in sorted(set(r["Group"] for r in records)):
    members = [r for r in records if r["Group"] == group]

    common = set.intersection(
        *(r["Contacts"] for r in members)
    )

    union = set.union(
        *(r["Contacts"] for r in members)
    )

    print(f"\nGROUP {group}")
    print("Placements:", len(members))
    print("Common residues:", ", ".join(sorted(common)) or "None")
    print("All contacted residues:", ", ".join(sorted(union)))

summary = pd.DataFrame([
    {
        "Group": r["Group"],
        "Combination": r["Combination"],
        "Role": r["Role"],
        "Ligand": r["Ligand"],
        "Contacts": ", ".join(sorted(r["Contacts"]))
    }
    for r in records
])

summary.to_csv(
    "/content/spatial_group_receptor_contacts.csv",
    index=False
)

print("\nSaved: /content/spatial_group_receptor_contacts.csv")

import os
from pathlib import Path

print("=== RECEPTOR SOURCE IDENTIFICATION ===")

print("\nAvailable PDB and CIF files:")
for folder in ["/content", "/content/drive/MyDrive"]:
    if not os.path.exists(folder):
        continue
    for path in sorted(Path(folder).glob("*.pdb")):
        print(path)
    for path in sorted(Path(folder).glob("*.cif")):
        print(path)

print("\n=== MINIMIZED RECEPTOR DETAILS ===")

path = "/content/ER_EE_minimized_complex.pdb"
residues = []
seen = set()

with open(path) as f:
    for line in f:
        if not line.startswith("ATOM  "):
            continue

        key = (
            line[21:22].strip(),
            line[22:26].strip(),
            line[26:27].strip()
        )

        if key not in seen:
            seen.add(key)
            residues.append((
                key[0],
                key[1],
                line[17:20].strip()
            ))

print("Number of protein residues:", len(residues))
print("First 10:", residues[:10])
print("Last 10:", residues[-10:])

print("\nHeader information from available PDB files:")

for path in sorted(Path("/content").glob("*.pdb")):
    if "minimized_complex" in path.name:
        continue

    headers = []
    with open(path, errors="ignore") as f:
        for line in f:
            if line.startswith((
                "HEADER", "TITLE ", "COMPND",
                "SOURCE", "DBREF"
            )):
                headers.append(line.strip())

```

### Source lines 7033–7121
```python
assert best is not None, "No suitable receptor alignment found."

estradiol = [
    r for r in best["structure_chain"]
    if r.resname == "EST"
]

assert estradiol, "No crystallographic EST found in matched chain."

est_xyz = np.array([
    atom.coord for atom in estradiol[0]
    if atom.element.upper() not in {"H", "D"}
])

# Transform experimental EST into prepared receptor coordinates.
est_aligned = (
    est_xyz @ best["rotation"] + best["translation"]
)

np.save("/content/experimental_orthosteric_EST_aligned.npy", est_aligned)

print("=== EXPERIMENTAL ORTHOSTERIC SITE ALIGNMENT ===")
print("Matched 1ERE chain:", best["chain"])
print("Matched protein residues:", best["matches"])
print("Protein alignment RMSD (Å):", round(best["rms"], 3))
print("Experimental EST heavy atoms:", len(est_aligned))
print("Aligned EST centroid:",
      np.round(est_aligned.mean(axis=0), 3))
print("\nSaved experimental orthosteric coordinates.")

# Commented out IPython magic to ensure Python compatibility.

# %pip install biopython -q

from Bio.PDB import PDBParser, Superimposer
from Bio.Align import PairwiseAligner

print("Biopython installed successfully!")

import numpy as np
from Bio.PDB import PDBParser, Superimposer
from scipy.spatial.distance import cdist
import urllib.request

url = "https://files.rcsb.org/download/1ERE.pdb"
reference_path = "/content/1ERE_experimental.pdb"
urllib.request.urlretrieve(url, reference_path)

parser = PDBParser(QUIET=True)
experimental = parser.get_structure("1ERE", reference_path)
prepared = parser.get_structure(
    "prepared", "/content/ER_EE_minimized_complex.pdb"
)

# The prepared receptor retains residue identities and order,
# but has gaps and altered chain identifiers.
def get_ca(structure, chain_id=None):
    atoms = []
    for chain in structure[0]:
        if chain_id is not None and chain.id != chain_id:
            continue
        for residue in chain:
            if residue.id[0] == " " and "CA" in residue:
                atoms.append(residue["CA"])
    return atoms

prepared_ca = get_ca(prepared)
prepared_seq = [a.get_parent().resname for a in prepared_ca]

best = None

for chain in experimental[0]:
    reference_ca = get_ca(experimental, chain.id)
    reference_seq = [a.get_parent().resname for a in reference_ca]

    # Align residue sequences, accommodating missing stretches.
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -8
    aligner.extend_gap_score = -0.5

    from Bio.SeqUtils import seq1
    s1 = "".join(seq1(x, custom_map={"UNK":"X"}) for x in prepared_seq)
    s2 = "".join(seq1(x, custom_map={"UNK":"X"}) for x in reference_seq)

    alignment = aligner.align(s1, s2)[0]
```

### Source lines 7152–7451
```python
assert best is not None, "No suitable receptor alignment found."

estradiol = [
    r for r in best["structure_chain"]
    if r.resname == "EST"
]

assert estradiol, "No crystallographic EST found in matched chain."

est_xyz = np.array([
    atom.coord for atom in estradiol[0]
    if atom.element.upper() not in {"H", "D"}
])

# Transform experimental EST into prepared receptor coordinates.
est_aligned = (
    est_xyz @ best["rotation"] + best["translation"]
)

np.save("/content/experimental_orthosteric_EST_aligned.npy", est_aligned)

print("=== EXPERIMENTAL ORTHOSTERIC SITE ALIGNMENT ===")
print("Matched 1ERE chain:", best["chain"])
print("Matched protein residues:", best["matches"])
print("Protein alignment RMSD (Å):", round(best["rms"], 3))
print("Experimental EST heavy atoms:", len(est_aligned))
print("Aligned EST centroid:",
      np.round(est_aligned.mean(axis=0), 3))
print("\nSaved experimental orthosteric coordinates.")

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

# Experimental orthosteric ligand coordinates
est = np.load(
    "/content/experimental_orthosteric_EST_aligned.npy"
)

ligands = ["EE", "LNG", "NGMN", "MPA"]
records = []

# Four initial ligand A placements
for A in ligands:
    xyz = get_ligand_a_coordinates(
        f"/content/ER_{A}_minimized_complex.pdb"
    )

    distances = cdist(xyz, est)

    records.append({
        "Placement": f"{A} (Ligand A)",
        "Ligand_A": A,
        "Ligand_B": "-",
        "Centroid_distance_to_EST_A":
            np.linalg.norm(xyz.mean(axis=0) - est.mean(axis=0)),
        "Minimum_heavy_atom_distance_to_EST_A":
            distances.min(),
        "Ligand_atoms_within_4A_of_EST":
            int(np.sum(distances.min(axis=1) <= 4.0))
    })

# Twelve sequential ligand B placements
for A in ligands:
    for B in ligands:
        if A == B:
            continue

        mode = 2 if (A == "EE" and B == "NGMN") else 1

        xyz = get_pose(
            f"/content/{B}_docked_ER_{A}.pdbqt",
            mode=mode
        )

        distances = cdist(xyz, est)

        records.append({
            "Placement": f"{A} → {B} (Ligand B)",
            "Ligand_A": A,
            "Ligand_B": B,
            "Centroid_distance_to_EST_A":
                np.linalg.norm(xyz.mean(axis=0) - est.mean(axis=0)),
            "Minimum_heavy_atom_distance_to_EST_A":
                distances.min(),
            "Ligand_atoms_within_4A_of_EST":
                int(np.sum(distances.min(axis=1) <= 4.0))
        })

df = pd.DataFrame(records)

df = df.sort_values(
    "Centroid_distance_to_EST_A"
).reset_index(drop=True)

df.to_csv(
    "/content/experimental_orthosteric_proximity.csv",
    index=False
)

print("=== EXPERIMENTAL ORTHOSTERIC POCKET COMPARISON ===")
print(df.round(3).to_string(index=False))
print("\nSaved: /content/experimental_orthosteric_proximity.csv")

import pandas as pd
import numpy as np

spatial = pd.read_csv(
    "/content/final_ligand_spatial_distances.csv"
)

orthosteric = pd.read_csv(
    "/content/experimental_orthosteric_proximity.csv"
)

# These groups were previously determined by
# complete-linkage clustering at a 6 Å threshold.
a_groups = {
    "EE": 5,
    "LNG": 5,
    "NGMN": 5,
    "MPA": 2
}

b_groups = {
    ("EE", "LNG"): 3,
    ("EE", "NGMN"): 3,
    ("EE", "MPA"): 3,
    ("LNG", "EE"): 4,
    ("LNG", "NGMN"): 3,
    ("LNG", "MPA"): 3,
    ("NGMN", "EE"): 4,
    ("NGMN", "LNG"): 3,
    ("NGMN", "MPA"): 3,
    ("MPA", "EE"): 1,
    ("MPA", "LNG"): 1,
    ("MPA", "NGMN"): 1
}

# Identify columns without assuming their exact spelling.
print("Spatial CSV columns:", spatial.columns.tolist())

a_col = next(c for c in spatial.columns
             if c.strip().lower() in ["a", "ligand_a"])
b_col = next(c for c in spatial.columns
             if c.strip().lower() in ["b", "ligand_b"])

distance_cols = [
    c for c in spatial.columns
    if "min" in c.lower() and
    ("heavy" in c.lower() or "distance" in c.lower())
]
assert distance_cols, "Minimum-distance column not found."
distance_col = distance_cols[0]

results = []

for _, row in spatial.iterrows():
    A = str(row[a_col]).strip()
    B = str(row[b_col]).strip()
    d = float(row[distance_col])

    if A == "MPA":
        category = "Nearby nonorthosteric placements"
    else:
        category = "Orthosteric A + distinct B"

    results.append({
        "Combination": f"{A} → {B}",
        "A_Group": a_groups[A],
        "B_Group": b_groups[(A, B)],
        "Minimum_A_B_distance_A": round(d, 3),
        "Geometric_category": category
    })

final = pd.DataFrame(results)

final.to_csv(
    "/content/final_12_combination_binding_modes.csv",
    index=False
)

print("\n=== 12-COMBINATION BINDING MODE CLASSIFICATION ===")
print(final.to_string(index=False))

print("\nSaved: /content/final_12_combination_binding_modes.csv")

import pandas as pd

spatial = pd.read_csv(
    "/content/final_ligand_spatial_distances.csv"
)

# Use the actual column names from your CSV
spatial.columns = spatial.columns.str.strip()

a_groups = {
    "EE": 5, "LNG": 5, "NGMN": 5, "MPA": 2
}

b_groups = {
    ("EE", "LNG"): 3,
    ("EE", "NGMN"): 3,
    ("EE", "MPA"): 3,
    ("LNG", "EE"): 4,
    ("LNG", "NGMN"): 3,
    ("LNG", "MPA"): 3,
    ("NGMN", "EE"): 4,
    ("NGMN", "LNG"): 3,
    ("NGMN", "MPA"): 3,
    ("MPA", "EE"): 1,
    ("MPA", "LNG"): 1,
    ("MPA", "NGMN"): 1
}

records = []

for _, row in spatial.iterrows():
    A = str(row["Ligand A"]).strip()
    B = str(row["Ligand B"]).strip()

    distance = float(row["Minimum A–B distance (Å)"])
    centroid = float(row["Centroid distance (Å)"])

    if A == "MPA":
        category = "Nearby nonorthosteric placements"
    else:
        category = "Orthosteric A + distinct B"

    records.append({
        "Combination": f"{A} → {B}",
        "Selected_mode": int(row["Selected mode"]),
        "A_Group": a_groups[A],
        "B_Group": b_groups[(A, B)],
        "Centroid_distance_A": round(centroid, 3),
        "Minimum_A_B_distance_A": round(distance, 3),
        "Geometric_category": category
    })

final = pd.DataFrame(records)

assert len(final) == 12, "Expected exactly 12 combinations."

final.to_csv(
    "/content/final_12_combination_binding_modes.csv",
    index=False
)

print("=== 12-COMBINATION BINDING MODE CLASSIFICATION ===")
print(final.to_string(index=False))

print("\nSaved: /content/final_12_combination_binding_modes.csv")

import pandas as pd

spatial = pd.read_csv(
    "/content/final_ligand_spatial_distances.csv"
)

# Use the verified CSV column order.
# 0 = Ligand A
# 1 = Ligand B
# 2 = Selected mode
# 3 = Centroid distance
# 4 = Minimum A-B distance
spatial = spatial.iloc[:, :5].copy()
spatial.columns = [
    "Ligand_A",
    "Ligand_B",
    "Selected_mode",
    "Centroid_distance_A",
    "Minimum_AB_distance_A"
]

a_groups = {
    "EE": 5, "LNG": 5, "NGMN": 5, "MPA": 2
}

b_groups = {
    ("EE", "LNG"): 3,
    ("EE", "NGMN"): 3,
    ("EE", "MPA"): 3,
    ("LNG", "EE"): 4,
    ("LNG", "NGMN"): 3,
    ("LNG", "MPA"): 3,
    ("NGMN", "EE"): 4,
    ("NGMN", "LNG"): 3,
    ("NGMN", "MPA"): 3,
    ("MPA", "EE"): 1,
    ("MPA", "LNG"): 1,
    ("MPA", "NGMN"): 1
}

records = []

for _, row in spatial.iterrows():
    A = str(row["Ligand_A"]).strip()
    B = str(row["Ligand_B"]).strip()

    category = (
```

### Source lines 7452–7584
```python
        "Nearby nonorthosteric placements"
        if A == "MPA"
        else "Orthosteric A + spatially distinct B"
    )

    records.append({
        "Combination": f"{A} → {B}",
        "Selected_mode": int(row["Selected_mode"]),
        "A_Group": a_groups[A],
        "B_Group": b_groups[(A, B)],
        "Centroid_distance_A": round(
            float(row["Centroid_distance_A"]), 3
        ),
        "Minimum_AB_distance_A": round(
            float(row["Minimum_AB_distance_A"]), 3
        ),
        "Geometric_category": category
    })

final = pd.DataFrame(records)

assert len(final) == 12
assert final["Combination"].nunique() == 12

final.to_csv(
    "/content/final_12_combination_binding_modes.csv",
    index=False
)

print("=== FINAL 12-COMBINATION CLASSIFICATION ===")
print(final.to_string(index=False))

print("\nSUCCESS: All 12 combinations classified.")
print("Saved: /content/final_12_combination_binding_modes.csv")

import os
import glob

print("=== AVAILABLE VINA CONFIGURATION FILES ===")

patterns = [
    "/content/*config*",
    "/content/*vina*",
    "/content/*.txt",
    "/content/*.log"
]

found = sorted(set(
    path
    for pattern in patterns
    for path in glob.glob(pattern)
    if os.path.isfile(path)
))

for path in found:
    print("\nFILE:", path)

    if os.path.getsize(path) < 20000:
        with open(path, errors="ignore") as f:
            for line in f.readlines()[:35]:
                print(line.rstrip())

print("\n=== VINA SETTINGS IN NOTEBOOK HISTORY ===")

ip = get_ipython()
history = list(ip.history_manager.get_range())

keywords = [
    "Vina(",
    "compute_vina_maps",
    "dock(",
    "set_receptor(",
    "set_ligand_from_file("
]

matches = []

for _, _, source in history:
    if any(word in source for word in keywords):
        matches.append(source)

for i, source in enumerate(matches[-8:], 1):
    print(f"\n--- Relevant cell {i} ---")
    print(source[:3500])

if not matches:
    print("No Vina configuration found in notebook history.")

import os
import shutil
import subprocess

print("=== VINA INSTALLATION ===")

vina_path = shutil.which("vina")
print("Vina executable:", vina_path or "Not found")

if vina_path:
    result = subprocess.run(
        [vina_path, "--version"],
        capture_output=True,
        text=True
    )
    print("Version:", (result.stdout or result.stderr).strip())

print("\n=== EXISTING DOCKING OUTPUT HEADERS ===")

files = [
    "/content/LNG_docked_ER_EE.pdbqt",
    "/content/EE_docked_ER_LNG.pdbqt",
    "/content/NGMN_docked_ER_MPA.pdbqt",
    "/content/ER_EE_receptor.pdbqt",
    "/content/EE_ligand_B.pdbqt"
]

for path in files:
    print("\nFILE:", os.path.basename(path))

    if not os.path.exists(path):
        print("Missing")
        continue

    with open(path, errors="ignore") as f:
        for _, line in zip(range(18), f):
            print(line.rstrip())

print("\n=== VINA PYTHON PACKAGE ===")

try:
    import vina
    print("Vina Python package available")
    print("Package location:", vina.__file__)
except ImportError:
```

### Source lines 8101–8262
```python

            atoms[key] = np.array([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ])

    return atoms

print("=== ORIGINAL VS MINIMIZED PROTEIN GEOMETRY ===")

for A, original_path in original_files.items():
    minimized_path = f"/content/ER_{A}_minimized_complex.pdb"

    original = protein_heavy_coordinates(original_path)
    minimized = protein_heavy_coordinates(minimized_path)

    common = sorted(set(original) & set(minimized))

    distances = np.array([
        np.linalg.norm(original[k] - minimized[k])
        for k in common
    ])

    rmsd = np.sqrt(np.mean(distances**2))

    print(f"\nER–{A}")
    print("Matched heavy atoms:", len(common))
    print("Direct coordinate RMSD (Å):", round(rmsd, 4))
    print("Mean displacement (Å):", round(distances.mean(), 4))
    print("Maximum displacement (Å):", round(distances.max(), 4))
    print("Atoms displaced > 0.5 Å:", int(np.sum(distances > 0.5)))
    print("Atoms displaced > 1.0 Å:", int(np.sum(distances > 1.0)))

print("\nNo files modified.")

import numpy as np

original_files = {
    "EE": "/content/estrogenreceptor_ee_out_1.-9.9.complex.pdb",
    "LNG": "/content/estrogenreceptor_lng_out_1.-9.4.complex.pdb",
    "NGMN": "/content/estrogenreceptor_ngmn_out_1.-10.4.complex.pdb",
    "MPA": "/content/estrogenreceptor_mpa_out_4.-7.0.complex.pdb"
}

def ligand_heavy_xyz(path):
    coords = []

    with open(path, errors="ignore") as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue

            element = line[76:78].strip().upper()

            if not element:
                element = line[12:16].strip()[0].upper()

            if element in ("H", "D"):
                continue

            coords.append([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ])

    return np.asarray(coords, dtype=float)

print("=== LIGAND A GEOMETRY COMPARISON ===")

for A, original_path in original_files.items():
    minimized_path = f"/content/ER_{A}_minimized_complex.pdb"

    original = ligand_heavy_xyz(original_path)
    minimized = ligand_heavy_xyz(minimized_path)

    print(f"\nER–{A}")
    print("Original heavy atoms:", len(original))
    print("Minimized heavy atoms:", len(minimized))

    if len(original) != len(minimized):
        print("Atom counts differ; comparison skipped.")
        continue

    original_centroid = original.mean(axis=0)
    minimized_centroid = minimized.mean(axis=0)

    centroid_shift = np.linalg.norm(
        original_centroid - minimized_centroid
    )

    print("Original centroid:", np.round(original_centroid, 3))
    print("Minimized centroid:", np.round(minimized_centroid, 3))
    print("Centroid displacement (Å):", round(centroid_shift, 4))

    print(
        "WARNING: Atom-by-atom RMSD not calculated;",
        "atom ordering has not been verified."
    )

print("\nNo files modified.")

from pathlib import Path
from collections import Counter

original_files = {
    "EE": "/content/estrogenreceptor_ee_out_1.-9.9.complex.pdb",
    "LNG": "/content/estrogenreceptor_lng_out_1.-9.4.complex.pdb",
    "NGMN": "/content/estrogenreceptor_ngmn_out_1.-10.4.complex.pdb",
    "MPA": "/content/estrogenreceptor_mpa_out_4.-7.0.complex.pdb"
}

def heavy_atom_key(line):
    return (
        line[22:26].strip(),
        line[17:20].strip(),
        line[12:16].strip()
    )

def is_heavy_atom(line):
    element = line[76:78].strip().upper()

    if not element:
        element = line[12:16].strip()[0].upper()

    return element not in ("H", "D")

print("=== RECONSTRUCTED CONTROL INTERMEDIATES ===")

for A, original_path in original_files.items():
    minimized_path = f"/content/ER_{A}_minimized_complex.pdb"
    output_path = f"/content/ER_{A}_reconstructed_control_heavy.pdb"

    with open(original_path) as f:
        original_lines = [
            line for line in f
            if line.startswith(("ATOM  ", "HETATM"))
            and is_heavy_atom(line)
        ]

    original_protein_keys = {
        heavy_atom_key(line)
        for line in original_lines
        if line.startswith("ATOM  ")
    }

    with open(minimized_path) as f:
        added_lines = [
            line for line in f
            if line.startswith("ATOM  ")
            and is_heavy_atom(line)
            and heavy_atom_key(line) not in original_protein_keys
        ]

    protein_lines = [
        line for line in original_lines
        if line.startswith("ATOM  ")
    ] + added_lines

    ligand_lines = [
        line for line in original_lines
```

### Source lines 9774–9996
```python
        else:
            ligand.append(xyz)

    return protein, np.array(ligand)

def kabsch_transform(moving, reference):
    moving_center = moving.mean(axis=0)
    reference_center = reference.mean(axis=0)

    p = moving - moving_center
    q = reference - reference_center

    u, s, vt = np.linalg.svd(p.T @ q)

    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)

    rotation = u @ correction @ vt

    return rotation, moving_center, reference_center

def transform(coords, rotation, moving_center, reference_center):
    return (coords - moving_center) @ rotation + reference_center

def rmsd(a, b):
    return float(np.sqrt(np.mean(np.sum((a - b)**2, axis=1))))

results = []

with zipfile.ZipFile(backup) as z:
    names = {Path(n).name: n for n in z.namelist()}

    for lig in ligands:
        original_file = (
            f"ER_{lig}_reconstructed_control_corrected.pdb"
        )
        minimized_file = f"ER_{lig}_minimized_complex.pdb"

        original_text = z.read(
            names[original_file]
        ).decode("utf-8", errors="replace")

        minimized_text = z.read(
            names[minimized_file]
        ).decode("utf-8", errors="replace")

        original_protein, original_ligand = read_structure(
            original_text
        )
        minimized_protein, minimized_ligand = read_structure(
            minimized_text
        )

        keys = sorted(original_protein)

        if set(keys) != set(minimized_protein):
            raise ValueError(f"Protein atom mismatch: {lig}")

        if len(original_ligand) != len(minimized_ligand):
            raise ValueError(f"Ligand count mismatch: {lig}")

        ca_keys = [k for k in keys if k[2] == "CA"]

        original_ca = np.array([
            original_protein[k] for k in ca_keys
        ])
        minimized_ca = np.array([
            minimized_protein[k] for k in ca_keys
        ])

        # Globally align minimized receptor to original receptor.
        rotation, mc, rc = kabsch_transform(
            minimized_ca, original_ca
        )

        aligned_ca = transform(
            minimized_ca, rotation, mc, rc
        )

        whole_ca_rmsd = rmsd(aligned_ca, original_ca)

        # Define original Ligand A pocket:
        # residues with any protein heavy atom within 5 Å.
        pocket_residues = set()

        for key in keys:
            xyz = original_protein[key]

            distances = np.linalg.norm(
                original_ligand - xyz,
                axis=1
            )

            if distances.min() <= 5.0:
                pocket_residues.add((key[0], key[1]))

        pocket_ca_keys = [
            k for k in ca_keys
            if (k[0], k[1]) in pocket_residues
        ]

        pocket_heavy_keys = [
            k for k in keys
            if (k[0], k[1]) in pocket_residues
        ]

        def subset_rmsd(subset_keys):
            if not subset_keys:
                return float("nan")

            original_xyz = np.array([
                original_protein[k] for k in subset_keys
            ])

            minimized_xyz = np.array([
                minimized_protein[k] for k in subset_keys
            ])

            aligned_xyz = transform(
                minimized_xyz, rotation, mc, rc
            )

            return rmsd(aligned_xyz, original_xyz)

        pocket_ca_rmsd = subset_rmsd(pocket_ca_keys)
        pocket_heavy_rmsd = subset_rmsd(pocket_heavy_keys)

        # Centroid shift does not require ligand atom ordering.
        aligned_ligand = transform(
            minimized_ligand, rotation, mc, rc
        )

        centroid_shift = float(np.linalg.norm(
            aligned_ligand.mean(axis=0)
            - original_ligand.mean(axis=0)
        ))

        results.append({
            "Ligand A": lig,
            "Whole-receptor CA RMSD (Å)": whole_ca_rmsd,
            "Pocket CA RMSD (Å)": pocket_ca_rmsd,
            "Pocket heavy-atom RMSD (Å)": pocket_heavy_rmsd,
            "Ligand A centroid shift (Å)": centroid_shift,
            "Pocket residues": len(pocket_residues),
            "Matched pocket CA": len(pocket_ca_keys),
            "Matched pocket heavy atoms": len(pocket_heavy_keys)
        })

df = pd.DataFrame(results)

print("=" * 95)
print("TABLE 4 — ERα STRUCTURAL DISPLACEMENT")
print("=" * 95)
print(df.round(4).to_string(index=False))

print("\nMETHOD:")
print("Alignment: global receptor Cα Kabsch fit")
print("Pocket: original Ligand A protein contacts within 5 Å")
print("RMSD: after global protein alignment")
print("Ligand displacement: heavy-atom centroid shift")

print("\nSTEP 7 COMPLETE — ORIGINAL FILES UNCHANGED")

# ERα — STEP 8: INSPECT MINIMIZATION METADATA
# Read-only: searches saved CSV, JSON, YAML and logs.

import zipfile
from pathlib import Path

backup = Path(
    "/content/drive/MyDrive/Molecular_Docking_Research/"
    "ERalpha_Backup_20261026_190013.zip"
)

terms = [
    "amber", "openff", "sage", "forcefield",
    "force field", "nocutoff", "cutoff",
    "minimiz", "tolerance", "iteration",
    "openmm", "solvent", "implicit", "vacuum",
    "energy", "platform"
]

with zipfile.ZipFile(backup) as z:
    candidates = [
        n for n in z.namelist()
        if Path(n).suffix.lower() in {
            ".csv", ".json", ".yaml", ".yml",
            ".log", ".txt", ".py"
        }
        and not n.startswith(("sample_data/", ".config/"))
    ]

    print("METADATA FILES:", len(candidates))

    for name in candidates:
        content = z.read(name).decode(
            "utf-8", errors="replace"
        )

        lines = content.splitlines()
        matches = [
            (i + 1, line)
            for i, line in enumerate(lines)
            if any(term in line.lower() for term in terms)
        ]

        if matches:
            print("\n" + "=" * 75)
            print("FILE:", name)
            print("=" * 75)

            for line_no, line in matches[:30]:
                print(f"Line {line_no}: {line[:250]}")

            if len(matches) > 30:
                print(
                    f"... {len(matches)-30} additional matches"
                )

print("\nSTEP 8 COMPLETE — NO FILES MODIFIED")

# ERα — STEP 9: EXPERIMENTAL ORTHOSTERIC POCKET CONTACTS
# Read-only: uses saved receptor contacts and experimental PDB.
```


---

# Androgen receptor (AR)

# AR — original Dock–EM–Dock code excerpts
**Source:** `ar_energyminimization.py` (6,598 lines).
**Provenance:** All code excerpts below are copied verbatim from the uploaded Colab-exported Python source, with original line ranges. Sections are excerpts, not a standalone runnable script. Original complete file is included in the ZIP.
**Scope:** The original notebook may contain debugging, corrections, alternate parameterizations, or non-executable Colab shell commands. These excerpts do not imply all code paths were executed successfully.

## Original code excerpts
### Source lines 224–312
```python
for ligand, path in files.items():
    assert Path(path).exists(), f"Missing file: {path}"

    ligand_xyz = get_heavy_atoms(path, ligand_only=True)

    if len(ligand_xyz) == 0:
        print(f"\n{ligand}: No HETATM ligand detected")
        continue

    ligand_center = ligand_xyz.mean(axis=0)

    centroid_distance = np.linalg.norm(
        ligand_center - dht_center
    )

    minimum_distance = np.min(
        np.linalg.norm(
            ligand_xyz[:, None, :] - dht_xyz[None, :, :],
            axis=2
        )
    )

    print(f"\nAR–{ligand}")
    print("Ligand heavy atoms:", len(ligand_xyz))
    print("Centroid distance from DHT:",
          f"{centroid_distance:.3f} Å")
    print("Minimum distance from DHT:",
          f"{minimum_distance:.3f} Å")

print("\nAnalysis complete. No files modified.")

from pathlib import Path

print("=== PDB FILES IN GOOGLE COLAB ===")

pdb_files = sorted(Path("/content").glob("*.pdb"))

if not pdb_files:
    print("No PDB files found in /content/")
else:
    for i, file in enumerate(pdb_files, 1):
        print(f"{i}. {file.name}")

print("\nTotal PDB files:", len(pdb_files))
print("No files modified.")

import numpy as np
from pathlib import Path

files = {
    "EE": "/content/AR_clean_receptor_ee_out_1.-9.2.complex.pdb",
    "LNG": "/content/AR_clean_receptor_lng_out_1.-10.5.complex.pdb",
    "NGMN": "/content/AR_clean_receptor_ngmn_out_1.-10.5.complex.pdb",
    "MPA": "/content/AR_clean_receptor_mpa_out_1.-12.1.complex.pdb"
}

dht_file = "/content/AR_DHT_reference.pdb"

def read_atoms(path, record_types):
    atoms = []

    with open(path) as f:
        for line in f:
            if line[:6].strip() not in record_types:
                continue

            element = line[76:78].strip().upper()
            if element in ("H", "D"):
                continue

            try:
                xyz = np.array([
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54])
                ])
            except ValueError:
                continue

            atoms.append((line, xyz))

    return atoms

assert Path(dht_file).exists(), "DHT reference missing."

dht_atoms = read_atoms(dht_file, {"HETATM", "ATOM"})
dht_xyz = np.array([xyz for _, xyz in dht_atoms])

assert len(dht_xyz) == 21, (
```

### Source lines 330–418
```python
    residue_counts = {}
    for line, xyz in atoms:
        name = line[17:20].strip()
        residue_counts[name] = residue_counts.get(name, 0) + 1

    ligand_xyz = np.array([xyz for _, xyz in atoms])

    ligand_center = ligand_xyz.mean(axis=0)

    centroid_distance = np.linalg.norm(
        ligand_center - dht_center
    )

    minimum_distance = np.min(
        np.linalg.norm(
            ligand_xyz[:, None, :] - dht_xyz[None, :, :],
            axis=2
        )
    )

    print(f"\nAR–{ligand}")
    print("HETATM residue types:", residue_counts)
    print("Ligand heavy atoms:", len(ligand_xyz))
    print(
        "Centroid distance from DHT:",
        f"{centroid_distance:.3f} Å"
    )
    print(
        "Minimum distance from DHT:",
        f"{minimum_distance:.3f} Å"
    )

print("\nNo files modified.")

from pathlib import Path
from collections import Counter

complexes = {
    "EE": "AR_clean_receptor_ee_out_1.-9.2.complex.pdb",
    "LNG": "AR_clean_receptor_lng_out_1.-10.5.complex.pdb",
    "NGMN": "AR_clean_receptor_ngmn_out_1.-10.5.complex.pdb",
    "MPA": "AR_clean_receptor_mpa_out_1.-12.1.complex.pdb"
}

expected_ligand_atoms = {
    "EE": 22,
    "LNG": 23,
    "NGMN": 24,
    "MPA": 28
}

print("=== AR COMPLEX PRE-MINIMIZATION CHECK ===")

for ligand, filename in complexes.items():
    path = Path("/content") / filename

    if not path.exists():
        print(f"\nAR–{ligand}: FILE MISSING")
        continue

    with open(path) as f:
        lines = f.readlines()

    protein = [x for x in lines if x.startswith("ATOM  ")]
    hetero = [x for x in lines if x.startswith("HETATM")]

    residues = {
        (x[21], x[22:26], x[26])
        for x in protein
    }

    ligand_residues = Counter(
        x[17:20].strip() for x in hetero
    )

    print(f"\nAR–{ligand}")
    print("Protein atoms:", len(protein))
    print("Protein residues:", len(residues))
    print("Ligand records:", dict(ligand_residues))
    print(
        "Expected ligand heavy atoms:",
        expected_ligand_atoms[ligand]
    )

    status = (
        len(residues) == 249
        and len(hetero) == expected_ligand_atoms[ligand]
        and set(ligand_residues) == {"UNL"}
    )
```

### Source lines 954–1253
```python
    print(
        "PASS" if Path(file).is_file() else "MISSING",
        Path(file).name
    )

from openmm import app
from pathlib import Path

print("=== AR PROTEIN PREPARATION ===")

receptor_path = "/content/AR_clean_receptor.pdb"
output_path = "/content/AR_prepared/AR_protein_hydrogenated.pdb"

# Load the cleaned AR receptor
pdb = app.PDBFile(receptor_path)

# AMBER14 protein parameters
forcefield = app.ForceField("amber14/protein.ff14SB.xml")

# Add missing hydrogen atoms at approximately physiological pH
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield, pH=7.4)

# Confirm that AMBER14 can parameterize the protein
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=None
)

# Save hydrogenated receptor
with open(output_path, "w") as f:
    app.PDBFile.writeFile(
        modeller.topology,
        modeller.positions,
        f,
        keepIds=True
    )

print("AMBER14 protein parameterization: PASS")
print("Protein atoms before:", pdb.topology.getNumAtoms())
print("Protein atoms after:", modeller.topology.getNumAtoms())
print("Protein residues:", modeller.topology.getNumResidues())
print("Saved:", output_path)
print("\nOriginal receptor and ligand files unchanged.")

from openmm import app
from pathlib import Path

print("=== AR PROTEIN PREPARATION — TERMINAL FIX ===")

pdb = app.PDBFile("/content/AR_clean_receptor.pdb")

forcefield = app.ForceField(
    "amber14/protein.ff14SB.xml"
)

modeller = app.Modeller(
    pdb.topology,
    pdb.positions
)

# Add hydrogens without force-field template matching
modeller.addHydrogens(pH=7.4)

print("Hydrogens added successfully.")
print("Protein atoms:", modeller.topology.getNumAtoms())

# Test AMBER14 compatibility
try:
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None
    )

    output = "/content/AR_prepared/AR_protein_hydrogenated.pdb"

    with open(output, "w") as f:
        app.PDBFile.writeFile(
            modeller.topology,
            modeller.positions,
            f,
            keepIds=True
        )

    print("AMBER14 parameterization: PASS")
    print("Saved:", output)

except Exception as e:
    print("AMBER14 parameterization: FAILED")
    print("Reason:", str(e))

print("\nOriginal files unchanged.")

# Commented out IPython magic to ensure Python compatibility.

# %pip -q install pdbfixer

from pdbfixer import PDBFixer
from openmm import app
from pathlib import Path

print("=== AR RECEPTOR TERMINAL REPAIR ===")

fixer = PDBFixer(
    filename="/content/AR_clean_receptor.pdb"
)

fixer.findMissingResidues()

# Do not reconstruct missing sequence regions.
fixer.missingResidues = {}

fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()

fixer.findMissingAtoms()

print("Missing heavy atoms by residue:")
for residue, atoms in fixer.missingAtoms.items():
    print(
        f"{residue.name} {residue.id}:",
        [atom.name for atom in atoms]
    )

print("Missing terminal atoms:")
for residue, atoms in fixer.missingTerminals.items():
    print(
        f"{residue.name} {residue.id}:",
        atoms
    )

fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.4)

forcefield = app.ForceField(
    "amber14/protein.ff14SB.xml"
)

try:
    system = forcefield.createSystem(
        fixer.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None
    )

    output = "/content/AR_prepared/AR_protein_hydrogenated.pdb"

    with open(output, "w") as f:
        app.PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            f,
            keepIds=True
        )

    print("\nAMBER14 parameterization: PASS")
    print("Saved:", output)

except Exception as e:
    print("\nAMBER14 parameterization: FAILED")
    print("Reason:", str(e))

print("\nOriginal receptor and ligand files unchanged.")

from openmm import app
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Molecule
from rdkit import Chem
from pathlib import Path

print("=== AR COMPLEX PARAMETERIZATION ===")

protein = app.PDBFile(
    "/content/AR_prepared/AR_protein_hydrogenated.pdb"
)

ligand_dir = Path("/content/AR_prepared")

for name in ["EE", "LNG", "NGMN", "MPA"]:
    print(f"\nAR–{name}")

    try:
        sdf_path = ligand_dir / f"AR_{name}_docked.sdf"

        # Read the docked ligand, including hydrogens
        rdkit_mol = Chem.SDMolSupplier(
            str(sdf_path),
            removeHs=False
        )[0]

        if rdkit_mol is None:
            raise ValueError("Unable to read ligand SDF")

        ligand = Molecule.from_rdkit(
            rdkit_mol,
            allow_undefined_stereo=False
        )

        # Create combined protein–ligand topology
        ligand_topology = ligand.to_topology().to_openmm()

        modeller = app.Modeller(
            protein.topology,
            protein.positions
        )

        from openff.units import unit as off_unit
        from openmm import unit as omm_unit

        coords = ligand.conformers[0].m_as(
            off_unit.nanometer
        )

        ligand_positions = [
            app.Vec3(*xyz) * omm_unit.nanometer
            for xyz in coords
        ]

        modeller.add(
            ligand_topology,
            ligand_positions
        )

        # AMBER14 protein + OpenFF 2.2.1 ligand
        forcefield = app.ForceField(
            "amber14/protein.ff14SB.xml"
        )

        generator = SMIRNOFFTemplateGenerator(
            molecules=[ligand],
            forcefield="openff-2.2.1"
        )

        forcefield.registerTemplateGenerator(
            generator.generator
        )

        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=None
        )

        print("Parameterization: PASS")
        print("Total atoms:", system.getNumParticles())
        print("Ligand heavy atoms:", ligand.n_atoms - sum(
            atom.symbol == "H" for atom in ligand.atoms
        ))

    except Exception as e:
        print("Parameterization: FAILED")
        print("Reason:", type(e).__name__, str(e)[:300])

print("\nNo structures modified.")

from pathlib import Path
from rdkit import Chem
from openff.toolkit import Molecule
from openmm import app, Vec3, unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.units import unit as off_unit

print("=== CORRECTED AR COMPLEX PARAMETERIZATION ===")

protein = app.PDBFile(
    "/content/AR_prepared/AR_protein_hydrogenated.pdb"
)

for name in ["EE", "LNG", "NGMN", "MPA"]:
    print(f"\nAR–{name}")

    try:
        path = Path(f"/content/AR_prepared/AR_{name}_docked.sdf")

        mol = Chem.SDMolSupplier(
            str(path), removeHs=False
        )[0]

        if mol is None:
            raise ValueError("Could not read prepared ligand")

        # Inspect stereochemistry; do not guess missing configurations.
        stereo = Chem.FindPotentialStereo(mol)
        unspecified = [
            str(s) for s in stereo
            if s.specified == Chem.StereoSpecified.Unspecified
        ]

        if unspecified:
            print("Unspecified stereochemistry:", unspecified)
            print("REVIEW REQUIRED — no configuration assigned")
            continue

        ligand = Molecule.from_rdkit(
            mol,
            allow_undefined_stereo=False
        )

```

### Source lines 1254–1349
```python
        coords = ligand.conformers[0].m_as(
            off_unit.nanometer
        )

        positions = [
            Vec3(float(x), float(y), float(z)) * unit.nanometer
            for x, y, z in coords
        ]

        modeller = app.Modeller(
            protein.topology, protein.positions
        )

        modeller.add(
            ligand.to_topology().to_openmm(),
            positions
        )

        forcefield = app.ForceField(
            "amber14/protein.ff14SB.xml"
        )

        generator = SMIRNOFFTemplateGenerator(
            molecules=[ligand],
            forcefield="openff-2.2.1"
        )

        forcefield.registerTemplateGenerator(
            generator.generator
        )

        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=None
        )

        print("Parameterization: PASS")
        print("Total atoms:", system.getNumParticles())

    except Exception as e:
        print("Parameterization: FAILED")
        print("Reason:", type(e).__name__, str(e)[:400])

print("\nNo files modified.")

from rdkit import Chem

print("=== NGMN STEREOCHEMISTRY CHECK ===")

paths = {
    "PubChem reference": "/content/AR_ligands/NGMN.sdf",
    "Prepared docked": "/content/AR_prepared/AR_NGMN_docked.sdf"
}

for label, path in paths.items():
    print(f"\n{label}")

    mol = Chem.SDMolSupplier(
        path, removeHs=False
    )[0]

    if mol is None:
        print("FAILED: Cannot read structure")
        continue

    Chem.AssignStereochemistry(
        mol, cleanIt=True, force=True
    )

    print("Heavy atoms:", mol.GetNumHeavyAtoms())
    print("Isomeric SMILES:", Chem.MolToSmiles(
        Chem.RemoveHs(mol), isomericSmiles=True
    ))

    for info in Chem.FindPotentialStereo(mol):
        print(
            "Stereo type:", info.type,
            "| Center:", info.centeredOn,
            "| Specified:", info.specified,
            "| Descriptor:", info.descriptor
        )

print("\nNo files modified.")

from rdkit import Chem
from openff.toolkit import Molecule
from openff.units import unit as off_unit
from openmm import app, Vec3, unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

print("=== AR–NGMN PARAMETERIZATION ===")

try:
    path = "/content/AR_prepared/AR_NGMN_docked.sdf"

```

### Source lines 1367–1666
```python

    coords = ligand.conformers[0].m_as(
        off_unit.nanometer
    )

    ligand_positions = [
        Vec3(float(x), float(y), float(z)) * unit.nanometer
        for x, y, z in coords
    ]

    modeller = app.Modeller(
        protein.topology,
        protein.positions
    )

    modeller.add(
        ligand.to_topology().to_openmm(),
        ligand_positions
    )

    forcefield = app.ForceField(
        "amber14/protein.ff14SB.xml"
    )

    generator = SMIRNOFFTemplateGenerator(
        molecules=[ligand],
        forcefield="openff-2.2.1"
    )

    forcefield.registerTemplateGenerator(
        generator.generator
    )

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None
    )

    print("Parameterization: PASS")
    print("Total atoms:", system.getNumParticles())
    print("NGMN atoms:", ligand.n_atoms)
    print("Oxime E/Z: unspecified")
    print("Original docked coordinates retained")

except Exception as e:
    print("Parameterization: FAILED")
    print("Reason:", type(e).__name__, str(e)[:600])

print("\nNo files modified.")

import numpy as np
from pathlib import Path
from rdkit import Chem

from openmm import (
    app, unit, Vec3, Platform,
    VerletIntegrator, Context,
    LocalEnergyMinimizer, CustomExternalForce
)
from openff.toolkit import Molecule
from openff.units import unit as off_unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

print("=== AR RESTRAINED ENERGY MINIMIZATION ===")

outdir = Path("/content/AR_minimized")
outdir.mkdir(exist_ok=True)

protein = app.PDBFile(
    "/content/AR_prepared/AR_protein_hydrogenated.pdb"
)

protein_atom_count = protein.topology.getNumAtoms()

for name in ["EE", "LNG", "NGMN", "MPA"]:
    print(f"\n--- AR–{name} ---")

    try:
        sdf_path = f"/content/AR_prepared/AR_{name}_docked.sdf"

        rdkit_mol = Chem.SDMolSupplier(
            sdf_path, removeHs=False
        )[0]

        if rdkit_mol is None:
            raise ValueError("Ligand SDF could not be read")

        ligand = Molecule.from_rdkit(
            rdkit_mol,
            allow_undefined_stereo=(name == "NGMN")
        )

        coords = ligand.conformers[0].m_as(
            off_unit.nanometer
        )

        ligand_positions = [
            Vec3(float(x), float(y), float(z)) * unit.nanometer
            for x, y, z in coords
        ]

        modeller = app.Modeller(
            protein.topology,
            protein.positions
        )

        modeller.add(
            ligand.to_topology().to_openmm(),
            ligand_positions
        )

        forcefield = app.ForceField(
            "amber14/protein.ff14SB.xml"
        )

        generator = SMIRNOFFTemplateGenerator(
            molecules=[ligand],
            forcefield="openff-2.2.1"
        )

        forcefield.registerTemplateGenerator(
            generator.generator
        )

        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=None
        )

        # Restrain protein heavy atoms near their starting positions.
        restraint = CustomExternalForce(
            "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
        )

        restraint.addGlobalParameter(
            "k", 1000.0
        )

        for parameter in ["x0", "y0", "z0"]:
            restraint.addPerParticleParameter(parameter)

        initial_nm = modeller.positions.value_in_unit(
            unit.nanometer
        )

        for atom in modeller.topology.atoms():
            if (
                atom.index < protein_atom_count
                and atom.element.symbol != "H"
            ):
                xyz = initial_nm[atom.index]

                restraint.addParticle(
                    atom.index,
                    [xyz.x, xyz.y, xyz.z]
                )

        system.addForce(restraint)

        integrator = VerletIntegrator(
            0.001 * unit.picoseconds
        )

        context = Context(
            system,
            integrator,
            Platform.getPlatformByName("CPU")
        )

        context.setPositions(modeller.positions)

        initial_energy = context.getState(
            getEnergy=True
        ).getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        )

        LocalEnergyMinimizer.minimize(
            context,
            tolerance=10,
            maxIterations=1000
        )

        state = context.getState(
            getPositions=True,
            getEnergy=True
        )

        final_energy = state.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        )

        final_positions = state.getPositions(
            asNumpy=True
        )

        # Save minimized protein–ligand complex.
        complex_path = outdir / f"AR_{name}_minimized_complex.pdb"

        with open(complex_path, "w") as f:
            app.PDBFile.writeFile(
                modeller.topology,
                final_positions,
                f,
                keepIds=True
            )

        # Save minimized ligand with original chemical connectivity.
        ligand_nm = final_positions.value_in_unit(
            unit.nanometer
        )[protein_atom_count:]

        conf = rdkit_mol.GetConformer()

        for i, xyz in enumerate(ligand_nm):
            conf.SetAtomPosition(
                i,
                (
                    float(xyz[0]) * 10,
                    float(xyz[1]) * 10,
                    float(xyz[2]) * 10
                )
            )

        ligand_path = outdir / f"AR_{name}_minimized_ligand.sdf"

        writer = Chem.SDWriter(str(ligand_path))
        writer.write(rdkit_mol)
        writer.close()

        print("Minimization: PASS")
        print(f"Initial energy: {initial_energy:.2f} kJ/mol")
        print(f"Final energy:   {final_energy:.2f} kJ/mol")
        print("Saved:", complex_path.name)
        print("Saved:", ligand_path.name)

        del context, integrator

    except Exception as e:
        print("Minimization: FAILED")
        print("Reason:", type(e).__name__, str(e)[:500])

print("\n=== MINIMIZATION BATCH COMPLETE ===")

# Commented out IPython magic to ensure Python compatibility.

# %pip -q install vina meeko

from pathlib import Path
import importlib.metadata

print("=== AR SEQUENTIAL DOCKING PREPARATION ===")

for package in ["vina", "meeko"]:
    try:
        version = importlib.metadata.version(package)
        print(f"{package}: PASS | version {version}")
    except Exception as e:
        print(f"{package}: FAILED | {e}")

print("\n=== MINIMIZED AR COMPLEXES ===")

folder = Path("/content/AR_minimized")

for name in ["EE", "LNG", "NGMN", "MPA"]:
    complex_path = folder / f"AR_{name}_minimized_complex.pdb"
    ligand_path = folder / f"AR_{name}_minimized_ligand.sdf"

    complex_ok = complex_path.is_file()
    ligand_ok = ligand_path.is_file()

    print(
        f"AR–{name}: "
        f"Complex {'PASS' if complex_ok else 'MISSING'} | "
        f"Ligand {'PASS' if ligand_ok else 'MISSING'}"
    )

print("\n=== SEQUENTIAL DOCKING PAIRS ===")

ligands = ["EE", "LNG", "NGMN", "MPA"]

pairs = [
    (a, b)
    for a in ligands
    for b in ligands
    if a != b
]

for a, b in pairs:
    print(f"{a} → {b}")

print(f"\nTotal ordered combinations: {len(pairs)}")

import subprocess
import importlib.metadata

print("=== INSTALLING AR DOCKING SOFTWARE ===")

```

### Source lines 1667–1710
```python
result = subprocess.run(
    [
        "mamba", "install", "-y",
        "-c", "conda-forge",
        "vina",
        "meeko"
    ],
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT
)

print(result.stdout[-5000:])

if result.returncode != 0:
    print("\nInstallation: FAILED")
else:
    print("\nInstallation command: PASS")

    print("\n=== SOFTWARE VERIFICATION ===")

    for package in ["vina", "meeko"]:
        try:
            version = importlib.metadata.version(package)
            print(f"{package}: PASS | version {version}")
        except Exception as e:
            print(f"{package}: NOT DETECTED | {e}")

    check = subprocess.run(
        ["vina", "--version"],
        text=True,
        capture_output=True
    )

    print(
        "Vina executable:",
        check.stdout.strip() or check.stderr.strip()
    )

from pathlib import Path
from rdkit import Chem
from openmm import app

print("=== AR SEQUENTIAL RECEPTOR PREPARATION ===")
```

### Source lines 2227–2483
```python

    lines = pdbqt.read_text().splitlines()

    atoms = [
        line for line in lines
        if line.startswith(("ATOM  ", "HETATM"))
    ]

    torsdof = [
        line for line in lines
        if line.startswith("TORSDOF")
    ]

    print("Preparation: PASS")
    print("PDBQT atom records:", len(atoms))
    print("Torsional degrees of freedom:",
          torsdof[0] if torsdof else "Not found")
    print("Saved:", pdbqt.name)

print("\n=== LIGAND B PREPARATION COMPLETE ===")

import numpy as np
from pathlib import Path

print("=== AR SEQUENTIAL DOCKING BOX ===")

reference = Path("/content/AR_DHT_reference.pdb")
output = Path("/content/AR_sequential")
output.mkdir(parents=True, exist_ok=True)

# Extract experimental DHT heavy-atom coordinates.
coords = []

for line in reference.read_text().splitlines():
    if line.startswith(("ATOM  ", "HETATM")):
        element = line[76:78].strip().upper()

        if element != "H":
            coords.append([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ])

coords = np.array(coords)

if len(coords) != 21:
    raise ValueError(
        f"Expected 21 DHT heavy atoms, found {len(coords)}"
    )

center = coords.mean(axis=0)

# Expanded box around the experimental binding site.
# Dimensions are in angstroms.
size = np.array([30.0, 30.0, 30.0])

config_path = output / "AR_docking_box.txt"

config_path.write_text(
    f"center_x = {center[0]:.3f}\n"
    f"center_y = {center[1]:.3f}\n"
    f"center_z = {center[2]:.3f}\n"
    f"size_x = {size[0]:.1f}\n"
    f"size_y = {size[1]:.1f}\n"
    f"size_z = {size[2]:.1f}\n"
)

print("Reference: Experimental DHT")
print("DHT heavy atoms:", len(coords))

print("\nDocking box center (Å):")
print(f"X = {center[0]:.3f}")
print(f"Y = {center[1]:.3f}")
print(f"Z = {center[2]:.3f}")

print("\nDocking box dimensions (Å):")
print(f"{size[0]:.1f} × {size[1]:.1f} × {size[2]:.1f}")

print("\nSaved:", config_path)
print("\n=== DOCKING BOX READY ===")

from pathlib import Path
import subprocess
import csv
import re

print("=== AR SEQUENTIAL AUTODOCK VINA DOCKING ===")

folder = Path("/content/AR_sequential")
results_dir = folder / "docking_results"
results_dir.mkdir(parents=True, exist_ok=True)

config = folder / "AR_docking_box.txt"

ligands = ["EE", "LNG", "NGMN", "MPA"]

pairs = [
    (a, b)
    for a in ligands
    for b in ligands
    if a != b
]

results = []

for index, (a, b) in enumerate(pairs, start=1):

    print(f"\n[{index}/12] {a} → {b}")

    receptor = folder / f"AR_{a}_fixed_receptor.pdbqt"
    ligand = folder / f"AR_{b}_ligandB.pdbqt"

    output = results_dir / f"AR_{a}_to_{b}_docked.pdbqt"
    log = results_dir / f"AR_{a}_to_{b}_vina.log"

    command = [
        "vina",
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--config", str(config),
        "--exhaustiveness", "16",
        "--num_modes", "10",
        "--energy_range", "5",
        "--seed", "2026",
        "--cpu", "2",
        "--out", str(output)
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True
    )

    log.write_text(result.stdout + "\n" + result.stderr)

    if result.returncode != 0:
        print("Docking: FAILED")
        print((result.stderr or result.stdout)[-600:])

        results.append([a, b, "FAILED", ""])
        continue

    # Read the best Vina score from the output poses.
    scores = []

    if output.exists():
        for line in output.read_text().splitlines():
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                scores.append(float(parts[3]))

    if scores:
        print("Docking: PASS")
        print(f"Best Vina score: {scores[0]:.3f} kcal/mol")
        print("Modes generated:", len(scores))

        results.append([
            a, b, "PASS", scores[0]
        ])
    else:
        print("Docking: NO SCORES FOUND")
        results.append([a, b, "NO SCORES", ""])

# Save results table.
csv_path = results_dir / "AR_sequential_docking_results.csv"

with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)

    writer.writerow([
        "Ligand_A",
        "Ligand_B",
        "Status",
        "Best_Vina_Score_kcal_mol"
    ])

    writer.writerows(results)

print("\n=== AR SEQUENTIAL DOCKING SUMMARY ===")

for a, b, status, score in results:
    print(f"{a:5s} → {b:5s} | {status:9s} | {score}")

print("\nSaved results:", csv_path)
print("Completed:", sum(r[2] == "PASS" for r in results), "/ 12")

from pathlib import Path
import numpy as np
import csv

print("=== AR SEQUENTIAL DOCKING SPATIAL ANALYSIS ===")

folder = Path("/content/AR_sequential")
results_dir = folder / "docking_results"
ligands = ["EE", "LNG", "NGMN", "MPA"]

def heavy_coords(lines):
    coords = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            atom_type = line.split()[-1].upper()
            if atom_type not in ("H", "HD"):
                coords.append([
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54])
                ])
    return np.array(coords, dtype=float)

def read_modes(path):
    modes = []
    current = []
    score = None

    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            current = []
            score = None
        elif line.startswith("REMARK VINA RESULT:"):
            score = float(line.split()[3])
        elif line.startswith(("ATOM  ", "HETATM")):
            current.append(line)
        elif line.startswith("ENDMDL"):
            if current and score is not None:
                modes.append((score, heavy_coords(current)))
            current = []

    if current and score is not None:
        modes.append((score, heavy_coords(current)))

    return modes

rows = []

for a in ligands:
    fixed_path = folder / f"AR_{a}_ligandA.pdbqt"
    a_coords = heavy_coords(
        fixed_path.read_text().splitlines()
    )

    for b in ligands:
        if a == b:
            continue

        path = results_dir / f"AR_{a}_to_{b}_docked.pdbqt"
        modes = read_modes(path)

        for mode_number, (score, b_coords) in enumerate(modes, 1):
            distances = np.linalg.norm(
                a_coords[:, None, :] -
                b_coords[None, :, :],
                axis=2
            )

            rows.append({
```

### Source lines 3058–3351
```python
print("\n=== ASSEMBLY SUMMARY ===")
print("Successful:", sum(r[2] == "PASS" for r in records), "/12")
print("Saved:", csv_path)

from pathlib import Path
import numpy as np
import csv

from rdkit import Chem
from openmm import (
    CustomExternalForce,
    LangevinMiddleIntegrator,
    unit
)
from openmm.app import (
    PDBFile,
    ForceField,
    Modeller,
    Simulation,
    NoCutoff
)
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Molecule

print("=== AR TERNARY COMPLEX ENERGY MINIMIZATION ===")

base = Path("/content")
folder = base / "AR_ternary"
output = folder / "minimized"
output.mkdir(parents=True, exist_ok=True)

ligands = ["EE", "LNG", "NGMN", "MPA"]
records = []

for a in ligands:
    for b in ligands:
        if a == b:
            continue

        label = f"{a}_to_{b}"
        print(f"\n--- AR {a} → {b} ---")

        try:
            protein_path = folder / f"AR_{label}_protein.pdb"

            ligand_a_path = (
                base / "AR_minimized" /
                f"AR_{a}_minimized_ligand.sdf"
            )

            ligand_b_path = (
                folder / f"AR_{label}_ligandB.sdf"
            )

            protein = PDBFile(str(protein_path))

            rd_a = Chem.SDMolSupplier(
                str(ligand_a_path), removeHs=False
            )[0]

            rd_b = Chem.SDMolSupplier(
                str(ligand_b_path), removeHs=False
            )[0]

            if rd_a is None or rd_b is None:
                raise ValueError("Invalid ligand structure")

            # Add missing hydrogens if necessary.
            rd_a = Chem.AddHs(rd_a, addCoords=True)
            rd_b = Chem.AddHs(rd_b, addCoords=True)

            # Separate OpenFF molecules.
            off_a = Molecule.from_rdkit(
                rd_a,
                allow_undefined_stereo=True
            )

            off_b = Molecule.from_rdkit(
                rd_b,
                allow_undefined_stereo=True
            )

            # Build topology and coordinates.
            modeller = Modeller(
                protein.topology,
                protein.positions
            )

            for rd_mol in [rd_a, rd_b]:
                conf = rd_mol.GetConformer()
                xyz = conf.GetPositions()

                # Create a temporary PDB representation.
                pdb_block = Chem.MolToPDBBlock(rd_mol)

                import io
                ligand_pdb = PDBFile(io.StringIO(pdb_block))

                positions = [
                    unit.Vec3(*coord) * unit.angstrom
                    for coord in xyz
                ]

                modeller.add(
                    ligand_pdb.topology,
                    positions
                )

            ff = ForceField("amber14-all.xml")

            generator = SMIRNOFFTemplateGenerator(
                molecules=[off_a, off_b],
                forcefield="openff-2.2.1"
            )

            ff.registerTemplateGenerator(
                generator.generator
            )

            system = ff.createSystem(
                modeller.topology,
                nonbondedMethod=NoCutoff,
                constraints=None
            )

            # Restrain protein heavy atoms.
            restraint = CustomExternalForce(
                "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
            )

            restraint.addGlobalParameter(
                "k", 1000.0
            )

            for parameter in ["x0", "y0", "z0"]:
                restraint.addPerParticleParameter(parameter)

            positions_nm = modeller.positions.value_in_unit(
                unit.nanometer
            )

            for i, atom in enumerate(modeller.topology.atoms()):
                if (
                    atom.residue.chain.index == 0
                    and atom.element.symbol != "H"
                ):
                    restraint.addParticle(
                        i, positions_nm[i]
                    )

            system.addForce(restraint)

            integrator = LangevinMiddleIntegrator(
                300 * unit.kelvin,
                1 / unit.picosecond,
                0.002 * unit.picoseconds
            )

            simulation = Simulation(
                modeller.topology,
                system,
                integrator
            )

            simulation.context.setPositions(
                modeller.positions
            )

            initial_energy = (
                simulation.context.getState(
                    getEnergy=True
                ).getPotentialEnergy().value_in_unit(
                    unit.kilojoule_per_mole
                )
            )

            simulation.minimizeEnergy(
                maxIterations=1000
            )

            state = simulation.context.getState(
                getEnergy=True,
                getPositions=True
            )

            final_energy = (
                state.getPotentialEnergy().value_in_unit(
                    unit.kilojoule_per_mole
                )
            )

            out_pdb = output / f"AR_{label}_minimized.pdb"

            with open(out_pdb, "w") as f:
                PDBFile.writeFile(
                    modeller.topology,
                    state.getPositions(),
                    f,
                    keepIds=True
                )

            records.append([
                a, b, "PASS",
                round(initial_energy, 2),
                round(final_energy, 2)
            ])

            print(f"Initial energy: {initial_energy:.2f} kJ/mol")
            print(f"Final energy:   {final_energy:.2f} kJ/mol")
            print("Minimization: PASS")

        except Exception as e:
            print("Minimization: FAILED")
            print(type(e).__name__, str(e)[:400])
            records.append([a, b, "FAILED", "", ""])

csv_path = output / "AR_ternary_minimization_energies.csv"

with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "Ligand_A",
        "Ligand_B",
        "Status",
        "Initial_Energy_kJ_mol",
        "Final_Energy_kJ_mol"
    ])
    writer.writerows(records)

print("\n=== AR TERNARY MINIMIZATION SUMMARY ===")

for row in records:
    print(
        f"{row[0]:4s} → {row[1]:4s} | "
        f"{row[2]:6s} | Final energy: {row[4]}"
    )

print("\nSuccessful:",
      sum(r[2] == "PASS" for r in records), "/12")
print("Saved:", csv_path)

from pathlib import Path
import io
import csv
import numpy as np

from rdkit import Chem
from openmm import (
    Vec3,
    CustomExternalForce,
    LangevinMiddleIntegrator,
    unit
)
from openmm.app import (
    PDBFile,
    ForceField,
    Modeller,
    Simulation,
    NoCutoff
)
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Molecule

print("=== AR TERNARY COMPLEX MINIMIZATION — CORRECTED ===")

base = Path("/content")
folder = base / "AR_ternary"
output = folder / "minimized"
output.mkdir(parents=True, exist_ok=True)

ligands = ["EE", "LNG", "NGMN", "MPA"]
records = []

for a in ligands:
    for b in ligands:
        if a == b:
            continue

        label = f"{a}_to_{b}"
        print(f"\n--- AR {a} → {b} ---")

        try:
            protein = PDBFile(
                str(folder / f"AR_{label}_protein.pdb")
            )

            rd_molecules = []

            for path in [
                base / "AR_minimized" /
                f"AR_{a}_minimized_ligand.sdf",
                folder / f"AR_{label}_ligandB.sdf"
            ]:
                mol = Chem.SDMolSupplier(
```

### Source lines 3367–3560
```python
            ]

            modeller = Modeller(
                protein.topology, protein.positions
            )

            for mol in rd_molecules:
                pdb_block = Chem.MolToPDBBlock(mol)
                ligand_pdb = PDBFile(io.StringIO(pdb_block))

                xyz = mol.GetConformer().GetPositions()

                positions = unit.Quantity(
                    [Vec3(*map(float, coord)) for coord in xyz],
                    unit.angstrom
                )

                modeller.add(
                    ligand_pdb.topology,
                    positions
                )

            ff = ForceField("amber14-all.xml")

            generator = SMIRNOFFTemplateGenerator(
                molecules=off_molecules,
                forcefield="openff-2.2.1"
            )

            ff.registerTemplateGenerator(generator.generator)

            system = ff.createSystem(
                modeller.topology,
                nonbondedMethod=NoCutoff,
                constraints=None
            )

            restraint = CustomExternalForce(
                "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
            )

            restraint.addGlobalParameter("k", 1000.0)

            for name in ["x0", "y0", "z0"]:
                restraint.addPerParticleParameter(name)

            xyz_nm = modeller.positions.value_in_unit(
                unit.nanometer
            )

            protein_atom_count = sum(
                1 for _ in protein.topology.atoms()
            )

            for i, atom in enumerate(modeller.topology.atoms()):
                if (
                    i < protein_atom_count
                    and atom.element.symbol != "H"
                ):
                    restraint.addParticle(i, xyz_nm[i])

            system.addForce(restraint)

            integrator = LangevinMiddleIntegrator(
                300 * unit.kelvin,
                1 / unit.picosecond,
                0.002 * unit.picoseconds
            )

            simulation = Simulation(
                modeller.topology, system, integrator
            )

            simulation.context.setPositions(modeller.positions)

            initial = simulation.context.getState(
                getEnergy=True
            ).getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole
            )

            simulation.minimizeEnergy(maxIterations=1000)

            state = simulation.context.getState(
                getEnergy=True, getPositions=True
            )

            final = state.getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole
            )

            with open(
                output / f"AR_{label}_minimized.pdb", "w"
            ) as f:
                PDBFile.writeFile(
                    modeller.topology,
                    state.getPositions(),
                    f,
                    keepIds=True
                )

            records.append([
                a, b, "PASS",
                round(initial, 2),
                round(final, 2)
            ])

            print(f"Initial energy: {initial:.2f} kJ/mol")
            print(f"Final energy:   {final:.2f} kJ/mol")
            print("Minimization: PASS")

            del simulation, integrator

        except Exception as e:
            print("Minimization: FAILED")
            print(type(e).__name__, str(e)[:500])
            records.append([a, b, "FAILED", "", ""])

csv_path = output / "AR_ternary_minimization_energies.csv"

with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "Ligand_A", "Ligand_B", "Status",
        "Initial_Energy_kJ_mol",
        "Final_Energy_kJ_mol"
    ])
    writer.writerows(records)

print("\n=== AR TERNARY MINIMIZATION SUMMARY ===")

for row in records:
    print(
        f"{row[0]:4s} → {row[1]:4s} | "
        f"{row[2]:6s} | Final energy: {row[4]}"
    )

print(
    "\nSuccessful:",
    sum(r[2] == "PASS" for r in records),
    "/12"
)
print("Saved:", csv_path)

from pathlib import Path
import numpy as np

print("=== AR POST-MINIMIZATION SPATIAL VALIDATION ===")

base = Path("/content/AR_ternary")
folder = base / "minimized"

ligands = ["EE", "LNG", "NGMN", "MPA"]

heavy_counts = {
    "EE": 22,
    "LNG": 23,
    "NGMN": 24,
    "MPA": 28
}

print(
    "\nPair          Min A–B distance   "
    "Contacts <2 Å   Status"
)

for a in ligands:
    for b in ligands:
        if a == b:
            continue

        path = folder / f"AR_{a}_to_{b}_minimized.pdb"

        try:
            lines = path.read_text().splitlines()

            # OpenMM writes ligand residues after protein.
            # Identify nonprotein residues using the
            # expected heavy-atom counts and residue order.
            from openmm.app import PDBFile

            pdb = PDBFile(str(path))
            atoms = list(pdb.topology.atoms())
            positions = pdb.positions.value_in_unit(
                __import__("openmm").unit.angstrom
            )

            residues = list(pdb.topology.residues())

            ligand_residues = [
                res for res in residues
                if res.name not in {
                    "ALA", "ARG", "ASN", "ASP", "CYS",
                    "GLU", "GLN", "GLY", "HIS", "ILE",
```

### Source lines 3599–3687
```python
            status = (
                "PASS" if clashes == 0 else "REVIEW"
            )

            print(
                f"{a:4s} → {b:4s}     "
                f"{minimum:6.2f} Å          "
                f"{clashes:3d}          "
                f"{status}"
            )

        except Exception as e:
            print(
                f"{a:4s} → {b:4s} | ERROR: "
                f"{str(e)[:120]}"
            )

print("\n=== VALIDATION COMPLETE ===")

from pathlib import Path
from rdkit import Chem
import numpy as np
import csv

print("=== AR SEQUENTIAL DOCKING CENTROID ANALYSIS ===")

base = Path("/content")
results = base / "AR_sequential/docking_results"
ligands = ["EE", "LNG", "NGMN", "MPA"]

rows = []

for a in ligands:
    ligand_a_path = (
        base / "AR_minimized" /
        f"AR_{a}_minimized_ligand.sdf"
    )

    mol_a = Chem.SDMolSupplier(
        str(ligand_a_path), removeHs=False
    )[0]

    if mol_a is None:
        raise ValueError(f"Invalid Ligand A: {a}")

    coords_a = mol_a.GetConformer().GetPositions()
    heavy_a = np.array([
        coords_a[atom.GetIdx()]
        for atom in mol_a.GetAtoms()
        if atom.GetAtomicNum() > 1
    ])

    centroid_a = heavy_a.mean(axis=0)

    for b in ligands:
        if a == b:
            continue

        ligand_b_path = (
            base / "AR_ternary" /
            f"AR_{a}_to_{b}_ligandB.sdf"
        )

        mol_b = Chem.SDMolSupplier(
            str(ligand_b_path), removeHs=False
        )[0]

        if mol_b is None:
            raise ValueError(f"Invalid Ligand B: {a} → {b}")

        coords_b = mol_b.GetConformer().GetPositions()
        heavy_b = np.array([
            coords_b[atom.GetIdx()]
            for atom in mol_b.GetAtoms()
            if atom.GetAtomicNum() > 1
        ])

        centroid_b = heavy_b.mean(axis=0)

        centroid_distance = float(
            np.linalg.norm(centroid_a - centroid_b)
        )

        min_distance = float(
            np.linalg.norm(
                heavy_a[:, None, :] -
                heavy_b[None, :, :],
                axis=2
            ).min()
```

### Source lines 4582–4758
```python

    prepared = read_pdbqt_atoms(
        base / "AR_sequential" /
        f"AR_{ligand}_ligandA.pdbqt"
    )

    original_elements = [a[1] for a in original]
    prepared_elements = [a[1] for a in prepared]

    same_order = original_elements == prepared_elements

    print(f"\n--- {ligand} ---")
    print("Original heavy atoms:", len(original))
    print("Prepared heavy atoms:", len(prepared))
    print("Element sequence matches:", same_order)

    if len(original) == len(prepared):
        original_xyz = np.array([a[2] for a in original])
        prepared_xyz = np.array([a[2] for a in prepared])

        distances = np.linalg.norm(
            original_xyz - prepared_xyz, axis=1
        )

        print("Index-matched coordinate RMSD:",
              round(float(np.sqrt(np.mean(distances**2))), 3),
              "Å")

        print("Maximum index-matched displacement:",
              round(float(distances.max()), 3), "Å")

        print("First five original atom names:",
              [a[0] for a in original[:5]])

        print("First five prepared atom names:",
              [a[0] for a in prepared[:5]])

print("\n=== CHECK COMPLETE ===")

from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment

print("=== AR LIGAND GEOMETRIC CORRESPONDENCE CHECK ===")

base = Path("/content")

def get_atoms(path, pdbqt=False):
    atoms = []

    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue

        if pdbqt:
            atom_type = line[77:79].strip().upper()
            if atom_type in ("H", "HD"):
                continue

            element = {
                "A": "C",
                "C": "C",
                "OA": "O",
                "NA": "N",
                "SA": "S"
            }.get(atom_type, atom_type)
        else:
            element = line[76:78].strip().upper()
            if element == "H":
                continue

        xyz = np.array([
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        ])

        atoms.append((element, xyz))

    return atoms

for ligand in ["EE", "LNG", "NGMN", "MPA"]:
    original = get_atoms(
        base / "AR_rigid_control" /
        f"AR_{ligand}_original_ligand.pdb"
    )

    prepared = get_atoms(
        base / "AR_sequential" /
        f"AR_{ligand}_ligandA.pdbqt",
        pdbqt=True
    )

    n = len(original)

    if n != len(prepared):
        print(f"{ligand}: ATOM COUNT MISMATCH")
        continue

    cost = np.full((n, n), 1e6)

    for i, (elem_i, xyz_i) in enumerate(original):
        for j, (elem_j, xyz_j) in enumerate(prepared):
            if elem_i == elem_j:
                cost[i, j] = np.linalg.norm(xyz_i - xyz_j)

    rows, cols = linear_sum_assignment(cost)

    distances = cost[rows, cols]

    print(f"\n--- {ligand} ---")
    print("Matched heavy atoms:", len(rows))
    print("RMSD after matching:",
          round(float(np.sqrt(np.mean(distances**2))), 3),
          "Å")
    print("Maximum matched displacement:",
          round(float(distances.max()), 3), "Å")
    print("Atoms displaced >2 Å:",
          int(np.sum(distances > 2)))

    if np.any(distances > 2):
        print("Status: REQUIRES FURTHER VALIDATION")
    else:
        print("Status: GEOMETRIC MATCH PLAUSIBLE")

print("\n=== CHECK COMPLETE ===")

from pathlib import Path
from collections import Counter
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

print("=== AR LIGAND ATOM MAPPING VALIDATION ===")

base = Path("/content")
ligands = ["EE", "LNG", "NGMN", "MPA"]

def read_atoms(path, pdbqt=False):
    atoms = []

    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue

        if pdbqt:
            atom_type = line[77:79].strip().upper()
            if atom_type in ("H", "HD"):
                continue
            element = {
                "A": "C", "C": "C",
                "OA": "O", "NA": "N",
                "SA": "S"
            }.get(atom_type, atom_type)
        else:
            element = line[76:78].strip().upper()
            if element == "H":
                continue

        xyz = np.array([
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        ])

        atoms.append((element, xyz))

    return atoms

# Approximate maximum covalent bond distances (Å).
# This is a geometric connectivity sanity check,
# not a replacement for explicit chemical bond topology.
bond_limits = {
    ("C", "C"): 1.85,
    ("C", "N"): 1.75,
    ("C", "O"): 1.70,
    ("N", "N"): 1.70,
```

### Source lines 6200–6431
```python
        f"maximum {nearest.max():.3f} Å | "
        f">1.6 Å: {np.sum(nearest > 1.6)}"
    )

    if np.any(nearest > 1.6):
        print(
            f"  WARNING: {ligand} still has hydrogen-distance "
            "outliers; file not written."
        )
        continue

    output_path.write_text(
        "\n".join(updated) + "\n"
    )

    print(f"  Saved: {output_path}")

print("\n=== HYDROGEN REPOSITIONING COMPLETE ===")

from pathlib import Path
import subprocess
import csv
import re

print("=== AR ORIGINAL RIGID SEQUENTIAL DOCKING ===")

base = Path("/content")
rigid = base / "AR_rigid_control"
sequential = base / "AR_sequential"
results_dir = rigid / "docking_results"
results_dir.mkdir(parents=True, exist_ok=True)

vina = base / "vina_1.2.7"

ligands = ["EE", "LNG", "NGMN", "MPA"]

records = []

for a in ligands:
    receptor = rigid / f"AR_{a}_final_rigid_receptor.pdbqt"

    if not receptor.exists():
        raise FileNotFoundError(receptor)

    for b in ligands:
        if a == b:
            continue

        ligand_b = sequential / f"AR_{b}_ligandB.pdbqt"

        output = results_dir / f"AR_{a}_to_{b}_rigid_docked.pdbqt"
        log = results_dir / f"AR_{a}_to_{b}_rigid_vina.log"

        if not ligand_b.exists():
            raise FileNotFoundError(ligand_b)

        command = [
            str(vina),
            "--receptor", str(receptor),
            "--ligand", str(ligand_b),
            "--center_x", "27.092",
            "--center_y", "2.279",
            "--center_z", "4.950",
            "--size_x", "30",
            "--size_y", "30",
            "--size_z", "30",
            "--exhaustiveness", "16",
            "--num_modes", "10",
            "--energy_range", "5",
            "--seed", "2026",
            "--cpu", "2",
            "--out", str(output),
            "--log", str(log)
        ]

        print(f"Docking {a} → {b} ...", flush=True)

        run = subprocess.run(
            command,
            capture_output=True,
            text=True
        )

        if run.returncode != 0:
            print((run.stdout + run.stderr)[-3000:])
            raise RuntimeError(f"Docking failed: {a} → {b}")

        log_text = log.read_text()

        match = re.search(
            r"^\s*1\s+(-?\d+\.\d+)\s+",
            log_text,
            flags=re.MULTILINE
        )

        if not match:
            raise RuntimeError(
                f"Could not extract docking score: {a} → {b}"
            )

        score = float(match.group(1))

        records.append({
            "Ligand_A": a,
            "Ligand_B": b,
            "Rigid_Vina_Score_kcal_mol": score,
            "Output_File": str(output)
        })

        print(f"  Best score: {score:.3f} kcal/mol")

csv_path = results_dir / "AR_rigid_sequential_results.csv"

with csv_path.open("w", newline="") as file:
    writer = csv.DictWriter(
        file,
        fieldnames=[
            "Ligand_A",
            "Ligand_B",
            "Rigid_Vina_Score_kcal_mol",
            "Output_File"
        ]
    )
    writer.writeheader()
    writer.writerows(records)

print("\n=== ALL 12 RIGID DOCKINGS COMPLETE ===")
print("Results saved:", csv_path)

from pathlib import Path
import subprocess
import csv
import re

print("=== AR ORIGINAL RIGID SEQUENTIAL DOCKING ===")

base = Path("/content")
rigid = base / "AR_rigid_control"
sequential = base / "AR_sequential"
results_dir = rigid / "docking_results"
results_dir.mkdir(parents=True, exist_ok=True)

vina = base / "vina_1.2.7"
ligands = ["EE", "LNG", "NGMN", "MPA"]

records = []

for a in ligands:
    receptor = rigid / f"AR_{a}_final_rigid_receptor.pdbqt"

    if not receptor.exists():
        raise FileNotFoundError(receptor)

    for b in ligands:
        if a == b:
            continue

        ligand_b = sequential / f"AR_{b}_ligandB.pdbqt"

        output = results_dir / f"AR_{a}_to_{b}_rigid_docked.pdbqt"
        log = results_dir / f"AR_{a}_to_{b}_rigid_vina.log"

        if not ligand_b.exists():
            raise FileNotFoundError(ligand_b)

        command = [
            str(vina),
            "--receptor", str(receptor),
            "--ligand", str(ligand_b),
            "--center_x", "27.092",
            "--center_y", "2.279",
            "--center_z", "4.950",
            "--size_x", "30",
            "--size_y", "30",
            "--size_z", "30",
            "--exhaustiveness", "16",
            "--num_modes", "10",
            "--energy_range", "5",
            "--seed", "2026",
            "--cpu", "2",
            "--out", str(output)
        ]

        print(f"Docking {a} → {b} ...", flush=True)

        run = subprocess.run(
            command,
            capture_output=True,
            text=True
        )

        log_text = run.stdout + "\n" + run.stderr
        log.write_text(log_text)

        if run.returncode != 0:
            print(log_text[-3000:])
            raise RuntimeError(f"Docking failed: {a} → {b}")

        # Extract the best mode score from Vina's result table.
        match = re.search(
            r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+",
            log_text,
            flags=re.MULTILINE
        )

        if not match:
            raise RuntimeError(
                f"Could not extract docking score: {a} → {b}. "
                f"Inspect {log}"
            )

        score = float(match.group(1))

        records.append({
            "Ligand_A": a,
            "Ligand_B": b,
            "Rigid_Vina_Score_kcal_mol": score,
            "Output_File": str(output)
        })

        print(f"  Best score: {score:.3f} kcal/mol", flush=True)

csv_path = results_dir / "AR_rigid_sequential_results.csv"

with csv_path.open("w", newline="") as file:
    writer = csv.DictWriter(
        file,
        fieldnames=[
            "Ligand_A",
            "Ligand_B",
            "Rigid_Vina_Score_kcal_mol",
            "Output_File"
```


---

# Cytochrome P450 3A4 (CYP3A4)

# CYP3A4 — original Dock–EM–Dock code excerpts
**Source:** `cypenergyminimization.py` (27,841 lines).
**Provenance:** All code excerpts below are copied verbatim from the uploaded Colab-exported Python source, with original line ranges. Sections are excerpts, not a standalone runnable script. Original complete file is included in the ZIP.
**Scope:** The original notebook may contain debugging, corrections, alternate parameterizations, or non-executable Colab shell commands. These excerpts do not imply all code paths were executed successfully.

## Original code excerpts
### Source lines 2548–2707
```python
        initial_nm = molecule.conformers[0].m_as(
            offunit.nanometer
        )

        integrator = openmm.VerletIntegrator(
            0.001 * unit.picoseconds
        )

        context = openmm.Context(
            system,
            integrator,
            openmm.Platform.getPlatformByName("Reference")
        )

        context.setPositions(
            initial_nm * unit.nanometer
        )

        initial_energy = (
            context.getState(getEnergy=True)
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
        )

        openmm.LocalEnergyMinimizer.minimize(
            context,
            tolerance=10.0,
            maxIterations=1000
        )

        final_state = context.getState(
            getEnergy=True,
            getPositions=True
        )

        final_energy = (
            final_state.getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
        )

        final_nm = final_state.getPositions(
            asNumpy=True
        ).value_in_unit(unit.nanometer)

        displacement_angstrom = (
            final_nm - initial_nm
        ) * 10.0

        rmsd = float(np.sqrt(
            np.mean(
                np.sum(displacement_angstrom**2, axis=1)
            )
        ))

        if not np.isfinite(final_energy):
            raise RuntimeError("Invalid final energy")

        if final_energy > initial_energy + 0.01:
            raise RuntimeError("Energy increased unexpectedly")

        # Save a separate minimized SDF.
        molecule._conformers = []
        molecule.add_conformer(
            final_nm * offunit.nanometer
        )

        output_path = (
            f"/content/CYP3A4_{code}_ligand_only_minimized.sdf"
        )

        molecule.to_file(
            output_path,
            file_format="SDF"
        )

        print(
            "Initial energy:",
            f"{initial_energy:.4f} kJ/mol",
            flush=True
        )

        print(
            "Final energy:",
            f"{final_energy:.4f} kJ/mol",
            flush=True
        )

        print(
            "Energy change:",
            f"{final_energy - initial_energy:.4f} kJ/mol",
            flush=True
        )

        print(
            "Coordinate RMS displacement:",
            f"{rmsd:.4f} Å",
            flush=True
        )

        print("Saved:", output_path, flush=True)
        print("STATUS: MINIMIZED", flush=True)

        del context
        del integrator

        success_count += 1

    except Exception:
        print("STATUS: FAILED", flush=True)
        traceback.print_exc()

print(
    f"\nSUCCESSFULLY MINIMIZED: {success_count}/4",
    flush=True
)

if success_count == 4:
    print("=== STEP 31 COMPLETE ===", flush=True)
else:
    print("=== STEP 31 INCOMPLETE ===", flush=True)
    sys.exit(1)
"""

script_path = Path("/content/cyp3a4_step31.py")
script_path.write_text(script)

env = os.environ.copy()
env["MPLBACKEND"] = "Agg"
env["USE_TF"] = "0"
env["USE_FLAX"] = "0"

result = subprocess.run(
    ["/usr/local/bin/python", str(script_path)],
    capture_output=True,
    text=True,
    env=env,
    timeout=600
)

print(result.stdout)

if result.stderr:
    print("\n=== ERROR DETAILS ===")
    print(result.stderr[-7000:])

print("\nProcess return code:", result.returncode)

from pathlib import Path
from collections import Counter
import numpy as np

print("=== CYP3A4 STEP 32: RECEPTOR–HEME STRUCTURAL VALIDATION ===")

files = {
    "EE": "/content/1TQN_ee_out_1.-10.6.complex.pdb",
    "LNG": "/content/1TQN_lng_out_2.-8.5.complex.pdb",
    "NGMN": "/content/1TQN_ngmn_out_1.-9.5.complex.pdb",
    "MPA": "/content/1TQN_mpa_out_2.-8.0.complex.pdb"
}

```

### Source lines 2903–3009
```python
            float(line[38:46]),
            float(line[46:54])
        ])

    pdb_coords = np.array(pdb_coords)

    sdf_elements = [
        atom.GetSymbol().upper()
        for atom in sdf_atoms
    ]

    print("SDF heavy atoms:", len(sdf_atoms))
    print("PDB heavy atoms:", len(pdb_elements))

    if len(sdf_atoms) != len(pdb_elements):
        print("STATUS: ATOM COUNT MISMATCH")
        continue

    if sdf_elements != pdb_elements:
        print("STATUS: ATOM ORDER/ELEMENT MISMATCH")
        continue

    differences = sdf_coords - pdb_coords

    rmsd = np.sqrt(
        np.mean(np.sum(differences**2, axis=1))
    )

    maximum = np.max(
        np.linalg.norm(differences, axis=1)
    )

    print(
        "Direct coordinate RMSD:",
        f"{rmsd:.6f} Å"
    )

    print(
        "Maximum atom displacement:",
        f"{maximum:.6f} Å"
    )

    if rmsd < 0.01:
        print("STATUS: ALIGNED")
    else:
        print("STATUS: ALIGNMENT REQUIRES REVIEW")

print("\n=== STEP 33 COMPLETE ===")

from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from rdkit import Chem

print("=== CYP3A4 STEP 34: INITIAL RECEPTOR–LIGAND CONTACT ANALYSIS ===")

complex_files = {
    "EE": "/content/1TQN_ee_out_1.-10.6.complex.pdb",
    "LNG": "/content/1TQN_lng_out_2.-8.5.complex.pdb",
    "NGMN": "/content/1TQN_ngmn_out_1.-9.5.complex.pdb",
    "MPA": "/content/1TQN_mpa_out_2.-8.0.complex.pdb"
}

def read_receptor_atoms(filename):
    protein = []
    heme = []
    iron = []

    for line in Path(filename).read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue

        element = line[76:78].strip().upper()

        if element in ("H", "D"):
            continue

        atom = {
            "name": line[12:16].strip(),
            "resname": line[17:20].strip(),
            "resid": line[22:26].strip(),
            "element": element,
            "xyz": np.array([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ])
        }

        if line.startswith("ATOM  "):
            protein.append(atom)

        elif atom["resname"] == "HEM":
            if element == "FE":
                iron.append(atom)
            else:
                heme.append(atom)

    return protein, heme, iron

for code, filename in complex_files.items():

    print(f"\n=== {code} ===")

    mol = Chem.SDMolSupplier(
        f"/content/CYP3A4_{code}_prepared.sdf",
        removeHs=False
```

### Source lines 6630–6888
```python
result = subprocess.run(
    [
        "/content/cyp3a4_openff_env/bin/python",
        "-u",
        str(script_path)
    ],
    capture_output=True,
    text=True,
    env=env,
    timeout=180
)

print(result.stdout)

if result.stderr:
    print("\n=== DIAGNOSTIC MESSAGES ===")
    print(result.stderr[-4000:])

print("Process return code:", result.returncode)

import os
import subprocess
from pathlib import Path

print("=== CYP3A4 STEP 57: LIGAND-ONLY ENERGY MINIMIZATION ===")

script = r'''
from pathlib import Path
import numpy as np
import openmm
from openmm import unit
from openff.toolkit import Molecule, ForceField
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from openff.nagl_models import list_available_nagl_models
from openff.interchange import Interchange
from openff.units import unit as offunit

print("=== FORCE FIELD SETUP ===", flush=True)

ff = ForceField("openff_unconstrained-2.2.1.offxml")

models = list(list_available_nagl_models())

model_path = next(
    path for path in models
    if path.name == "openff-gnn-am1bcc-1.0.0.pt"
)

wrapper = NAGLToolkitWrapper()

print("Force field: OpenFF Sage 2.2.1", flush=True)
print("Charge model:", model_path.name, flush=True)
print("OpenMM:", openmm.__version__, flush=True)

ligands = {
    "EE": "/content/CYP3A4_EE_prepared.sdf",
    "LNG": "/content/CYP3A4_LNG_prepared.sdf",
    "NGMN": "/content/CYP3A4_NGMN_E_geometry_prepared.sdf",
    "MPA": "/content/CYP3A4_MPA_prepared.sdf"
}

results = []

for code, filename in ligands.items():

    print(f"\n=== {code} MINIMIZATION ===", flush=True)

    mol = Molecule.from_file(
        filename,
        allow_undefined_stereo=False
    )

    mol.assign_partial_charges(
        partial_charge_method=str(model_path),
        toolkit_registry=wrapper
    )

    interchange = Interchange.from_smirnoff(
        force_field=ff,
        topology=mol.to_topology(),
        charge_from_molecules=[mol]
    )

    system = interchange.to_openmm()

    initial_nm = np.asarray(
        mol.conformers[0].m_as(offunit.nanometer)
    )

    integrator = openmm.VerletIntegrator(
        0.001 * unit.picoseconds
    )

    platform = openmm.Platform.getPlatformByName(
        "Reference"
    )

    context = openmm.Context(
        system,
        integrator,
        platform
    )

    context.setPositions(
        initial_nm * unit.nanometer
    )

    initial_energy = (
        context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )

    openmm.LocalEnergyMinimizer.minimize(
        context,
        tolerance=10.0,
        maxIterations=1000
    )

    state = context.getState(
        getEnergy=True,
        getPositions=True
    )

    final_energy = (
        state.getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )

    final_nm = state.getPositions(
        asNumpy=True
    ).value_in_unit(unit.nanometer)

    displacement_angstrom = (
        final_nm - initial_nm
    ) * 10.0

    rms_displacement = float(np.sqrt(
        np.mean(
            np.sum(displacement_angstrom**2, axis=1)
        )
    ))

    energy_change = final_energy - initial_energy

    mol._conformers = [
        final_nm * offunit.nanometer
    ]

    output = Path(
        f"/content/CYP3A4_{code}_ligand_only_minimized.sdf"
    )

    mol.to_file(
        str(output),
        file_format="SDF"
    )

    print(
        f"Initial energy: {initial_energy:.4f} kJ/mol",
        flush=True
    )

    print(
        f"Final energy: {final_energy:.4f} kJ/mol",
        flush=True
    )

    print(
        f"Energy change: {energy_change:.4f} kJ/mol",
        flush=True
    )

    print(
        f"All-atom RMS displacement: "
        f"{rms_displacement:.4f} Å",
        flush=True
    )

    print("Saved:", output, flush=True)

    results.append(
        (
            code,
            initial_energy,
            final_energy,
            energy_change,
            rms_displacement
        )
    )

    del context
    del integrator

print("\n=== MINIMIZATION SUMMARY ===", flush=True)

print(
    f"{'Ligand':<8}"
    f"{'Initial':>13}"
    f"{'Final':>13}"
    f"{'Change':>13}"
    f"{'RMS (Å)':>13}",
    flush=True
)

for code, initial, final, change, rms in results:
    print(
        f"{code:<8}"
        f"{initial:>13.4f}"
        f"{final:>13.4f}"
        f"{change:>13.4f}"
        f"{rms:>13.4f}",
        flush=True
    )

print("\nSTATUS: MINIMIZATION COMPLETE", flush=True)
print("=== STEP 57 COMPLETE ===", flush=True)
'''

script_path = Path("/content/cyp3a4_step57.py")
script_path.write_text(script)

env = os.environ.copy()
env["MPLBACKEND"] = "Agg"
env["OMP_NUM_THREADS"] = "2"
env["OPENBLAS_NUM_THREADS"] = "2"

try:
    result = subprocess.run(
        [
            "/content/cyp3a4_openff_env/bin/python",
            "-u",
            str(script_path)
        ],
        capture_output=True,
        text=True,
        env=env,
        timeout=900
    )

    print(result.stdout)

    if result.stderr:
        print("\n=== DIAGNOSTIC MESSAGES ===")
        print(result.stderr[-6000:])

    print(
        "\nProcess return code:",
        result.returncode
    )

except subprocess.TimeoutExpired as exc:
    if exc.stdout:
        print(exc.stdout.decode(errors="replace")
              if isinstance(exc.stdout, bytes)
              else exc.stdout)

    print("STATUS: MINIMIZATION TIMED OUT")

```

### Source lines 8454–8542
```python
        original_nm,
        receptor_nm
    ])

    integrator = openmm.VerletIntegrator(
        0.001 * unit.picoseconds
    )

    context = openmm.Context(
        system,
        integrator,
        openmm.Platform.getPlatformByName("Reference")
    )

    context.setPositions(
        positions_nm * unit.nanometer
    )

    initial_energy = (
        context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )

    openmm.LocalEnergyMinimizer.minimize(
        context,
        tolerance=10.0,
        maxIterations=1000
    )

    state = context.getState(
        getEnergy=True,
        getPositions=True
    )

    final_energy = (
        state.getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )

    final_positions_nm = state.getPositions(
        asNumpy=True
    ).value_in_unit(unit.nanometer)

    final_ligand_nm = final_positions_nm[
        :mol.n_atoms
    ]

    displacement = (
        final_ligand_nm - original_nm
    ) * 10.0

    heavy_rms = np.sqrt(np.mean(
        np.sum(
            displacement[heavy_indices]**2,
            axis=1
        )
    ))

    mol._conformers = [
        final_ligand_nm * offunit.nanometer
    ]

    output = (
        f"/content/CYP3A4_{code}_"
        f"fixed_receptor_refined.sdf"
    )

    mol.to_file(
        output,
        file_format="SDF"
    )

    print(
        "Receptor heavy atoms:",
        len(receptor),
        flush=True
    )

    print(
        "Ligand-receptor steric pairs:",
        pair_count,
        flush=True
    )

    print(
        f"Initial model energy: "
        f"{initial_energy:.4f} kJ/mol",
        flush=True
```

### Source lines 11217–11305
```python
print("Vina version:", metadata.version("vina"))

for method in [
    "__init__",
    "set_receptor",
    "set_ligand_from_file",
    "compute_vina_maps",
    "score",
    "optimize",
    "dock"
]:
    obj = getattr(Vina, method, None)
    if obj is None:
        print(f"{method}: NOT AVAILABLE")
    else:
        try:
            signature = inspect.signature(obj)
        except (TypeError, ValueError):
            signature = "signature unavailable"
        print(f"{method}: {signature}")

print("\nVina scoring-function initialization:")
for scoring in ["vina", "vinardo", "ad4"]:
    try:
        engine = Vina(sf_name=scoring, cpu=1, seed=42)
        print(f"{scoring}: INITIALIZATION SUCCESSFUL")
    except Exception as error:
        print(f"{scoring}: {type(error).__name__}: {error}")
"""

env = os.environ.copy()
env["MPLBACKEND"] = "Agg"

result = subprocess.run(
    [python_exe, "-c", inspection],
    capture_output=True,
    text=True,
    env=env,
    timeout=120
)

print(result.stdout)

if result.stderr:
    print("Diagnostics:")
    print(result.stderr[-3000:])

print("Inspection return code:", result.returncode)

print("\n=== REFINED LIGAND DISTANCES TO HEME ===")

receptor = Path(
    "/content/CYP3A4_protein_heme_unprepared.pdb"
)

if not receptor.exists():
    raise FileNotFoundError(receptor)

heme_atoms = []

for line in receptor.read_text().splitlines():
    if not line.startswith("HETATM"):
        continue

    if line[17:20].strip() != "HEM":
        continue

    heme_atoms.append({
        "name": line[12:16].strip(),
        "element": line[76:78].strip().upper(),
        "xyz": np.array([
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        ])
    })

heme_noniron = [
    a for a in heme_atoms
    if a["element"] != "FE"
]

iron = [
    a for a in heme_atoms
    if a["element"] == "FE"
]

if len(heme_atoms) != 43 or len(iron) != 1:
    raise RuntimeError("Unexpected heme composition.")
```

### Source lines 11953–12127
```python
        int(ligand_inside.sum()),
        "/",
        len(ligand_coords)
    )
    print(
        "Heme atoms inside box:",
        int(heme_inside.sum()),
        "/ 43"
    )
    print(
        "Box volume:",
        f"{np.prod(size):.1f} Å³"
    )

    if not np.all(ligand_inside):
        print("WARNING: Some ligand atoms fall outside the box.")

print("\n=== VINA RECEPTOR LOADING TEST ===")

python_exe = "/content/cyp3a4_openff_env/bin/python"

test_code = f"""
from vina import Vina

engine = Vina(
    sf_name="vina",
    cpu=1,
    seed=42,
    verbosity=0
)

engine.set_receptor(
    rigid_pdbqt_filename={str(receptor_path)!r}
)

print("Protein-only receptor loaded successfully.")
"""

env = os.environ.copy()
env["MPLBACKEND"] = "Agg"

result = subprocess.run(
    [python_exe, "-c", test_code],
    capture_output=True,
    text=True,
    env=env,
    timeout=120
)

print(result.stdout)

if result.stderr:
    print("Diagnostics:")
    print(result.stderr[-4000:])

print("Vina receptor loading return code:", result.returncode)

print("\n=== STEP 81 COMPLETE ===")
print("No docking or scoring was performed.")
print("The receptor remains protein-only.")

from pathlib import Path
import subprocess
import os
import json

print("=== CYP3A4 STEP 82: PROTEIN-ONLY LOCAL VINA OPTIMIZATION ===")

python_exe = Path("/content/cyp3a4_openff_env/bin/python")
config_path = Path("/content/CYP3A4_protein_only_docking_boxes.json")
receptor_path = Path("/content/CYP3A4_protein_diagnostic.pdbqt")

for path in [python_exe, config_path, receptor_path]:
    if not path.exists():
        raise FileNotFoundError(path)

configs = json.loads(config_path.read_text())

# Restrict the analysis to ligands whose search boxes
# contain no original heme atoms.
selected = ["LNG", "MPA"]

for code in selected:
    if configs[code]["heme_atoms_inside_box"] != 0:
        raise RuntimeError(
            f"{code}: heme atoms intersect the search box."
        )

script = r"""
import json
import numpy as np
from pathlib import Path
from vina import Vina

configs = json.loads(
    Path(
        "/content/CYP3A4_protein_only_docking_boxes.json"
    ).read_text()
)

results = {}

for code in ["LNG", "MPA"]:

    config = configs[code]

    print(f"\n{'=' * 55}")
    print(f"LIGAND: {code}")
    print("=" * 55)

    engine = Vina(
        sf_name="vina",
        cpu=1,
        seed=42,
        verbosity=0
    )

    engine.set_receptor(
        rigid_pdbqt_filename=config["receptor"]
    )

    engine.set_ligand_from_file(
        config["ligand_pdbqt"]
    )

    engine.compute_vina_maps(
        center=config["center"],
        box_size=config["box_size"],
        spacing=0.375
    )

    initial = np.asarray(
        engine.score(),
        dtype=float
    )

    optimized = np.asarray(
        engine.optimize(max_steps=1000),
        dtype=float
    )

    output = (
        f"/content/CYP3A4_{code}_"
        "protein_only_local_optimized.pdbqt"
    )

    engine.write_pose(
        output,
        overwrite=True
    )

    print(
        "Initial Vina score:",
        f"{initial[0]:.4f} kcal/mol"
    )

    print(
        "Optimized Vina score:",
        f"{optimized[0]:.4f} kcal/mol"
    )

    print(
        "Score change:",
        f"{optimized[0] - initial[0]:+.4f} kcal/mol"
    )

    print(
        "Optimized pose:",
        output
    )

    results[code] = {
        "initial_score_kcal_mol": float(initial[0]),
        "optimized_score_kcal_mol": float(optimized[0]),
        "score_change_kcal_mol": float(
```

### Source lines 12359–12447
```python

    print("\nMinimum distances (Å):")

    distances = {}

    for label, environment in [
        ("Protein", protein),
        ("Heme excluding Fe", heme),
        ("Iron", iron)
    ]:
        before = minimum_distance(original_xyz, environment)
        after = minimum_distance(optimized_xyz, environment)

        distances[label] = {
            "before": before,
            "after": after
        }

        print(
            f"{label}: {before:.3f} → {after:.3f}"
        )

    print("\nVina score components:")

    engine = Vina(
        sf_name="vina",
        cpu=1,
        seed=42,
        verbosity=0
    )

    engine.set_receptor(
        rigid_pdbqt_filename=str(
            root / "CYP3A4_protein_diagnostic.pdbqt"
        )
    )

    engine.compute_vina_maps(
        center=configs[code]["center"],
        box_size=configs[code]["box_size"]
    )

    scores = {}

    for label, ligand_path in [
        ("Starting", initial_path),
        ("Optimized", optimized_path)
    ]:
        engine.set_ligand_from_file(str(ligand_path))
        components = np.asarray(engine.score(), dtype=float)

        scores[label] = components.tolist()

        print(
            f"{label}:",
            np.round(components, 4).tolist()
        )

    results[code] = {
        "heavy_atom_count": len(original),
        "coordinate_rms_displacement_A": rms,
        "maximum_atom_displacement_A": max_displacement,
        "minimum_distances_A": distances,
        "vina_score_components": scores,
        "receptor_model": "protein_only"
    }

output = root / "CYP3A4_step83_pose_score_audit.json"
output.write_text(json.dumps(results, indent=2))

print("\n=== STEP 83 COMPLETE ===")
print("Results saved:", output)
print("No new optimization or docking was performed.")
"""

env = os.environ.copy()
env["MPLBACKEND"] = "Agg"

result = subprocess.run(
    [str(python_exe), "-c", script],
    capture_output=True,
    text=True,
    env=env,
    timeout=300
)

print(result.stdout)

if result.stderr:
```

### Source lines 12470–12558
```python
import numpy as np
import json

root = Path("/content")

configs = json.loads(
    (root / "CYP3A4_protein_only_docking_boxes.json").read_text()
)

results = {}

for code in ["LNG", "MPA"]:

    print(f"\n{'=' * 55}")
    print(f"LIGAND: {code}")
    print("=" * 55)

    config = configs[code]

    files = {
        "Starting": root / f"CYP3A4_{code}_refined_meeko.pdbqt",
        "Optimized": root / f"CYP3A4_{code}_protein_only_local_optimized.pdbqt"
    }

    engine = Vina(
        sf_name="vina",
        cpu=1,
        seed=42,
        verbosity=0
    )

    engine.set_receptor(
        rigid_pdbqt_filename=config["receptor"]
    )

    engine.compute_vina_maps(
        center=config["center"],
        box_size=config["box_size"],
        spacing=0.375
    )

    results[code] = {}

    for label, path in files.items():

        if not path.exists():
            raise FileNotFoundError(path)

        engine.set_ligand_from_file(str(path))

        repeated_scores = []

        for repeat in range(5):
            score = np.asarray(
                engine.score(),
                dtype=float
            )
            repeated_scores.append(score)

        repeated_scores = np.array(repeated_scores)

        total_scores = repeated_scores[:, 0]

        print(f"\n{label}")
        print("File:", path.name)
        print(
            "Repeated total scores:",
            np.round(total_scores, 6).tolist()
        )
        print(
            "Score range:",
            f"{np.ptp(total_scores):.8f} kcal/mol"
        )
        print(
            "Full score components:",
            np.round(repeated_scores[0], 6).tolist()
        )

        atom_lines = [
            line for line in path.read_text().splitlines()
            if line.startswith(("ATOM  ", "HETATM"))
        ]

        print("PDBQT atom records:", len(atom_lines))

        if atom_lines:
            print(
                "Example coordinate fields:",
                repr(atom_lines[0][30:54])
```

### Source lines 12619–12707
```python
    (root / "CYP3A4_protein_only_docking_boxes.json").read_text()
)

step_limits = [50, 200, 1000, 2000]

results = {}

for code in ["LNG", "MPA"]:

    print(f"\n{'=' * 60}")
    print(f"LIGAND: {code}")
    print("=" * 60)

    config = configs[code]

    if config["heme_atoms_inside_box"] != 0:
        raise RuntimeError(
            f"{code}: heme overlaps search box."
        )

    results[code] = {}

    for max_steps in step_limits:

        engine = Vina(
            sf_name="vina",
            cpu=1,
            seed=42,
            verbosity=0
        )

        engine.set_receptor(
            rigid_pdbqt_filename=config["receptor"]
        )

        engine.set_ligand_from_file(
            config["ligand_pdbqt"]
        )

        engine.compute_vina_maps(
            center=config["center"],
            box_size=config["box_size"],
            spacing=0.375
        )

        starting_score = float(engine.score()[0])

        optimized_energy = np.asarray(
            engine.optimize(max_steps=max_steps),
            dtype=float
        )

        pose_path = root / (
            f"CYP3A4_{code}_"
            f"convergence_{max_steps}_steps.pdbqt"
        )

        engine.write_pose(
            str(pose_path),
            overwrite=True
        )

        # Independently score the saved pose.
        engine.set_ligand_from_file(str(pose_path))

        saved_pose_score = float(engine.score()[0])

        results[code][str(max_steps)] = {
            "starting_score_kcal_mol": starting_score,
            "optimizer_score_kcal_mol": float(
                optimized_energy[0]
            ),
            "saved_pose_score_kcal_mol": saved_pose_score,
            "pose_file": str(pose_path)
        }

        print(
            f"max_steps={max_steps:4d} | "
            f"starting={starting_score:8.3f} | "
            f"optimizer={optimized_energy[0]:8.3f} | "
            f"saved_pose={saved_pose_score:8.3f}"
        )

    saved_scores = [
        results[code][str(n)]["saved_pose_score_kcal_mol"]
        for n in step_limits
    ]

    print(
```

### Source lines 14656–14783
```python
        "Refinement_final_model_energy_kJ_mol":
            data["final_energy"],
        "Refinement_model_energy_change_kJ_mol":
            round(
                data["final_energy"] - data["initial_energy"],
                4
            ),
        "Minimum_protein_distance_A":
            round(float(d_protein.min()), 3),
        "Minimum_heme_distance_A":
            round(float(d_heme.min()), 3),
        "Minimum_iron_distance_A":
            round(float(d_iron.min()), 3),
        "Ligand_heme_pairs_below_4A":
            int(np.sum(d_heme < 4)),
        "Ligand_heme_pairs_below_6A":
            int(np.sum(d_heme < 6)),
        "Ligand_heavy_atom_count": len(xyz)
    })

df = pd.DataFrame(rows)

output = root / "CYP3A4_FINAL_CONSOLIDATED_RESULTS.csv"

df.to_csv(output, index=False)

print("\n=== CONSOLIDATED CYP3A4 RESULTS ===")
print(df.to_string(index=False))

print("\n=== IMPORTANT METHOD LABELS ===")
print("Docking scores: original heme-inclusive CB-Dock3")
print("Refinement: fixed-receptor ligand refinement")
print("Distances: refined ligand coordinates to original receptor")
print("Model energies: NOT binding free energies")
print("Full protein-heme minimization: NOT performed")

print("\nSaved:", output)
print("\n=== STEP 93 COMPLETE ===")
"""

result = subprocess.run(
    [str(py), "-c", script],
    capture_output=True,
    text=True,
    env=os.environ.copy(),
    timeout=180
)

print(result.stdout)

if result.stderr:
    print("Diagnostics:")
    print(result.stderr[-4000:])

print("Process return code:", result.returncode)

if result.returncode != 0:
    raise RuntimeError("Results consolidation failed.")

from pathlib import Path
import shutil
import numpy as np

print("=== CYP3A4 STEP 94: DOCK–EM–DOCK COMPLEX PREPARATION ===")

root = Path("/content")
output_dir = root / "CYP3A4_DOCK_EM_DOCK"
output_dir.mkdir(exist_ok=True)

complexes = {
    "EE": "1TQN_ee_out_1.-10.6.complex.pdb",
    "LNG": "1TQN_lng_out_2.-8.5.complex.pdb",
    "NGMN": "1TQN_ngmn_out_1.-9.5.complex.pdb",
    "MPA": "1TQN_mpa_out_2.-8.0.complex.pdb",
}

for ligand, filename in complexes.items():

    source = root / filename

    if not source.exists():
        raise FileNotFoundError(source)

    destination = output_dir / f"CYP3A4_{ligand}_original_docked.pdb"

    shutil.copy2(source, destination)

    protein = []
    heme = []
    ligand_atoms = []

    fe = None
    sg = None

    for line in destination.read_text().splitlines():

        if not line.startswith(("ATOM  ", "HETATM")):
            continue

        residue = line[17:20].strip()
        resid = line[22:26].strip()
        atom = line[12:16].strip()

        xyz = np.array([
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        ])

        if residue == "HEM":
            heme.append(xyz)

            if atom == "FE":
                fe = xyz

        elif residue == "UNL":
            ligand_atoms.append(xyz)

        elif line.startswith("ATOM"):
            protein.append(xyz)

            if residue == "CYS" and resid == "442" and atom == "SG":
                sg = xyz

    if fe is None or sg is None:
        raise RuntimeError(
            f"{ligand}: Missing heme Fe or CYS442 SG"
        )
```

### Source lines 16249–16430
```python

print("\n=== STEP 108 COMPLETE ===")
print("Fe–SG bond and finite energy verified.")
print("No minimization or redocking performed.")
"""

result = subprocess.run(
    [python_exe, "-c", script],
    capture_output=True,
    text=True,
    timeout=180
)

print(result.stdout)

if result.stderr:
    print("\nDiagnostics:")
    print(result.stderr[-4000:])

print("Process return code:", result.returncode)

from pathlib import Path
import subprocess

print("=== CYP3A4 STEP 109: RESTRAINED RECEPTOR ENERGY MINIMIZATION ===")

python_exe = "/content/cyp3a4_openff_env/bin/python"

script = r"""
from pathlib import Path
import numpy as np
from openmm import (
    CustomExternalForce,
    LangevinMiddleIntegrator,
    Context,
    Platform,
    LocalEnergyMinimizer,
    unit
)
from openmm import app

work = Path("/content/CYP3A4_DOCK_EM_DOCK")

top_file = work / "CYP3A4_IC6_receptor.prmtop"
crd_file = work / "CYP3A4_IC6_receptor.inpcrd"

top = app.AmberPrmtopFile(str(top_file))
crd = app.AmberInpcrdFile(str(crd_file))

system = top.createSystem(
    nonbondedMethod=app.NoCutoff,
    constraints=None
)

positions = crd.positions

print("\n=== 1. ADD POSITION RESTRAINTS ===")

restraint = CustomExternalForce(
    "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
)

restraint.addGlobalParameter(
    "k",
    1000.0 * unit.kilojoule_per_mole / unit.nanometer**2
)

for parameter in ["x0", "y0", "z0"]:
    restraint.addPerParticleParameter(parameter)

xyz = positions.value_in_unit(unit.nanometer)

restrained = 0

for atom in top.topology.atoms():

    is_backbone = (
        atom.name in {"N", "CA", "C", "O"}
        and atom.residue.name not in {"HEM", "CYP"}
    )

    is_heme = (
        atom.residue.name == "HEM"
        and atom.element.symbol != "H"
    )

    if is_backbone or is_heme:
        restraint.addParticle(atom.index, xyz[atom.index])
        restrained += 1

system.addForce(restraint)

print("Restrained atoms:", restrained)
print("Total atoms:", system.getNumParticles())

print("\n=== 2. INITIAL ENERGY ===")

integrator = LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1 / unit.picosecond,
    0.002 * unit.picoseconds
)

platform = Platform.getPlatformByName("CPU")

context = Context(system, integrator, platform)
context.setPositions(positions)

initial = context.getState(getEnergy=True)
initial_energy = initial.getPotentialEnergy().value_in_unit(
    unit.kilojoule_per_mole
)

print(f"Initial energy: {initial_energy:.3f} kJ/mol")

print("\n=== 3. MINIMIZATION ===")

LocalEnergyMinimizer.minimize(
    context,
    tolerance=10,
    maxIterations=1000
)

final_state = context.getState(
    getEnergy=True,
    getPositions=True
)

final_energy = final_state.getPotentialEnergy().value_in_unit(
    unit.kilojoule_per_mole
)

print(f"Final energy: {final_energy:.3f} kJ/mol")
print(f"Energy change: {final_energy-initial_energy:.3f} kJ/mol")

if not np.isfinite(final_energy):
    raise RuntimeError("Nonfinite minimized energy.")

print("\n=== 4. SAVE MINIMIZED RECEPTOR ===")

output = work / "CYP3A4_IC6_receptor_minimized.pdb"

with output.open("w") as handle:
    app.PDBFile.writeFile(
        top.topology,
        final_state.getPositions(),
        handle,
        keepIds=True
    )

print("Saved:", output)

print("\n=== STEP 109 COMPLETE ===")
print("Receptor-only minimization completed.")
print("Full ligand-containing Dock–EM–Dock still pending.")
"""

result = subprocess.run(
    [python_exe, "-c", script],
    capture_output=True,
    text=True,
    timeout=600
)

print(result.stdout)

if result.stderr:
    print("\nDiagnostics:")
    print(result.stderr[-4000:])

print("Process return code:", result.returncode)

from pathlib import Path
import subprocess
import os

print("=== CYP3A4 STEP 110: GAFF2 LIGAND PARAMETERIZATION ===")

work = Path("/content/CYP3A4_DOCK_EM_DOCK")
amber = Path("/content/cyp3a4_amber_env")
ligand_dir = work / "GAFF2_ligands"
ligand_dir.mkdir(exist_ok=True)
```

### Source lines 16978–17167
```python

        integrator = LangevinMiddleIntegrator(
            300 * unit.kelvin,
            1 / unit.picosecond,
            0.002 * unit.picoseconds
        )

        context = Context(system, integrator, platform)
        context.setPositions(positions)

        initial_energy = (
            context.getState(getEnergy=True)
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
        )

        print("Total atoms:", system.getNumParticles(), flush=True)
        print("Ligand heavy atoms:", len(ligand_heavy), flush=True)
        print("Restrained atoms:", restrained, flush=True)
        print(f"Initial energy: {initial_energy:.3f} kJ/mol", flush=True)

        if not np.isfinite(initial_energy):
            raise RuntimeError("Initial energy is nonfinite")

        LocalEnergyMinimizer.minimize(
            context,
            tolerance=10,
            maxIterations=1000
        )

        state = context.getState(
            getEnergy=True,
            getPositions=True
        )

        final_energy = (
            state.getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
        )

        xyz_final = np.array(
            state.getPositions().value_in_unit(unit.nanometer)
        )

        heavy_indices = [a.index for a in ligand_heavy]

        ligand_rmsd = np.sqrt(
            np.mean(
                np.sum(
                    (
                        xyz_final[heavy_indices]
                        - xyz_initial[heavy_indices]
                    )**2,
                    axis=1
                )
            )
        ) * 10.0

        fe_sg_distance = np.linalg.norm(
            xyz_final[fe.index] - xyz_final[sg.index]
        ) * 10.0

        print(f"Final energy: {final_energy:.3f} kJ/mol", flush=True)
        print(
            f"Energy change: {final_energy-initial_energy:.3f} kJ/mol",
            flush=True
        )
        print(
            f"Ligand heavy-atom displacement RMSD: {ligand_rmsd:.3f} Å",
            flush=True
        )
        print(
            f"Final Fe–SG distance: {fe_sg_distance:.3f} Å",
            flush=True
        )

        if not np.isfinite(final_energy):
            raise RuntimeError("Final energy is nonfinite")

        if not np.all(np.isfinite(xyz_final)):
            raise RuntimeError("Nonfinite minimized coordinates")

        with out_path.open("w") as handle:
            app.PDBFile.writeFile(
                topology,
                state.getPositions(),
                handle,
                keepIds=True
            )

        print("Saved:", out_path, flush=True)
        print("STATUS: SUCCESS", flush=True)

        results[name] = {
            "initial": initial_energy,
            "final": final_energy,
            "rmsd": ligand_rmsd,
            "fe_sg": fe_sg_distance
        }

    except Exception as exc:
        print("STATUS: FAILED", flush=True)
        print(type(exc).__name__ + ":", exc, flush=True)
        results[name] = None

    finally:
        for variable in [
            "context", "integrator", "system",
            "top", "crd", "state"
        ]:
            if variable in locals():
                del locals()[variable]
        gc.collect()

print("\n=== STEP 113 SUMMARY ===", flush=True)

for name in names:
    result = results.get(name)

    if result is None:
        print(f"{name}: FAILED", flush=True)
    else:
        print(
            f"{name}: SUCCESS | "
            f"Final energy {result['final']:.2f} kJ/mol | "
            f"Ligand RMSD {result['rmsd']:.3f} Å | "
            f"Fe–SG {result['fe_sg']:.3f} Å",
            flush=True
        )

print("\n=== STEP 113 COMPLETE ===", flush=True)
print("Full-complex minimization attempted for all four ligands.", flush=True)
print("Redocking not yet performed.", flush=True)
"""

result = subprocess.run(
    [python_exe, "-u", "-c", script],
    capture_output=True,
    text=True,
    timeout=1800
)

print(result.stdout)

if result.stderr:
    print("\nDiagnostics:")
    print(result.stderr[-5000:])

print("Process return code:", result.returncode)

from pathlib import Path
import subprocess
import shutil

print("=== CYP3A4 STEP 114: PREPARE HEME-INCLUSIVE REDOCKING ===")

work = Path("/content/CYP3A4_DOCK_EM_DOCK")
outdir = work / "redocking"
outdir.mkdir(exist_ok=True)

names = ["EE", "LNG", "NGMN", "MPA"]

print("\n=== AVAILABLE PDBQT CONVERSION TOOLS ===")

candidates = [
    "/content/cyp3a4_amber_env/bin/obabel",
    "/content/cyp3a4_openff_env/bin/obabel",
    "/usr/bin/obabel",
    "/content/cyp3a4_openff_env/bin/mk_prepare_receptor.py",
    "/content/cyp3a4_openff_env/bin/mk_prepare_ligand.py",
]

for candidate in candidates:
    path = Path(candidate)
    print(f"{candidate}: {'AVAILABLE' if path.exists() else 'NOT FOUND'}")

for command in ["obabel", "prepare_receptor4.py", "mk_prepare_receptor.py"]:
    print(f"PATH lookup {command}: {shutil.which(command)}")

print("\n=== EXTRACT MINIMIZED RECEPTORS AND LIGANDS ===")

for name in names:
    print(f"\n--- {name} ---")

    source = work / f"CYP3A4_{name}_full_complex_minimized.pdb"

    if not source.exists():
        print("MISSING:", source)
        continue

```

### Source lines 17550–17638
```python
        raise RuntimeError(
            f"{name}: heme validation failed"
        )

    output.write_text("\n".join(lines) + "\n")

    print("Saved:", output, flush=True)

    # Confirm iron is physically present in the output.
    iron_lines = [
        line for line in lines
        if line[12:16].strip().upper() == "FE"
        and line[17:20].strip() == "HEM"
    ]

    assert len(iron_lines) == 1

    print("Iron record:", iron_lines[0], flush=True)

    print("STATUS: PDBQT WRITTEN", flush=True)

    try:
        from vina import Vina

        v = Vina(sf_name="vina", cpu=1)
        v.set_receptor(str(output))

        print("Vina receptor loading: SUCCESS", flush=True)

    except Exception as exc:
        print(
            "Vina receptor loading: FAILED",
            type(exc).__name__,
            str(exc)[:500],
            flush=True
        )

print("\n=== STEP 116 COMPLETE ===", flush=True)
print("No redocking performed.", flush=True)
"""

result = subprocess.run(
    [python_exe, "-u", "-c", script],
    capture_output=True,
    text=True,
    timeout=300
)

print(result.stdout)

if result.stderr:
    print("\nDiagnostics:")
    print(result.stderr[-4000:])

print("Process return code:", result.returncode)

from pathlib import Path
import subprocess

print("=== CYP3A4 STEP 117: HEME-INCLUSIVE REDOCKING ===")

python_exe = "/content/cyp3a4_openff_env/bin/python"

script = r"""
from pathlib import Path
import numpy as np
from rdkit import Chem
from vina import Vina

work = Path("/content/CYP3A4_DOCK_EM_DOCK")
redock_dir = work / "redocking"
redock_dir.mkdir(exist_ok=True)

names = ["EE", "LNG", "NGMN", "MPA"]

results = {}

for name in names:

    print(f"\n{'='*60}", flush=True)
    print(f"REDOCKING CYP3A4–HEME–{name}", flush=True)
    print(f"{'='*60}", flush=True)

    receptor = (
        redock_dir /
        f"CYP3A4_{name}_minimized_receptor_heme.pdbqt"
    )

    ligand = Path(
```

### Source lines 17669–17757
```python

        heavy_indices = [
            atom.GetIdx()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() > 1
        ]

        heavy_xyz = xyz[heavy_indices]

        center = heavy_xyz.mean(axis=0)

        ligand_span = (
            heavy_xyz.max(axis=0)
            - heavy_xyz.min(axis=0)
        )

        box_size = np.maximum(
            ligand_span + 16.0,
            24.0
        )

        print("Search center:", center.round(3).tolist(), flush=True)
        print("Box size:", box_size.round(3).tolist(), flush=True)

        v = Vina(
            sf_name="vina",
            cpu=2,
            seed=42,
            verbosity=0
        )

        v.set_receptor(str(receptor))
        v.set_ligand_from_file(str(ligand))

        v.compute_vina_maps(
            center=center.tolist(),
            box_size=box_size.tolist()
        )

        v.dock(
            exhaustiveness=16,
            n_poses=10
        )

        energies = v.energies()

        if len(energies) == 0:
            raise RuntimeError("Vina returned no poses")

        v.write_poses(
            str(output),
            n_poses=10,
            overwrite=True
        )

        best_score = float(energies[0][0])

        print(
            f"Best redocking score: {best_score:.3f} kcal/mol",
            flush=True
        )

        print(
            "Number of poses:", len(energies),
            flush=True
        )

        print("Saved:", output, flush=True)

        results[name] = best_score

        print("STATUS: SUCCESS", flush=True)

    except Exception as exc:
        print("STATUS: FAILED", flush=True)
        print(
            type(exc).__name__ + ":",
            str(exc)[:1000],
            flush=True
        )

        results[name] = None

print("\n=== STEP 117 REDOCKING SUMMARY ===", flush=True)

for name in names:
    score = results[name]

    if score is None:
```

### Source lines 18170–18288
```python
        original_heavy = [
            a.GetIdx()
            for a in mol.GetAtoms()
            if a.GetAtomicNum() > 1
        ]

        original_xyz = np.array(
            mol.GetConformer().GetPositions()
        )[original_heavy]

        minimized_xyz = minimized_ligand_xyz(minimized)

        redocked_lines = read_first_pdbqt_model(redocked)
        redocked_xyz = pdbqt_heavy_xyz(redocked_lines)

        if not (
            len(original_xyz)
            == len(minimized_xyz)
            == len(redocked_xyz)
        ):
            raise RuntimeError(
                "Heavy-atom count mismatch"
            )

        # Coordinate RMSD without superposition.
        # Meaningful only if atom order is consistent.
        direct_rmsd = np.sqrt(
            np.mean(
                np.sum(
                    (redocked_xyz - minimized_xyz)**2,
                    axis=1
                )
            )
        )

        centroid_shift = np.linalg.norm(
            redocked_xyz.mean(axis=0)
            - minimized_xyz.mean(axis=0)
        )

        iron = get_iron(minimized)

        fe_distance = np.min(
            np.linalg.norm(
                redocked_xyz - iron,
                axis=1
            )
        )

        print(
            f"Heavy atoms: {len(redocked_xyz)}"
        )

        print(
            f"Direct coordinate RMSD: "
            f"{direct_rmsd:.3f} Å"
        )

        print(
            f"Centroid shift: "
            f"{centroid_shift:.3f} Å"
        )

        print(
            f"Nearest ligand–Fe distance: "
            f"{fe_distance:.3f} Å"
        )

        if direct_rmsd > 3.0:
            print(
                "FLAG: Substantial pose movement"
            )

        rows.append({
            "Ligand": name,
            "Heavy_atoms": len(redocked_xyz),
            "Direct_RMSD_A": round(direct_rmsd, 3),
            "Centroid_shift_A": round(centroid_shift, 3),
            "Nearest_Fe_distance_A": round(fe_distance, 3),
            "Status": "CHECK" if direct_rmsd > 3 else "PASS"
        })

        print("STATUS: COMPLETE")

    except Exception as exc:
        print(
            "STATUS: FAILED",
            type(exc).__name__,
            str(exc)
        )

        rows.append({
            "Ligand": name,
            "Heavy_atoms": "",
            "Direct_RMSD_A": "",
            "Centroid_shift_A": "",
            "Nearest_Fe_distance_A": "",
            "Status": "FAILED"
        })

output = redock / "CYP3A4_final_pose_validation.csv"

with output.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "Ligand",
            "Heavy_atoms",
            "Direct_RMSD_A",
            "Centroid_shift_A",
            "Nearest_Fe_distance_A",
            "Status"
        ]
    )

    writer.writeheader()
    writer.writerows(rows)

print("\n=== STEP 119 SUMMARY ===")
```

### Source lines 18403–18597
```python
min_xyz = heavy_coordinates(minimized)
red_xyz = heavy_coordinates(redocked)

print("Reference heavy atoms:", len(ref_xyz))
print("Minimized heavy atoms:", len(min_xyz))
print("Redocked heavy atoms:", len(red_xyz))

if len(ref_xyz) != len(min_xyz) or len(ref_xyz) != len(red_xyz):
    raise RuntimeError("Heavy-atom count mismatch")

print("\n=== COORDINATE COMPARISON ===")

direct_rmsd = np.sqrt(
    np.mean(
        np.sum((red_xyz - min_xyz)**2, axis=1)
    )
)

centroid_shift = np.linalg.norm(
    red_xyz.mean(axis=0)
    - min_xyz.mean(axis=0)
)

print(
    f"Direct coordinate RMSD: {direct_rmsd:.3f} Å"
)

print(
    f"Centroid shift: {centroid_shift:.3f} Å"
)

# RDKit symmetry-aware comparison after alignment.
# This removes overall translation and rotation, so it is
# a conformational comparison, not a binding-pose RMSD.
try:
    reference_noh = Chem.RemoveHs(reference)
    redocked_noh = Chem.RemoveHs(redocked)

    symmetry_rmsd = rdMolAlign.GetBestRMS(
        Chem.Mol(redocked_noh),
        Chem.Mol(reference_noh)
    )

    print(
        f"Symmetry-aware aligned RMSD "
        f"(redocked vs original): {symmetry_rmsd:.3f} Å"
    )

except Exception as exc:
    print(
        "Symmetry-aware comparison unavailable:",
        str(exc)[:400]
    )

print("\n=== INTERPRETATION ===")

if centroid_shift > 3:
    print(
        "NGMN changed binding-site position substantially."
    )
    print(
        "This cannot be explained by atom ordering alone, "
        "because centroid position is independent of atom order."
    )
else:
    print(
        "NGMN centroid remained close to the minimized position."
    )

print("\n=== STEP 120 COMPLETE ===")
"""

result = subprocess.run(
    [python_exe, "-u", "-c", script],
    capture_output=True,
    text=True,
    timeout=300
)

print(result.stdout)

if result.stderr:
    print("\nDiagnostics:")
    print(result.stderr[-5000:])

print("Process return code:", result.returncode)

from pathlib import Path
import csv

print("=== CYP3A4 STEP 121: CONSOLIDATED DOCK–EM–DOCK RESULTS ===")

work = Path("/content/CYP3A4_DOCK_EM_DOCK")
redock = work / "redocking"

data = [
    {
        "Ligand": "EE",
        "Original_CB_Dock3_kcal_mol": -10.6,
        "Initial_EM_energy_kJ_mol": 466185.766,
        "Final_EM_energy_kJ_mol": -45881.115,
        "EM_ligand_displacement_A": 0.560,
        "Redocking_Vina_kcal_mol": -10.898,
        "Direct_pose_RMSD_A": 0.697,
        "Centroid_shift_A": 0.547,
        "Nearest_Fe_distance_A": 9.524,
        "Pose_interpretation": "Small positional change"
    },
    {
        "Ligand": "LNG",
        "Original_CB_Dock3_kcal_mol": -8.5,
        "Initial_EM_energy_kJ_mol": 466383.984,
        "Final_EM_energy_kJ_mol": -45844.797,
        "EM_ligand_displacement_A": 0.710,
        "Redocking_Vina_kcal_mol": -7.778,
        "Direct_pose_RMSD_A": 1.123,
        "Centroid_shift_A": 0.583,
        "Nearest_Fe_distance_A": 21.162,
        "Pose_interpretation": "Small centroid displacement"
    },
    {
        "Ligand": "NGMN",
        "Original_CB_Dock3_kcal_mol": -9.5,
        "Initial_EM_energy_kJ_mol": 466319.132,
        "Final_EM_energy_kJ_mol": -45811.024,
        "EM_ligand_displacement_A": 0.590,
        "Redocking_Vina_kcal_mol": -10.430,
        "Direct_pose_RMSD_A": 7.025,
        "Centroid_shift_A": 3.904,
        "Nearest_Fe_distance_A": 8.780,
        "Pose_interpretation": "Repositioned binding pose; aligned RMSD 0.100 A"
    },
    {
        "Ligand": "MPA",
        "Original_CB_Dock3_kcal_mol": -8.0,
        "Initial_EM_energy_kJ_mol": 467032.296,
        "Final_EM_energy_kJ_mol": -45881.554,
        "EM_ligand_displacement_A": 0.625,
        "Redocking_Vina_kcal_mol": -7.597,
        "Direct_pose_RMSD_A": 1.850,
        "Centroid_shift_A": 0.209,
        "Nearest_Fe_distance_A": 21.240,
        "Pose_interpretation": "Small centroid displacement"
    }
]

for row in data:
    name = row["Ligand"]

    receptor = redock / (
        f"CYP3A4_{name}_minimized_receptor_heme.pdbqt"
    )

    poses = redock / (
        f"CYP3A4_{name}_heme_redocked.pdbqt"
    )

    minimized = work / (
        f"CYP3A4_{name}_full_complex_minimized.pdb"
    )

    if not all(p.exists() for p in [receptor, poses, minimized]):
        raise FileNotFoundError(
            f"Missing CYP3A4 result file for {name}"
        )

    receptor_lines = receptor.read_text().splitlines()

    heme_lines = [
        line for line in receptor_lines
        if line.startswith(("ATOM", "HETATM"))
        and line[17:20].strip() == "HEM"
    ]

    iron_lines = [
        line for line in heme_lines
        if line[12:16].strip().upper() == "FE"
    ]

    if len(heme_lines) != 43 or len(iron_lines) != 1:
        raise RuntimeError(
            f"{name}: heme/iron validation failed"
        )

    row["Heme_included"] = "YES"
    row["Full_complex_EM"] = "YES"
    row["Redocking_receptor"] = "Minimized protein-heme complex"
    row["Redocking_ligand_input"] = "Previously prepared ligand PDBQT"

    print(
        f"{name}: verified | "
        f"Vina {row['Redocking_Vina_kcal_mol']:.3f} kcal/mol | "
        f"Heme retained"
    )

```

### Source lines 18610–18769
```python

from pathlib import Path
import subprocess

print("=== CYP3A4 STEP 122: TABLE 4 + SEQUENTIAL FILE INVENTORY ===")

script = r"""
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from openmm import app, unit

work = Path("/content/CYP3A4_DOCK_EM_DOCK")
names = ["EE", "LNG", "NGMN", "MPA"]

print("\n=== TABLE 4: STRUCTURAL DISPLACEMENT ===")
print("Pocket definition: protein residues within 5 Å of original ligand heavy atoms")
print("RMSDs: direct coordinate comparisons without receptor superposition\n")

def atoms_and_xyz(topology, positions):
    atoms = list(topology.atoms())
    xyz = np.asarray(positions.value_in_unit(unit.angstrom))
    return atoms, xyz

def rmsd(a, b):
    return float(np.sqrt(np.mean(np.sum((a-b)**2, axis=1))))

for name in names:
    top_path = work / f"CYP3A4_{name}_full_complex.prmtop"
    crd_path = work / f"CYP3A4_{name}_full_complex.inpcrd"
    min_path = work / f"CYP3A4_{name}_full_complex_minimized.pdb"

    try:
        top = app.AmberPrmtopFile(str(top_path))
        crd = app.AmberInpcrdFile(str(crd_path))
        minimized = app.PDBFile(str(min_path))

        atoms, before = atoms_and_xyz(top.topology, crd.positions)
        after = np.asarray(
            minimized.positions.value_in_unit(unit.angstrom)
        )

        if len(atoms) != len(after):
            raise ValueError("Atom count mismatch")

        protein_residues = {
            r.index for r in top.topology.residues()
            if r.name not in {"HEM", "LIG"}
        }

        protein_heavy = [
            a.index for a in atoms
            if a.residue.index in protein_residues
            and a.element is not None
            and a.element.symbol != "H"
        ]

        ca = [
            a.index for a in atoms
            if a.residue.index in protein_residues
            and a.name == "CA"
        ]

        ligand = [
            a.index for a in atoms
            if a.residue.name == "LIG"
            and a.element is not None
            and a.element.symbol != "H"
        ]

        distances = cdist(before[protein_heavy], before[ligand])

        pocket_residues = {
            atoms[protein_heavy[i]].residue.index
            for i in np.where(np.any(distances <= 5.0, axis=1))[0]
        }

        pocket_ca = [
            i for i in ca
            if atoms[i].residue.index in pocket_residues
        ]

        pocket_heavy = [
            i for i in protein_heavy
            if atoms[i].residue.index in pocket_residues
        ]

        centroid_shift = np.linalg.norm(
            after[ligand].mean(axis=0)
            - before[ligand].mean(axis=0)
        )

        print(
            f"{name}: whole Cα={rmsd(before[ca],after[ca]):.4f} Å | "
            f"pocket Cα={rmsd(before[pocket_ca],after[pocket_ca]):.4f} Å | "
            f"pocket heavy={rmsd(before[pocket_heavy],after[pocket_heavy]):.4f} Å | "
            f"ligand centroid={centroid_shift:.4f} Å | "
            f"pocket residues={len(pocket_residues)}",
            flush=True
        )

    except Exception as e:
        print(f"{name}: ERROR: {type(e).__name__}: {e}", flush=True)

print("\n=== INVENTORY: POTENTIAL SEQUENTIAL STRUCTURES ===")

roots = [
    Path("/content/CYP3A4_DOCK_EM_DOCK"),
    Path("/content")
]

extensions = {".pdb", ".pdbqt", ".sdf", ".mol2", ".csv", ".zip"}

found = []

for root in roots:
    if not root.exists():
        continue

    for path in root.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in extensions:
            continue

        name = path.name.upper()

        # Avoid duplicate listings of files inside the working directory.
        if root == Path("/content") and work in path.parents:
            continue

        pair_tokens = [
            "EE_LNG", "EE_NGMN", "EE_MPA",
            "LNG_EE", "LNG_NGMN", "LNG_MPA",
            "NGMN_EE", "NGMN_LNG", "NGMN_MPA",
            "MPA_EE", "MPA_LNG", "MPA_NGMN"
        ]

        relevant = (
            any(token in name for token in pair_tokens)
            or ("CYP" in name and any(
                token in name
                for token in ["SEQ", "PAIR", "BACKUP", "RESULT"]
            ))
        )

        if relevant:
            found.append(path)

print("Potential sequential files found:", len(found))

for path in sorted(set(found))[:100]:
    print(path)

if len(found) > 100:
    print(f"... and {len(found)-100} additional files")

print("\n=== STEP 122 COMPLETE ===")
"""

result = subprocess.run(
    ["/content/cyp3a4_openff_env/bin/python", "-u", "-c", script],
```

### Source lines 19019–19107
```python

    output = (
        seq_dir /
        f"CYP3A4_{ligand_a}_receptor_heme_ligandA.pdbqt"
    )

    output.write_text(
        "".join(line + "\n" for line in receptor_records)
        + "".join(ligand_records)
    )

    built[ligand_a] = output

    print(
        f"{ligand_a}: "
        f"{len(receptor_records)} receptor/heme heavy atoms + "
        f"{len(ligand_records)} Ligand A heavy atoms | "
        f"Saved: {output.name}"
    )

df["sequential_receptor_pdbqt"] = df["ligand_A"].map(
    lambda a: str(built[a])
)

df.to_csv(manifest_path, index=False)

print("\nSequential combinations:", len(df))
print("Sequential receptor files:", len(built))
print("Heme retained in original receptor records.")
print("\n=== STEP 124 COMPLETE ===")

from pathlib import Path
import pandas as pd
import numpy as np
from vina import Vina
from openmm import app, unit

print("=== CYP3A4 STEP 125: DOCK ALL 12 SEQUENTIAL COMBINATIONS ===")

work = Path("/content/CYP3A4_DOCK_EM_DOCK")
seq_dir = work / "sequential"
manifest_path = seq_dir / "CYP3A4_12_sequential_manifest.csv"

df = pd.read_csv(manifest_path)

if len(df) != 12:
    raise RuntimeError("Expected exactly 12 combinations")

# Search-box dimensions, in Å.
# These boxes are centered on the original Ligand B docking pose.
BOX_SIZE = [26.0, 26.0, 26.0]

EXHAUSTIVENESS = 16
SEED = 42
CPU = 2

def get_ligand_centroid(ligand):
    """
    Obtain the heavy-atom centroid of the original ligand pose.
    This defines the center of the sequential docking search box.
    """
    sdf_candidates = list(
        Path("/content").glob(
            f"*CYP3A4*{ligand}*.sdf"
        )
    )

    if not sdf_candidates:
        sdf_candidates = list(
            work.rglob(f"*{ligand}*.sdf")
        )

    from rdkit import Chem

    for path in sdf_candidates:
        supplier = Chem.SDMolSupplier(
            str(path),
            removeHs=False
        )

        mol = next(
            (m for m in supplier if m is not None),
            None
        )

        if mol is None or not mol.GetNumConformers():
            continue

        conf = mol.GetConformer()
```

### Source lines 19128–19427
```python
        [round(v, 3) for v in centers[ligand]]
    )

results = []

for idx, row in df.iterrows():
    a = row["ligand_A"]
    b = row["ligand_B"]

    sequence = f"{a}_to_{b}"

    receptor = Path(row["sequential_receptor_pdbqt"])
    ligand = Path(row["ligand_B_pdbqt"])

    if not receptor.exists():
        raise FileNotFoundError(receptor)

    if not ligand.exists():
        raise FileNotFoundError(ligand)

    output = seq_dir / f"CYP3A4_{sequence}_sequential_redocked.pdbqt"

    print(f"\n[{idx+1}/12] Docking {a} -> {b}", flush=True)

    v = Vina(
        sf_name="vina",
        cpu=CPU,
        seed=SEED,
        verbosity=0
    )

    v.set_receptor(str(receptor))
    v.set_ligand_from_file(str(ligand))

    v.compute_vina_maps(
        center=centers[b],
        box_size=BOX_SIZE
    )

    v.dock(
        exhaustiveness=EXHAUSTIVENESS,
        n_poses=10
    )

    energies = v.energies()

    if len(energies) == 0:
        raise RuntimeError(f"No docking poses for {sequence}")

    best_score = float(energies[0][0])

    v.write_poses(
        str(output),
        n_poses=10,
        overwrite=True
    )

    original = float(
        row["original_sequential_score_kcal_mol"]
    )

    delta = best_score - original

    if delta < 0:
        interpretation = "More negative numerical score"
    elif delta > 0:
        interpretation = "Less negative numerical score"
    else:
        interpretation = "Unchanged numerical score"

    results.append({
        "sequence": f"{a} -> {b}",
        "ligand_A": a,
        "ligand_B": b,
        "original_rigid_sequential_score": original,
        "post_EM_sequential_score": best_score,
        "delta_score": delta,
        "interpretation": interpretation,
        "pose_file": str(output)
    })

    print(
        f"Original: {original:.3f} | "
        f"Post-EM: {best_score:.3f} | "
        f"Delta: {delta:+.3f} kcal/mol",
        flush=True
    )

results_df = pd.DataFrame(results)

output_csv = (
    seq_dir /
    "CYP3A4_12_post_EM_sequential_docking_scores.csv"
)

results_df.to_csv(output_csv, index=False)

print("\n=== TABLE 5: SEQUENTIAL DOCKING SCORES ===")
print(
    results_df[
        [
            "sequence",
            "original_rigid_sequential_score",
            "post_EM_sequential_score",
            "delta_score",
            "interpretation"
        ]
    ].to_string(index=False)
)

print("\nSaved:", output_csv)
print("Completed:", len(results_df), "/ 12")
print("\n=== STEP 125 COMPLETE ===")

import sys
import subprocess
from pathlib import Path

print("=== CYP3A4 STEP 125A: FIX VINA ENVIRONMENT ===")

env_python = Path("/content/cyp3a4_openff_env/bin/python")

if env_python.exists():
    print("Previous CYP3A4 environment found.")

    result = subprocess.run(
        [str(env_python), "-c",
         "from vina import Vina; print('Vina import successful')"],
        capture_output=True,
        text=True
    )

    if result.returncode == 0:
        print(result.stdout.strip())
        print("VINA_PYTHON =", env_python)
        print("No installation required.")
    else:
        print("Vina missing from previous environment.")
        print("Installing Vina there...")

        subprocess.check_call([
            str(env_python), "-m", "pip",
            "install", "vina"
        ])

        subprocess.check_call([
            str(env_python), "-c",
            "from vina import Vina; print('Vina import successful')"
        ])

else:
    print("Previous environment not found.")
    print("Installing Vina in the current Colab environment...")

    subprocess.check_call([
        sys.executable, "-m", "pip",
        "install", "vina"
    ])

    subprocess.check_call([
        sys.executable, "-c",
        "from vina import Vina; print('Vina import successful')"
    ])

print("\n=== STEP 125A COMPLETE ===")

import subprocess
from pathlib import Path

print("=== CYP3A4 STEP 125B: 12 SEQUENTIAL DOCKINGS ===", flush=True)

python_exe = "/content/cyp3a4_openff_env/bin/python"

script = r'''
from pathlib import Path
import numpy as np
import pandas as pd
from vina import Vina
from openmm import app, unit

work = Path("/content/CYP3A4_DOCK_EM_DOCK")
seq = work / "sequential"

manifest = pd.read_csv(
    seq / "CYP3A4_12_sequential_manifest.csv"
)

assert len(manifest) == 12

# Original Ligand B positions define the search centers.
centers = {}

for ligand in ["EE", "LNG", "NGMN", "MPA"]:
    top_path = work / f"CYP3A4_{ligand}_full_complex.prmtop"
    crd_path = work / f"CYP3A4_{ligand}_full_complex.inpcrd"

    topology = app.AmberPrmtopFile(str(top_path))
    coordinates = app.AmberInpcrdFile(str(crd_path))

    xyz = np.asarray(
        coordinates.positions.value_in_unit(unit.angstrom)
    )

    indices = [
        atom.index
        for atom in topology.topology.atoms()
        if atom.residue.name == "LIG"
        and atom.element is not None
        and atom.element.symbol != "H"
    ]

    if not indices:
        raise RuntimeError(f"{ligand}: original ligand not found")

    centers[ligand] = xyz[indices].mean(axis=0).tolist()

    print(
        f"{ligand} original pose center:",
        [round(x, 3) for x in centers[ligand]],
        flush=True
    )

results_path = (
    seq / "CYP3A4_12_post_EM_sequential_docking_scores.csv"
)

results = []

for number, row in enumerate(manifest.itertuples(index=False), 1):

    a = row.ligand_A
    b = row.ligand_B
    sequence = f"{a}_to_{b}"

    receptor = Path(row.sequential_receptor_pdbqt)
    ligand_b = Path(row.ligand_B_pdbqt)

    pose_path = seq / f"CYP3A4_{sequence}_sequential_redocked.pdbqt"

    if not receptor.exists():
        raise FileNotFoundError(receptor)

    if not ligand_b.exists():
        raise FileNotFoundError(ligand_b)

    print(
        f"\n[{number}/12] {a} -> {b}",
        flush=True
    )

    vina = Vina(
        sf_name="vina",
        cpu=2,
        seed=42,
        verbosity=0
    )

    vina.set_receptor(str(receptor))
    vina.set_ligand_from_file(str(ligand_b))

    vina.compute_vina_maps(
        center=centers[b],
        box_size=[26.0, 26.0, 26.0]
    )

    vina.dock(
        exhaustiveness=16,
        n_poses=10
    )

    energies = vina.energies()

    if len(energies) == 0:
        raise RuntimeError(f"{sequence}: no docking poses")

    score = float(energies[0][0])

    vina.write_poses(
        str(pose_path),
        n_poses=10,
        overwrite=True
    )

    original = float(row.original_sequential_score_kcal_mol)
    delta = score - original

    results.append({
        "sequence": f"{a} -> {b}",
        "ligand_A": a,
        "ligand_B": b,
        "original_rigid_sequential_score": original,
        "post_EM_sequential_score": score,
        "delta_score": delta,
        "interpretation": (
            "More negative numerical score" if delta < -0.0005
            else "Less negative numerical score" if delta > 0.0005
            else "Unchanged numerical score"
        ),
        "pose_file": str(pose_path)
    })
```

### Source lines 19428–19493
```python

    pd.DataFrame(results).to_csv(results_path, index=False)

    print(
        f"Score: {score:.3f} kcal/mol | "
        f"Original: {original:.3f} | "
        f"Difference: {delta:+.3f}",
        flush=True
    )

print("\n=== TABLE 5: ALL 12 SEQUENCES ===")

df = pd.DataFrame(results)

print(
    df[
        [
            "sequence",
            "original_rigid_sequential_score",
            "post_EM_sequential_score",
            "delta_score",
            "interpretation"
        ]
    ].to_string(index=False)
)

print("\nSaved:", results_path)
print("Completed:", len(results), "/ 12")
print("=== STEP 125B COMPLETE ===")
'''

result = subprocess.run(
    [python_exe, "-u", "-c", script],
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT
)

print(result.stdout, flush=True)

if result.returncode != 0:
    raise RuntimeError(
        f"Sequential docking failed (exit code {result.returncode})"
    )

import subprocess
from pathlib import Path

print("=== CYP3A4 STEP 126: STRUCTURAL ANALYSIS OF 12 PAIRS ===", flush=True)

script = r'''
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from openmm import app, unit

work = Path("/content/CYP3A4_DOCK_EM_DOCK")
seq = work / "sequential"

manifest = pd.read_csv(
    seq / "CYP3A4_12_sequential_manifest.csv"
)

scores = pd.read_csv(
    seq / "CYP3A4_12_post_EM_sequential_docking_scores.csv"
```

### Source lines 19665–19767
```python
        "sequence": sequence,
        "ligand_A_location": "Original Ligand A pocket",
        "ligand_B_location": b_location,
        "observed_binding_configuration": configuration
    })

    results8.append({
        "sequence": sequence,
        "X_A": round(float(centroid_b[0]), 3),
        "Y_A": round(float(centroid_b[1]), 3),
        "Z_A": round(float(centroid_b[2]), 3),
        "heavy_atoms": len(ligand_b)
    })

    results9.append({
        "sequence": sequence,
        "centroid_distance_A": round(centroid_distance, 4),
        "minimum_heavy_atom_distance_A": round(min_distance, 4)
    })

    print(
        f"[{number}/12] {sequence} | "
        f"A pocket contacts={len(a_pocket)} | "
        f"B pocket contacts={len(b_pocket)} | "
        f"centroid distance={centroid_distance:.3f} Å | "
        f"minimum distance={min_distance:.3f} Å",
        flush=True
    )

tables = {
    6: pd.DataFrame(results6),
    7: pd.DataFrame(results7),
    8: pd.DataFrame(results8),
    9: pd.DataFrame(results9)
}

for number, df in tables.items():
    path = seq / f"CYP3A4_Table_{number}.csv"
    df.to_csv(path, index=False)

    print(f"\n=== TABLE {number} ===")
    print(df.to_string(index=False))
    print("Saved:", path)

print("\n=== STEP 126 COMPLETE ===")
'''

result = subprocess.run(
    ["/content/cyp3a4_openff_env/bin/python", "-u", "-c", script],
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    text=True
)

print(result.stdout, flush=True)

if result.returncode != 0:
    raise RuntimeError(
        f"Step 126 failed with exit code {result.returncode}"
    )

import sys
import subprocess
from pathlib import Path

print("=== CYP3A4 STEP 126A: ENVIRONMENT RECOVERY ===")

work = Path("/content/CYP3A4_DOCK_EM_DOCK")
seq = work / "sequential"

old_python = Path("/content/cyp3a4_openff_env/bin/python")

print("\n1. PYTHON ENVIRONMENT")
print("Current Python:", sys.executable)
print("Previous environment exists:", old_python.exists())

print("\n2. PROJECT FILES")
print("Project directory exists:", work.exists())
print("Sequential directory exists:", seq.exists())

manifest = seq / "CYP3A4_12_sequential_manifest.csv"
scores = seq / "CYP3A4_12_post_EM_sequential_docking_scores.csv"

print("12-combination manifest exists:", manifest.exists())
print("Sequential docking scores exist:", scores.exists())

print("\n3. DOCKED POSES")

ligands = ["EE", "LNG", "NGMN", "MPA"]

existing_poses = []
missing_poses = []

for a in ligands:
    for b in ligands:
        if a == b:
            continue

        pose = seq / f"CYP3A4_{a}_to_{b}_sequential_redocked.pdbqt"

        if pose.exists():
            existing_poses.append(pose)
            print(f"FOUND: {a} -> {b}")
```

### Source lines 26467–26619
```python
    report[ligand] = {
        "protein_atoms": len(protein_lines),
        "heme_heavy_atoms": len(heme_lines),
        "ligand_A_heavy_atoms": len(ligand_lines),
        "file": str(output),
        "charge_model": "zero charge for heme and Ligand A",
        "status": "provisional reconstructed receptor"
    }

    print(
        f"{ligand}: protein={len(protein_lines)}, "
        f"heme=43, Ligand A={len(ligand_lines)} — SAVED"
    )

report_path = OUT / "Step178_Receptor_Report.json"
report_path.write_text(json.dumps(report, indent=2))
shutil.copy2(report_path, DRIVE / report_path.name)

print("\n=== STEP 178 COMPLETE ===")
print("Four reconstructed receptors saved to Google Drive.")

from pathlib import Path
import subprocess, shutil, os, csv, json

print("=== STEP 179: RUN 12 SEQUENTIAL DOCKINGS ===")

BASE = Path("/content/CYP3A4_DOCK_EM_DOCK")
ENV = Path("/content/cyp3a4_parameterization_env")
OUT = BASE / "sequential/reconstructed_poses"
DRIVE = Path("/content/drive/MyDrive/Dock_EM_Dock_Research/CYP3A4/sequential/reconstructed_poses")
OUT.mkdir(parents=True, exist_ok=True)
DRIVE.mkdir(parents=True, exist_ok=True)

RECEPTORS = BASE / "sequential/reconstructed_receptors"
LIGANDS = BASE / "sequential/ligand_pdbqt"

# Historical box centers for Ligand B.
CENTERS = {
    "EE": (-24.301, -12.517, -13.995),
    "LNG": (-23.900, -27.553, -34.601),
    "NGMN": (-21.003, -15.678, -11.734),
    "MPA": (-24.652, -27.063, -33.990)
}

script = r'''
import csv, json, shutil
from pathlib import Path
from vina import Vina

BASE = Path("/content/CYP3A4_DOCK_EM_DOCK")
OUT = BASE / "sequential/reconstructed_poses"
DRIVE = Path("/content/drive/MyDrive/Dock_EM_Dock_Research/CYP3A4/sequential/reconstructed_poses")
RECEPTORS = BASE / "sequential/reconstructed_receptors"
LIGANDS = BASE / "sequential/ligand_pdbqt"

CENTERS = {
    "EE": (-24.301, -12.517, -13.995),
    "LNG": (-23.900, -27.553, -34.601),
    "NGMN": (-21.003, -15.678, -11.734),
    "MPA": (-24.652, -27.063, -33.990)
}

names = ["EE", "LNG", "NGMN", "MPA"]
results = []

for a in names:
    receptor = RECEPTORS / f"CYP3A4_{a}_fixedA.pdbqt"

    for b in names:
        if a == b:
            continue

        pair = f"{a}_to_{b}"
        ligand = LIGANDS / f"{b}_B.pdbqt"
        pose = OUT / f"{pair}_docked.pdbqt"

        print(f"\nDocking {a} → {b}", flush=True)

        row = {
            "sequence": f"{a}→{b}",
            "ligand_A": a,
            "ligand_B": b,
            "score_kcal_mol": "",
            "pose_file": str(pose),
            "status": "FAILED"
        }

        try:
            v = Vina(sf_name="vina", cpu=2, seed=42, verbosity=0)
            v.set_receptor(str(receptor))
            v.set_ligand_from_file(str(ligand))
            v.compute_vina_maps(center=list(CENTERS[b]), box_size=[26, 26, 26])
            v.dock(exhaustiveness=16, n_poses=10)
            energies = v.energies(n_poses=1)
            v.write_poses(str(pose), n_poses=1, overwrite=True)

            if not pose.exists() or pose.stat().st_size == 0:
                raise RuntimeError("No docking pose written")

            row["score_kcal_mol"] = float(energies[0][0])
            row["status"] = "SAVED"

            shutil.copy2(pose, DRIVE / pose.name)
            print("Score:", row["score_kcal_mol"], "kcal/mol", flush=True)
            print("Pose backed up to Drive.", flush=True)

        except Exception as exc:
            row["error"] = str(exc)
            print("ERROR:", exc, flush=True)

        results.append(row)

        # Update CSV and back it up after every run.
        csv_path = OUT / "Step179_Docking_Results.csv"
        with csv_path.open("w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "sequence", "ligand_A", "ligand_B",
                    "score_kcal_mol", "pose_file", "status", "error"
                ],
                extrasaction="ignore"
            )
            writer.writeheader()
            writer.writerows(results)

        shutil.copy2(csv_path, DRIVE / csv_path.name)

print("\n=== STEP 179 COMPLETE ===")
print("Successful dockings:", sum(r["status"] == "SAVED" for r in results), "/ 12")
'''

script_path = OUT / "Step179_Run_Docking.py"
script_path.write_text(script)
shutil.copy2(script_path, DRIVE / script_path.name)

env = os.environ.copy()
env["CONDA_PREFIX"] = str(ENV)
env["PATH"] = str(ENV / "bin") + os.pathsep + env.get("PATH", "")
env.pop("PYTHONPATH", None)
env.pop("PYTHONHOME", None)

result = subprocess.run(
    [str(ENV / "bin/python"), "-u", str(script_path)],
    env=env,
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT
)

print(result.stdout)
print("\nExit code:", result.returncode)

```

### Source lines 26679–26791
```python

# Run a single receptor-loading test before launching all 12 jobs.
script = r'''
from vina import Vina
from pathlib import Path
import csv, shutil

BASE = Path("/content/CYP3A4_DOCK_EM_DOCK")
REC = BASE / "sequential/reconstructed_receptors"
LIG = BASE / "sequential/ligand_pdbqt"
POSE = BASE / "sequential/reconstructed_poses"
DRIVE = Path("/content/drive/MyDrive/Dock_EM_Dock_Research/CYP3A4/sequential/reconstructed_poses")

names = ["EE", "LNG", "NGMN", "MPA"]

centers = {
    "EE": [-24.301, -12.517, -13.995],
    "LNG": [-23.900, -27.553, -34.601],
    "NGMN": [-21.003, -15.678, -11.734],
    "MPA": [-24.652, -27.063, -33.990]
}

print("\nTesting receptor compatibility...", flush=True)

test = Vina(sf_name="vina", cpu=2, seed=42, verbosity=0)
test.set_receptor(str(REC / "CYP3A4_EE_fixedA.pdbqt"))
print("Receptor accepted by Vina.", flush=True)

results = []

for a in names:
    for b in names:
        if a == b:
            continue

        sequence = f"{a}→{b}"
        output = POSE / f"{a}_to_{b}_docked.pdbqt"

        print("\nDocking", sequence, flush=True)

        row = {
            "sequence": sequence,
            "score_kcal_mol": "",
            "status": "FAILED",
            "error": ""
        }

        try:
            v = Vina(sf_name="vina", cpu=2, seed=42, verbosity=0)
            v.set_receptor(str(REC / f"CYP3A4_{a}_fixedA.pdbqt"))
            v.set_ligand_from_file(str(LIG / f"{b}_B.pdbqt"))
            v.compute_vina_maps(center=centers[b], box_size=[26,26,26])
            v.dock(exhaustiveness=16, n_poses=10)

            energy = float(v.energies(n_poses=1)[0][0])
            v.write_poses(str(output), n_poses=1, overwrite=True)

            shutil.copy2(output, DRIVE / output.name)

            row["score_kcal_mol"] = energy
            row["status"] = "SAVED"

            print("Score:", energy, "kcal/mol", flush=True)

        except Exception as e:
            row["error"] = str(e)
            print("FAILED:", e, flush=True)

        results.append(row)

        csv_file = POSE / "Step180_Docking_Results.csv"
        with csv_file.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(row.keys()))
            writer.writeheader()
            writer.writerows(results)

        shutil.copy2(csv_file, DRIVE / csv_file.name)

print("\nSUCCESSFUL:", sum(r["status"] == "SAVED" for r in results), "/12")
'''

script_path = POSE / "Step180_Run_Docking.py"
script_path.write_text(script)
shutil.copy2(script_path, DRIVE / "reconstructed_poses" / script_path.name)

environment = os.environ.copy()
environment["CONDA_PREFIX"] = str(ENV)
environment["PATH"] = str(ENV / "bin") + os.pathsep + environment.get("PATH", "")
environment.pop("PYTHONPATH", None)
environment.pop("PYTHONHOME", None)

result = subprocess.run(
    [str(ENV / "bin/python"), "-u", str(script_path)],
    env=environment,
    capture_output=True,
    text=True
)

print(result.stdout)
if result.stderr:
    print(result.stderr[-4000:])
print("Exit code:", result.returncode)

# ============================================================
# STEP 181 — GENERATE CYP3A4 TABLES 6–9
# Corrected heavy-atom handling and docking pose validation
# ============================================================

from pathlib import Path
import numpy as np
import csv
import shutil
import subprocess
```

### Source lines 27213–27303
```python
        })

        # ----------------------------------------------------
        # TABLE 8: LIGAND B CENTROIDS
        # ----------------------------------------------------

        table8.append({
            "Sequence": sequence,
            "Ligand B": b,
            "Heavy atom count": len(ligand_b),
            "Centroid X (Å)":
                round(float(centroid_b[0]), 3),
            "Centroid Y (Å)":
                round(float(centroid_b[1]), 3),
            "Centroid Z (Å)":
                round(float(centroid_b[2]), 3)
        })

        # ----------------------------------------------------
        # TABLE 9: INTERLIGAND DISTANCES
        # ----------------------------------------------------

        table9.append({
            "Sequence": sequence,
            "Interligand centroid distance (Å)":
                round(centroid_distance, 3),
            "Minimum heavy atom distance (Å)":
                round(minimum_distance, 3)
        })

        validation.append({
            "Sequence": sequence,
            "Ligand A heavy atoms": len(ligand_a),
            "Ligand B heavy atoms": len(ligand_b),
            "Minimum interligand distance (Å)":
                round(minimum_distance, 3),
            "Warning": warning
        })

        print(
            f"{sequence:12s} | "
            f"centroid = {centroid_distance:8.3f} Å | "
            f"minimum = {minimum_distance:7.3f} Å | "
            f"shared residues = {len(shared):2d} | "
            f"{warning}",
            flush=True
        )


# ------------------------------------------------------------
# SAVE ALL CSV TABLES
# ------------------------------------------------------------

print("\nSaving CSV files...")

paths = [
    save_csv(
        "Table_6_Protein_Contacts.csv",
        table6
    ),
    save_csv(
        "Table_7_Binding_Configurations.csv",
        table7
    ),
    save_csv(
        "Table_8_Ligand_B_Centroids.csv",
        table8
    ),
    save_csv(
        "Table_9_Interligand_Distances.csv",
        table9
    ),
    save_csv(
        "Step181_Validation.csv",
        validation
    )
]


# ------------------------------------------------------------
# CREATE FORMATTED EXCEL WORKBOOK
# ------------------------------------------------------------

print("Creating Excel workbook...")

wb = Workbook()
wb.remove(wb.active)

sheets = [
    ("Table 6", table6),
    ("Table 7", table7),
```


---

# Glutathione peroxidase (GPx)

# GPx — original Dock–EM–Dock code excerpts
**Source:** `gpx_energyminimization.py` (4,301 lines).
**Provenance:** All code excerpts below are copied verbatim from the uploaded Colab-exported Python source, with original line ranges. Sections are excerpts, not a standalone runnable script. Original complete file is included in the ZIP.
**Scope:** The original notebook may contain debugging, corrections, alternate parameterizations, or non-executable Colab shell commands. These excerpts do not imply all code paths were executed successfully.

## Original code excerpts
### Source lines 71–159
```python
    print(type(e).__name__ + ":", str(e))

# Commented out IPython magic to ensure Python compatibility.

# %pip uninstall -y torch torchvision

# %pip install --no-cache-dir \
    "torch==2.11.0" \
    "torchvision==0.26.0"

import torch
import torchvision
import openmm
from rdkit import Chem
from openff.toolkit import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

print("PyTorch:", torch.__version__)
print("torchvision:", torchvision.__version__)
print("OpenMM:", openmm.__version__)
print("RDKit: OK")
print("OpenFF Toolkit: OK")
print("OpenMM ForceFields: OK")

print("\n=== READY FOR GPX1 ENERGY MINIMIZATION ===")

from pathlib import Path
from rdkit import Chem

files = {
    "EE": {
        "complex": "gpxreceptor_ee_out_2.-6.5.complex.pdb",
        "ligand": "gpxreceptor_ee_out_2.-6.5.pdb"
    },
    "LNG": {
        "complex": "gpxreceptor_lng_out_1.-6.3.complex.pdb",
        "ligand": "gpxreceptor_lng_out_1.-6.3.pdb"
    },
    "NGMN": {
        "complex": "gpxreceptor_ngmn_out_1.-6.1.complex.pdb",
        "ligand": "gpxreceptor_ngmn_out_1.-6.1.pdb"
    },
    "MPA": {
        "complex": "gpxreceptor_mpa_out_1.-6.2.complex.pdb",
        "ligand": "gpxreceptor_mpa_out_1.-6.2.pdb"
    }
}

print("=== GPX1 STRUCTURE VALIDATION ===")

for name, paths in files.items():
    print(f"\n--- GPX1–{name} ---")

    for kind, filename in paths.items():
        path = Path("/content") / filename

        if not path.exists():
            print(f"{kind}: FILE NOT FOUND")
            continue

        lines = path.read_text(errors="replace").splitlines()

        atom_count = sum(
            line.startswith(("ATOM  ", "HETATM"))
            for line in lines
        )

        print(f"{kind}: {atom_count} atom records")

    ligand_path = Path("/content") / paths["ligand"]

    if ligand_path.exists():
        mol = Chem.MolFromPDBFile(
            str(ligand_path),
            removeHs=False,
            sanitize=False
        )

        if mol is None:
            print("Ligand: RDKit could not read structure")
        else:
            print("Ligand RDKit atoms:", mol.GetNumAtoms())
            print("Ligand conformers:", mol.GetNumConformers())
            print(
                "Elements:",
                sorted(set(
                    atom.GetSymbol()
                    for atom in mol.GetAtoms()
                ))
```

### Source lines 959–1102
```python
    print("Hydrogen addition: SUCCESS")
    print("Total atoms:", modeller.topology.getNumAtoms())
    print("Total residues:", modeller.topology.getNumResidues())
    print("Saved:", output_path.name)

except Exception as e:
    print("Hydrogen addition: FAILED")
    print(type(e).__name__, str(e)[:1500])

print("\n=== HYDROGEN PREPARATION COMPLETE ===")

from openmm import app

print("=== GPX1 AMBER14 FORCE-FIELD CHECK ===")

try:
    pdb = app.PDBFile(
        "/content/GPX1_prepared/GPX1_protein_H.pdb"
    )

    forcefield = app.ForceField(
        "amber14/protein.ff14SB.xml"
    )

    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None
    )

    print("Force-field parameterization: SUCCESS")
    print("Particles:", system.getNumParticles())
    print("Forces:", system.getNumForces())
    print("Protein ready for complex preparation!")

except Exception as e:
    print("Force-field parameterization: FAILED")
    print(type(e).__name__, str(e)[:1500])

print("\n=== FORCE-FIELD CHECK COMPLETE ===")

from pathlib import Path
from openff.toolkit import Molecule
from openff.units import unit as off_unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openmm import app

base_dir = Path("/content/GPX1_prepared")
ligand_dir = base_dir / "ligands_hydrogenated"

off_molecules = {}

print("=== GPX1 OPENFF LIGAND VALIDATION ===")

for name in ["EE", "LNG", "NGMN", "MPA"]:

    ligand_path = ligand_dir / f"GPX1_{name}_H.sdf"

    try:
        molecule = Molecule.from_file(
            str(ligand_path),
            file_format="SDF",
            allow_undefined_stereo=False
        )

        generator = SMIRNOFFTemplateGenerator(
            molecules=[molecule],
            forcefield="openff-2.2.1"
        )

        # Test parameterization using a ligand-only topology
        topology = molecule.to_topology()
        omm_topology = topology.to_openmm()

        forcefield = app.ForceField()
        forcefield.registerTemplateGenerator(
            generator.generator
        )

        system = forcefield.createSystem(
            omm_topology,
            nonbondedMethod=app.NoCutoff,
            constraints=None
        )

        off_molecules[name] = molecule

        print(f"\n{name}: SUCCESS")
        print("Atoms:", molecule.n_atoms)
        print("OpenFF parameterization: SUCCESS")

    except Exception as e:
        print(f"\n{name}: FAILED")
        print(type(e).__name__, str(e)[:1000])

print("\n=== LIGAND VALIDATION COMPLETE ===")

from openmm import app, unit
from pathlib import Path
import numpy as np

base_dir = Path("/content/GPX1_prepared")

protein = app.PDBFile(
    str(base_dir / "GPX1_protein_H.pdb")
)

protein_xyz = np.array(
    protein.positions.value_in_unit(unit.angstrom)
)

print("=== GPX1 PROTEIN–LIGAND COORDINATE CHECK ===")

for name, molecule in off_molecules.items():

    ligand_xyz = molecule.conformers[0].m_as(
        off_unit.angstrom
    )

    distances = np.linalg.norm(
        protein_xyz[:, None, :] -
        ligand_xyz[None, :, :],
        axis=2
    )

    minimum_distance = distances.min()

    print(f"\nGPX1–{name}")
    print("Protein atoms:", len(protein_xyz))
    print("Ligand atoms:", len(ligand_xyz))
    print(
        "Minimum protein–ligand distance:",
        round(float(minimum_distance), 3),
        "Å"
    )

print("\n=== COORDINATE CHECK COMPLETE ===")

from rdkit import Chem
from pathlib import Path

base = Path("/content/GPX1_prepared")

paths = {
```

### Source lines 1250–1549
```python
print("C=N stereochemistry:", cn_bonds[0].GetStereo())

if cn_bonds[0].GetStereo() not in (
    Chem.BondStereo.STEREOE,
    Chem.BondStereo.STEREOZ
):
    raise ValueError("C=N stereochemistry remains undefined")

# Save corrected structure.
writer = Chem.SDWriter(str(corrected_path))
writer.write(mol)
writer.close()

print("Saved:", corrected_path.name)

print("\n=== OPENFF PARAMETERIZATION TEST ===")

try:
    molecule = Molecule.from_file(
        str(corrected_path),
        file_format="SDF",
        allow_undefined_stereo=False
    )

    generator = SMIRNOFFTemplateGenerator(
        molecules=[molecule],
        forcefield="openff-2.2.1"
    )

    forcefield = app.ForceField()
    forcefield.registerTemplateGenerator(
        generator.generator
    )

    system = forcefield.createSystem(
        molecule.to_topology().to_openmm(),
        nonbondedMethod=app.NoCutoff,
        constraints=None
    )

    off_molecules["NGMN"] = molecule

    print("NGMN OpenFF parameterization: SUCCESS")
    print("Atoms:", molecule.n_atoms)
    print("Particles:", system.getNumParticles())

except Exception as e:
    print("NGMN OpenFF parameterization: FAILED")
    print(type(e).__name__, str(e)[:1500])

print("\n=== STEP 15C COMPLETE ===")

from pathlib import Path
import numpy as np
from rdkit import Chem
from openmm import app, unit

base_dir = Path("/content/GPX1_prepared")

protein = app.PDBFile(
    str(base_dir / "GPX1_protein_H.pdb")
)

protein_atoms = list(protein.topology.atoms())

protein_xyz = np.array(
    protein.positions.value_in_unit(unit.angstrom)
)

protein_heavy = np.array([
    atom.element.symbol != "H"
    for atom in protein_atoms
])

print("=== GPX1 PROTEIN–LIGAND CLASH ANALYSIS ===")

for name in ["EE", "LNG", "NGMN", "MPA"]:

    filename = (
        f"GPX1_{name}_H_E.sdf"
        if name == "NGMN"
        else f"GPX1_{name}_H.sdf"
    )

    ligand_path = (
        base_dir / "ligands_hydrogenated" / filename
    )

    ligand = Chem.SDMolSupplier(
        str(ligand_path),
        removeHs=False
    )[0]

    if ligand is None:
        print(f"{name}: FAILED to load")
        continue

    ligand_xyz = np.array(
        ligand.GetConformer().GetPositions()
    )

    ligand_heavy = np.array([
        atom.GetAtomicNum() > 1
        for atom in ligand.GetAtoms()
    ])

    all_distances = np.linalg.norm(
        protein_xyz[:, None, :] -
        ligand_xyz[None, :, :],
        axis=2
    )

    heavy_distances = all_distances[
        np.ix_(protein_heavy, ligand_heavy)
    ]

    print(f"\n--- GPX1–{name} ---")

    print(
        "Minimum all-atom distance:",
        round(float(all_distances.min()), 3),
        "Å"
    )

    print(
        "Minimum heavy-atom distance:",
        round(float(heavy_distances.min()), 3),
        "Å"
    )

    print(
        "Heavy-atom pairs below 2.0 Å:",
        int(np.sum(heavy_distances < 2.0))
    )

    print(
        "Heavy-atom pairs below 2.5 Å:",
        int(np.sum(heavy_distances < 2.5))
    )

print("\n=== CLASH ANALYSIS COMPLETE ===")

from pathlib import Path
from openmm import app, unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Molecule

base_dir = Path("/content/GPX1_prepared")
ligand_dir = base_dir / "ligands_hydrogenated"

protein = app.PDBFile(
    str(base_dir / "GPX1_protein_H.pdb")
)

combined_systems = {}

print("=== GPX1 COMPLEX SYSTEM CONSTRUCTION ===")

for name in ["EE", "LNG", "NGMN", "MPA"]:

    print(f"\n--- GPX1–{name} ---")

    try:
        ligand_filename = (
            "GPX1_NGMN_H_E.sdf"
            if name == "NGMN"
            else f"GPX1_{name}_H.sdf"
        )

        ligand = Molecule.from_file(
            str(ligand_dir / ligand_filename),
            file_format="SDF",
            allow_undefined_stereo=False
        )

        ligand_topology = ligand.to_topology().to_openmm()

        ligand_positions = (
            ligand.conformers[0].m_as(
                __import__("openff.units").units.unit.nanometer
            ) * unit.nanometer
        )

        modeller = app.Modeller(
            protein.topology,
            protein.positions
        )

        modeller.add(
            ligand_topology,
            ligand_positions
        )

        forcefield = app.ForceField(
            "amber14/protein.ff14SB.xml"
        )

        generator = SMIRNOFFTemplateGenerator(
            molecules=[ligand],
            forcefield="openff-2.2.1"
        )

        forcefield.registerTemplateGenerator(
            generator.generator
        )

        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=None
        )

        combined_systems[name] = {
            "system": system,
            "topology": modeller.topology,
            "positions": modeller.positions,
            "protein_atoms": protein.topology.getNumAtoms()
        }

        print("System construction: SUCCESS")
        print("Total atoms:", modeller.topology.getNumAtoms())
        print("System particles:", system.getNumParticles())
        print("Forces:", system.getNumForces())

    except Exception as e:
        print("System construction: FAILED")
        print(type(e).__name__, str(e)[:1500])

print("\n=== SYSTEM CONSTRUCTION COMPLETE ===")

import openmm as mm
from openmm import app, unit
from pathlib import Path
import pandas as pd
import numpy as np

base_dir = Path("/content/GPX1_prepared")
min_dir = base_dir / "minimized_complexes"
min_dir.mkdir(exist_ok=True)

results = []

print("=== GPX1 RESTRAINED ENERGY MINIMIZATION ===")

for name, data in combined_systems.items():

    print(f"\n--- GPX1–{name} ---")

    try:
        # Copy the system to preserve the original.
        system = mm.XmlSerializer.deserialize(
            mm.XmlSerializer.serialize(data["system"])
        )

        topology = data["topology"]
        positions = data["positions"]
        n_protein = data["protein_atoms"]

        # Add harmonic positional restraints
        # to protein atoms only.
        restraint = mm.CustomExternalForce(
            "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
        )

        restraint.addGlobalParameter(
            "k", 1000.0 * unit.kilojoule_per_mole /
                 unit.nanometer**2
        )

        for parameter in ["x0", "y0", "z0"]:
            restraint.addPerParticleParameter(parameter)

        for i in range(n_protein):
            pos = positions[i].value_in_unit(
                unit.nanometer
            )
            restraint.addParticle(i, list(pos))

        system.addForce(restraint)

        integrator = mm.VerletIntegrator(
            0.001 * unit.picoseconds
        )

        platform = mm.Platform.getPlatformByName("CPU")

        simulation = app.Simulation(
            topology,
            system,
            integrator,
            platform
        )

        simulation.context.setPositions(positions)

        # Initial potential energy
        initial_state = simulation.context.getState(
            getEnergy=True
        )

```

### Source lines 1550–1764
```python
        initial_energy = (
            initial_state.getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
        )

        print(
            "Initial energy:",
            round(initial_energy, 2),
            "kJ/mol"
        )

        # Energy minimization
        simulation.minimizeEnergy(
            tolerance=10.0 *
                      unit.kilojoule_per_mole /
                      unit.nanometer,
            maxIterations=1000
        )

        # Final potential energy and coordinates
        final_state = simulation.context.getState(
            getEnergy=True,
            getPositions=True
        )

        final_energy = (
            final_state.getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
        )

        delta_energy = final_energy - initial_energy

        print(
            "Final energy:",
            round(final_energy, 2),
            "kJ/mol"
        )

        print(
            "Energy change:",
            round(delta_energy, 2),
            "kJ/mol"
        )

        # Save minimized protein-ligand complex.
        output_pdb = (
            min_dir / f"GPX1_{name}_minimized.pdb"
        )

        with open(output_pdb, "w") as f:
            app.PDBFile.writeFile(
                topology,
                final_state.getPositions(),
                f,
                keepIds=True
            )

        results.append({
            "Ligand": name,
            "Initial Energy (kJ/mol)": initial_energy,
            "Final Energy (kJ/mol)": final_energy,
            "Delta E (kJ/mol)": delta_energy
        })

        print("Saved:", output_pdb.name)

        del simulation
        del integrator

    except Exception as e:
        print("MINIMIZATION FAILED")
        print(type(e).__name__, str(e)[:1500])

print("\n=== MINIMIZATION RESULTS ===")

df = pd.DataFrame(results)

if not df.empty:
    print(df.round(2).to_string(index=False))

    csv_path = (
        base_dir / "GPX1_energy_minimization_results.csv"
    )

    df.to_csv(csv_path, index=False)

    print("\nResults saved:", csv_path)

print("\n=== ENERGY MINIMIZATION COMPLETE ===")

from pathlib import Path
from openmm import app, unit
import numpy as np

min_dir = Path(
    "/content/GPX1_prepared/minimized_complexes"
)

expected_atoms = {
    "EE": 2927,
    "LNG": 2932,
    "NGMN": 2934,
    "MPA": 2943
}

print("=== GPX1 MINIMIZED STRUCTURE VALIDATION ===")

for ligand, expected in expected_atoms.items():

    path = min_dir / f"GPX1_{ligand}_minimized.pdb"

    print(f"\n--- GPX1–{ligand} ---")

    if not path.exists():
        print("FAILED: File not found")
        continue

    try:
        pdb = app.PDBFile(str(path))

        atom_count = pdb.topology.getNumAtoms()

        coordinates = np.array(
            pdb.positions.value_in_unit(unit.angstrom)
        )

        finite = np.isfinite(coordinates).all()

        print("Atoms:", atom_count)
        print("Expected atoms:", expected)
        print("Coordinates finite:", finite)

        if atom_count == expected and finite:
            print("Validation: PASS")
        else:
            print("Validation: NEEDS REVIEW")

    except Exception as e:
        print("Validation: FAILED")
        print(type(e).__name__, str(e)[:500])

print("\n=== VALIDATION COMPLETE ===")

from pathlib import Path
from collections import Counter

min_dir = Path(
    "/content/GPX1_prepared/minimized_complexes"
)

print("=== GPX1 SEQUENTIAL DOCKING PREPARATION ===")

for name in ["EE", "LNG", "NGMN", "MPA"]:

    path = min_dir / f"GPX1_{name}_minimized.pdb"
    lines = path.read_text(errors="replace").splitlines()

    atom_records = [
        line for line in lines
        if line.startswith("ATOM  ")
    ]

    hetatm_records = [
        line for line in lines
        if line.startswith("HETATM")
    ]

    hetatm_residues = Counter(
        (
            line[17:20].strip(),
            line[21:22].strip(),
            line[22:26].strip()
        )
        for line in hetatm_records
    )

    print(f"\n--- GPX1–{name} ---")
    print("ATOM records:", len(atom_records))
    print("HETATM records:", len(hetatm_records))
    print("HETATM residues:", dict(hetatm_residues))

    if hetatm_records:
        print("First HETATM record:")
        print(hetatm_records[0])

print("\n=== INSPECTION COMPLETE ===")

from pathlib import Path
from collections import Counter

base_dir = Path("/content/GPX1_prepared")
min_dir = base_dir / "minimized_complexes"

seq_dir = base_dir / "sequential_docking"
seq_dir.mkdir(exist_ok=True)

ligand_codes = {
    "EE": "EES",
    "LNG": "LNG",
    "NGMN": "NGM",
    "MPA": "MPA"
}

expected_ligand_heavy = {
    "EE": 22,
    "LNG": 23,
    "NGMN": 24,
    "MPA": 28
}

print("=== GPX1 SEQUENTIAL RECEPTOR PREPARATION ===")

for name in ["EE", "LNG", "NGMN", "MPA"]:

    input_path = min_dir / f"GPX1_{name}_minimized.pdb"
```

### Source lines 1863–2162
```python
    centroid = ligand_coords.mean(axis=0)

    residue_names = set(
        line[17:20].strip()
        for line in ligand_lines
    )

    valid = (
        len(protein_lines) == 1446
        and len(ligand_lines) == expected_heavy_atoms[name]
        and len(residue_names) == 1
        and np.isfinite(ligand_coords).all()
    )

    print(f"\n--- GPX1–{name} ---")
    print("Protein heavy atoms:", len(protein_lines))
    print("Ligand heavy atoms:", len(ligand_lines))
    print("Ligand residue:", residue_names)
    print(
        "Ligand centroid (Å):",
        np.round(centroid, 3)
    )
    print("Validation:", "PASS" if valid else "FAILED")

print("\n=== SEQUENTIAL DOCKING COMBINATIONS ===")

combinations = list(permutations(ligands, 2))

for i, (first, second) in enumerate(combinations, 1):
    print(
        f"{i:02d}. GPX1–{first} + {second}"
    )

print("\nTotal combinations:", len(combinations))
print("\n=== STEP 21A COMPLETE ===")

from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd

seq_dir = Path(
    "/content/GPX1_prepared/sequential_docking"
)

ligands = ["EE", "LNG", "NGMN", "MPA"]

ligand_coordinates = {}

for name in ligands:

    path = seq_dir / f"GPX1_{name}_retained_receptor.pdb"

    coordinates = []

    for line in path.read_text().splitlines():

        if line.startswith("HETATM"):
            coordinates.append([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ])

    ligand_coordinates[name] = np.array(coordinates)

results = []

print("=== GPX1 LIGAND SPATIAL COMPARISON ===")

for first, second in combinations(ligands, 2):

    coords1 = ligand_coordinates[first]
    coords2 = ligand_coordinates[second]

    centroid1 = coords1.mean(axis=0)
    centroid2 = coords2.mean(axis=0)

    centroid_distance = np.linalg.norm(
        centroid1 - centroid2
    )

    distances = np.linalg.norm(
        coords1[:, None, :] -
        coords2[None, :, :],
        axis=2
    )

    minimum_distance = distances.min()

    close_contacts = np.sum(distances < 4.0)

    results.append({
        "Ligand Pair": f"{first}–{second}",
        "Centroid Distance (Å)": centroid_distance,
        "Minimum Atom Distance (Å)": minimum_distance,
        "Atom Pairs <4 Å": int(close_contacts)
    })

df_spatial = pd.DataFrame(results)

print(df_spatial.round(3).to_string(index=False))

output_path = (
    seq_dir / "GPX1_ligand_spatial_comparison.csv"
)

df_spatial.to_csv(output_path, index=False)

print("\nSaved:", output_path)
print("\n=== STEP 21B COMPLETE ===")

from pathlib import Path
import numpy as np
from openmm import app, unit

min_dir = Path(
    "/content/GPX1_prepared/minimized_complexes"
)

ligands = ["EE", "LNG", "NGMN", "MPA"]

backbone_coordinates = {}
backbone_identifiers = {}

print("=== GPX1 RECEPTOR ALIGNMENT VERIFICATION ===")

for name in ligands:

    pdb = app.PDBFile(
        str(min_dir / f"GPX1_{name}_minimized.pdb")
    )

    coords = np.array(
        pdb.positions.value_in_unit(unit.angstrom)
    )

    atoms = list(pdb.topology.atoms())

    selected = [
        atom for atom in atoms
        if atom.residue.chain.id == "A"
        and atom.name in ("N", "CA", "C")
        and atom.element.symbol != "H"
    ]

    identifiers = [
        (
            atom.residue.id,
            atom.residue.name,
            atom.name
        )
        for atom in selected
    ]

    selected_coords = np.array([
        coords[atom.index]
        for atom in selected
    ])

    backbone_coordinates[name] = selected_coords
    backbone_identifiers[name] = identifiers

reference = backbone_coordinates["EE"]
reference_ids = backbone_identifiers["EE"]

for name in ligands:

    mobile = backbone_coordinates[name]

    if backbone_identifiers[name] != reference_ids:
        print(f"{name}: Backbone atom mismatch")
        continue

    # RMSD without superposition:
    # checks whether structures share the same
    # coordinate frame.
    direct_rmsd = np.sqrt(
        np.mean(
            np.sum((mobile - reference)**2, axis=1)
        )
    )

    # Kabsch alignment:
    # checks structural similarity after
    # optimal rigid-body superposition.
    ref_centered = reference - reference.mean(axis=0)
    mob_centered = mobile - mobile.mean(axis=0)

    covariance = mob_centered.T @ ref_centered

    U, S, Vt = np.linalg.svd(covariance)

    correction = np.eye(3)
    correction[-1, -1] = np.sign(
        np.linalg.det(U @ Vt)
    )

    rotation = U @ correction @ Vt

    aligned = mob_centered @ rotation

    aligned_rmsd = np.sqrt(
        np.mean(
            np.sum(
                (aligned - ref_centered)**2,
                axis=1
            )
        )
    )

    print(f"\nGPX1–{name} vs GPX1–EE")
    print("Backbone atoms:", len(reference_ids))
    print("Direct RMSD:", round(direct_rmsd, 4), "Å")
    print("Aligned RMSD:", round(aligned_rmsd, 4), "Å")

print("\n=== STEP 22 COMPLETE ===")

from pathlib import Path
from itertools import permutations
import pandas as pd

base_dir = Path("/content/GPX1_prepared")
seq_dir = base_dir / "sequential_docking"

ligand_dir = base_dir / "ligands_hydrogenated"

ligands = ["EE", "LNG", "NGMN", "MPA"]

records = []

for first, second in permutations(ligands, 2):

    receptor_path = (
        seq_dir /
        f"GPX1_{first}_retained_receptor.pdb"
    )

    ligand_filename = (
        "GPX1_NGMN_H_E.sdf"
        if second == "NGMN"
        else f"GPX1_{second}_H.sdf"
    )

    ligand_path = ligand_dir / ligand_filename

    records.append({
        "First Ligand (Retained)": first,
        "Second Ligand (Docked)": second,
        "Receptor File": receptor_path.name,
        "Ligand File": ligand_filename,
        "Receptor Exists": receptor_path.exists(),
        "Ligand Exists": ligand_path.exists()
    })

manifest = pd.DataFrame(records)

print("=== GPX1 SEQUENTIAL DOCKING MANIFEST ===")

print(manifest.to_string(index=False))

manifest_path = (
    seq_dir / "GPX1_sequential_docking_manifest.csv"
)

manifest.to_csv(manifest_path, index=False)

print("\nTotal docking experiments:", len(manifest))

print(
    "All required files exist:",
    bool(
        manifest["Receptor Exists"].all()
        and manifest["Ligand Exists"].all()
    )
)

print("\nSaved:", manifest_path)
print("\n=== STEP 23A COMPLETE ===")

# Commented out IPython magic to ensure Python compatibility.

# Install tools for sequential docking
# %pip -q install vina meeko gemmi

import vina
import meeko
import gemmi

print("AutoDock Vina:", vina.__version__)
print("Meeko:", meeko.__version__)
print("Gemmi:", gemmi.__version__)
print("\n=== STEP 24 COMPLETE ===")

import subprocess
import shutil

print("=== INSTALLING AUTODOCK VINA ===")

conda = shutil.which("mamba") or shutil.which("conda")
```

### Source lines 2163–2192
```python

if conda is None:
    print("Conda/Mamba not found.")
else:
    subprocess.run(
        [
            conda, "install", "-y",
            "-c", "conda-forge",
            "vina", "meeko", "gemmi"
        ],
        check=True
    )

    print("\n=== VINA INSTALLATION CHECK ===")

    vina_exe = shutil.which("vina")

    if vina_exe:
        subprocess.run(
            [vina_exe, "--version"],
            check=False
        )
    else:
        print("Vina executable not found on PATH.")

print("\n=== STEP 24 COMPLETE ===")

import shutil
import subprocess
import importlib.util
```

### Source lines 2879–3178
```python

ligands = ["EE", "LNG", "NGMN", "MPA"]

# Determine a common search box from protein coordinates.
protein_lines = [
    line for line in (pdbqt / "GPX1_protein.pdbqt").read_text().splitlines()
    if line.startswith(("ATOM  ", "HETATM"))
]

coords = np.array([
    [
        float(line[30:38]),
        float(line[38:46]),
        float(line[46:54])
    ]
    for line in protein_lines
])

minimum = coords.min(axis=0)
maximum = coords.max(axis=0)

center = (minimum + maximum) / 2
size = maximum - minimum + 8.0

print("=== GPX1 SEQUENTIAL DOCKING ===")
print("Box center:", np.round(center, 2))
print("Box size:", np.round(size, 2))

records = []

for first, second in permutations(ligands, 2):

    receptor = pdbqt / f"GPX1_{first}_sequential_receptor.pdbqt"
    ligand = pdbqt / f"{second}.pdbqt"

    output = results_dir / f"GPX1_{first}_then_{second}.pdbqt"
    log = results_dir / f"GPX1_{first}_then_{second}.log"

    command = [
        "vina",
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--exhaustiveness", "16",
        "--num_modes", "9",
        "--cpu", "2",
        "--seed", "42",
        "--out", str(output)
    ]

    print(f"\nDocking: GPX1–{first} + {second}")

    result = subprocess.run(
        command,
        capture_output=True,
        text=True
    )

    log.write_text(result.stdout + "\n" + result.stderr)

    if result.returncode == 0 and output.exists():

        scores = []

        for line in output.read_text().splitlines():
            if line.startswith("REMARK VINA RESULT:"):
                scores.append(float(line.split()[3]))

        best = min(scores) if scores else None

        print("SUCCESS | Best score:", best, "kcal/mol")

        records.append({
            "Retained Ligand": first,
            "Docked Ligand": second,
            "Best Score (kcal/mol)": best,
            "Status": "SUCCESS"
        })

    else:
        print("FAILED:", result.stderr[-500:])

        records.append({
            "Retained Ligand": first,
            "Docked Ligand": second,
            "Best Score (kcal/mol)": None,
            "Status": "FAILED"
        })

df = pd.DataFrame(records)

print("\n=== SEQUENTIAL DOCKING RESULTS ===")
print(df.to_string(index=False))

csv = results_dir / "GPX1_sequential_docking_scores.csv"
df.to_csv(csv, index=False)

print("\nSaved:", csv)
print("\n=== STEP 33 COMPLETE ===")

from pathlib import Path
from itertools import permutations
import numpy as np
import pandas as pd

base = Path("/content/GPX1_prepared/sequential_docking")
receptor_dir = base / "pdbqt"
results_dir = base / "docking_results"

ligands = ["EE", "LNG", "NGMN", "MPA"]
heavy_counts = {"EE": 22, "LNG": 23, "NGMN": 24, "MPA": 28}

def atom_coordinates(lines):
    return np.array([
        [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        for line in lines
    ])

def pdbqt_atoms(path, first_model=False):
    lines = path.read_text().splitlines()

    if first_model and "MODEL" in "\n".join(lines):
        start = next(
            i for i, line in enumerate(lines)
            if line.startswith("MODEL")
        )
        end = next(
            i for i in range(start + 1, len(lines))
            if lines[i].startswith("ENDMDL")
        )
        lines = lines[start:end]

    return [
        line for line in lines
        if line.startswith(("ATOM  ", "HETATM"))
        and line.split()[-1].upper() not in ("H", "HD")
    ]

records = []

print("=== GPX1 SEQUENTIAL SPATIAL ANALYSIS ===")

for first, second in permutations(ligands, 2):

    receptor = (
        receptor_dir /
        f"GPX1_{first}_sequential_receptor.pdbqt"
    )
    docked = (
        results_dir /
        f"GPX1_{first}_then_{second}.pdbqt"
    )

    receptor_atoms = pdbqt_atoms(receptor)
    incoming_atoms = pdbqt_atoms(docked, first_model=True)

    # Retained ligand atoms are at the end of the receptor.
    retained_atoms = receptor_atoms[-heavy_counts[first]:]

    retained_xyz = atom_coordinates(retained_atoms)
    incoming_xyz = atom_coordinates(incoming_atoms)

    distances = np.linalg.norm(
        retained_xyz[:, None, :] -
        incoming_xyz[None, :, :],
        axis=2
    )

    centroid_distance = np.linalg.norm(
        retained_xyz.mean(axis=0) -
        incoming_xyz.mean(axis=0)
    )

    minimum_distance = distances.min()
    contacts_4A = int(np.sum(distances < 4.0))
    clashes_2A = int(np.sum(distances < 2.0))

    records.append({
        "Retained": first,
        "Docked": second,
        "Centroid Distance (Å)": centroid_distance,
        "Minimum Distance (Å)": minimum_distance,
        "Contacts <4 Å": contacts_4A,
        "Clashes <2 Å": clashes_2A
    })

df_spatial = pd.DataFrame(records)

print(df_spatial.round(3).to_string(index=False))

output = results_dir / "GPX1_sequential_spatial_analysis.csv"
df_spatial.to_csv(output, index=False)

print("\nSaved:", output)
print("\n=== STEP 34 COMPLETE ===")

from pathlib import Path
from itertools import permutations
import numpy as np
import pandas as pd

base = Path("/content/GPX1_prepared/sequential_docking")
pdbqt = base / "pdbqt"
results = base / "docking_results"

ligands = ["EE", "LNG", "NGMN", "MPA"]
records = []

def get_atoms(path, first_model=False):
    lines = path.read_text().splitlines()

    if first_model:
        start = next(i for i, x in enumerate(lines) if x.startswith("MODEL"))
        end = next(i for i in range(start + 1, len(lines))
                   if lines[i].startswith("ENDMDL"))
        lines = lines[start:end]

    return [
        x for x in lines
        if x.startswith(("ATOM  ", "HETATM"))
        and x.split()[-1].upper() not in ("H", "HD")
    ]

def xyz(lines):
    return np.array([
        [float(x[30:38]), float(x[38:46]), float(x[46:54])]
        for x in lines
    ])

print("=== GPX1 RESIDUE INTERACTION ANALYSIS ===")

for first, second in permutations(ligands, 2):

    receptor = get_atoms(
        pdbqt / f"GPX1_{first}_sequential_receptor.pdbqt"
    )

    docked = get_atoms(
        results / f"GPX1_{first}_then_{second}.pdbqt",
        first_model=True
    )

    # Protein atoms only: retained ligand records are HETATM.
    protein = [
        x for x in receptor
        if x.startswith("ATOM  ")
    ]

    protein_xyz = xyz(protein)
    ligand_xyz = xyz(docked)

    distances = np.linalg.norm(
        protein_xyz[:, None, :] -
        ligand_xyz[None, :, :],
        axis=2
    )

    contacting_indices = np.where(
        np.any(distances < 4.0, axis=1)
    )[0]

    residues = sorted({
        (
            protein[i][21:22].strip(),
            int(protein[i][22:26]),
            protein[i][17:20].strip()
        )
        for i in contacting_indices
    }, key=lambda x: (x[0], x[1]))

    residue_names = [
        f"{res}{num}" for chain, num, res in residues
    ]

    records.append({
        "Retained": first,
        "Docked": second,
        "Protein Residues <4 Å": ", ".join(residue_names),
        "Residue Count": len(residues)
    })

df = pd.DataFrame(records)

print(df.to_string(index=False))

output = results / "GPX1_sequential_residue_interactions.csv"
df.to_csv(output, index=False)

print("\nSaved:", output)
print("\n=== STEP 35 COMPLETE ===")

from pathlib import Path
from itertools import permutations
import numpy as np
import pandas as pd
```

### Source lines 3179–3478
```python

base = Path("/content/GPX1_prepared/sequential_docking")
pdbqt = base / "pdbqt"
results = base / "docking_results"

ligands = ["EE", "LNG", "NGMN", "MPA"]
records = []

def read_atoms(path, first_model=False):
    lines = path.read_text().splitlines()

    if first_model:
        start = next(i for i, line in enumerate(lines)
                     if line.startswith("MODEL"))
        end = next(i for i in range(start + 1, len(lines))
                   if lines[i].startswith("ENDMDL"))
        lines = lines[start:end]

    return [
        line for line in lines
        if line.startswith(("ATOM  ", "HETATM"))
        and line.split()[-1].upper() not in ("H", "HD")
    ]

def coordinates(lines):
    return np.array([
        [float(x[30:38]), float(x[38:46]), float(x[46:54])]
        for x in lines
    ])

def atom_type(line):
    return line.split()[-1].upper()

print("=== GPX1 INTERACTION TYPE ANALYSIS ===")

for first, second in permutations(ligands, 2):

    receptor = read_atoms(
        pdbqt / f"GPX1_{first}_sequential_receptor.pdbqt"
    )

    protein = [
        line for line in receptor
        if line.startswith("ATOM  ")
    ]

    ligand = read_atoms(
        results / f"GPX1_{first}_then_{second}.pdbqt",
        first_model=True
    )

    p_xyz = coordinates(protein)
    l_xyz = coordinates(ligand)

    distances = np.linalg.norm(
        p_xyz[:, None, :] - l_xyz[None, :, :],
        axis=2
    )

    potential_hbonds = set()
    hydrophobic_contacts = set()

    for i, j in zip(*np.where(distances <= 4.0)):

        p_type = atom_type(protein[i])
        l_type = atom_type(ligand[j])

        residue = (
            protein[i][17:20].strip()
            + protein[i][22:26].strip()
        )

        distance = distances[i, j]

        # Potential polar contacts (distance-based screening).
        if (
            p_type in ("N", "NA", "OA", "SA")
            and l_type in ("N", "NA", "OA", "SA")
            and distance <= 3.5
        ):
            potential_hbonds.add(residue)

        # Carbon-carbon proximity as a hydrophobic-contact proxy.
        if (
            p_type in ("C", "A")
            and l_type in ("C", "A")
            and distance <= 4.0
        ):
            hydrophobic_contacts.add(residue)

    records.append({
        "Retained": first,
        "Docked": second,
        "Potential H-bond Residues": ", ".join(
            sorted(potential_hbonds)
        ) or "None",
        "Hydrophobic Contact Residues": ", ".join(
            sorted(hydrophobic_contacts)
        ) or "None",
        "Potential H-bond Residue Count": len(potential_hbonds),
        "Hydrophobic Residue Count": len(hydrophobic_contacts)
    })

df = pd.DataFrame(records)

print(df.to_string(index=False))

output = results / "GPX1_sequential_interaction_types.csv"
df.to_csv(output, index=False)

print("\nSaved:", output)
print("\n=== STEP 36 COMPLETE ===")

from pathlib import Path
import pandas as pd

results = Path(
    "/content/GPX1_prepared/sequential_docking/docking_results"
)

scores = pd.read_csv(
    results / "GPX1_sequential_docking_scores.csv"
).rename(columns={
    "Retained Ligand": "Retained",
    "Docked Ligand": "Docked"
})

spatial = pd.read_csv(
    results / "GPX1_sequential_spatial_analysis.csv"
)

residues = pd.read_csv(
    results / "GPX1_sequential_residue_interactions.csv"
)

interactions = pd.read_csv(
    results / "GPX1_sequential_interaction_types.csv"
)

keys = ["Retained", "Docked"]

combined = (
    scores
    .merge(spatial, on=keys, validate="one_to_one")
    .merge(residues, on=keys, validate="one_to_one")
    .merge(interactions, on=keys, validate="one_to_one")
)

assert len(combined) == 12, "Expected exactly 12 combinations."

combined = combined.sort_values(
    "Best Score (kcal/mol)"
).reset_index(drop=True)

combined.insert(
    0, "Rank",
    range(1, len(combined) + 1)
)

combined["Hormone Proximity"] = combined[
    "Minimum Distance (Å)"
].apply(
    lambda d: "Close contact (<4 Å)"
    if d < 4 else "Spatially separated"
)

csv_path = results / "GPX1_FINAL_sequential_results.csv"
xlsx_path = results / "GPX1_FINAL_sequential_results.xlsx"

combined.to_csv(csv_path, index=False)

with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
    combined.to_excel(
        writer, sheet_name="All Results", index=False
    )
    scores.to_excel(
        writer, sheet_name="Docking Scores", index=False
    )
    spatial.to_excel(
        writer, sheet_name="Spatial Analysis", index=False
    )
    residues.to_excel(
        writer, sheet_name="Protein Residues", index=False
    )
    interactions.to_excel(
        writer, sheet_name="Interaction Types", index=False
    )

print("=== GPX1 FINAL SEQUENTIAL DOCKING SUMMARY ===")
print(combined[
    [
        "Rank",
        "Retained",
        "Docked",
        "Best Score (kcal/mol)",
        "Minimum Distance (Å)",
        "Contacts <4 Å",
        "Hormone Proximity"
    ]
].round(3).to_string(index=False))

print("\nTotal combinations:", len(combined))
print("Successful dockings:", (combined["Status"] == "SUCCESS").sum())
print("Close-contact combinations:",
      (combined["Minimum Distance (Å)"] < 4).sum())

print("\nSaved CSV:", csv_path)
print("Saved Excel:", xlsx_path)
print("\n=== STEP 37 COMPLETE ===")

from pathlib import Path
from itertools import permutations
import numpy as np
import pandas as pd

results = Path(
    "/content/GPX1_prepared/sequential_docking/docking_results"
)

ligands = ["EE", "LNG", "NGMN", "MPA"]
expected_atoms = {"EE": 22, "LNG": 23, "NGMN": 24, "MPA": 28}

records = []

print("=== STEP 38: LIGAND B CENTROID COORDINATES ===")

for first, second in permutations(ligands, 2):
    path = results / f"GPX1_{first}_then_{second}.pdbqt"

    lines = path.read_text().splitlines()

    # Use only the first (best-scoring) docking pose.
    if any(line.startswith("MODEL") for line in lines):
        start = next(i for i, x in enumerate(lines)
                     if x.startswith("MODEL"))
        end = next(i for i in range(start + 1, len(lines))
                   if lines[i].startswith("ENDMDL"))
        lines = lines[start:end]

    heavy_atoms = [
        line for line in lines
        if line.startswith(("ATOM  ", "HETATM"))
        and line.split()[-1].upper() not in ("H", "HD")
    ]

    coordinates = np.array([
        [
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        ]
        for line in heavy_atoms
    ])

    assert len(coordinates) == expected_atoms[second], (
        f"{first} → {second}: unexpected heavy-atom count"
    )

    centroid = coordinates.mean(axis=0)

    records.append({
        "Ligand A": first,
        "Ligand B": second,
        "Centroid X (Å)": round(float(centroid[0]), 4),
        "Centroid Y (Å)": round(float(centroid[1]), 4),
        "Centroid Z (Å)": round(float(centroid[2]), 4),
        "Heavy Atoms": len(coordinates)
    })

df = pd.DataFrame(records)

print(df.to_string(index=False))

output = results / "GPX1_Table8_LigandB_Centroids.csv"
df.to_csv(output, index=False)

print("\nSaved:", output)
print("=== STEP 38 COMPLETE ===")

from pathlib import Path

base = Path("/content/GPX1_prepared")

print("=== STEP 39: ORIGINAL STRUCTURE FILE CHECK ===")

for name in ["EE", "LNG", "NGMN", "MPA"]:
    print(f"\n--- {name} ---")

    matches = [
        p for p in base.rglob("*")
        if p.is_file()
        and p.suffix.lower() in (".pdb", ".sdf")
        and name.lower() in p.name.lower()
        and "minimized" not in p.name.lower()
        and "sequential_receptor" not in p.name.lower()
    ]

    for path in matches:
        print(path.relative_to(base))

```

### Source lines 3479–3517
```python
    if not matches:
        print("No original structure files found.")

print("\n=== STEP 39 COMPLETE ===")

from pathlib import Path

base = Path("/content/GPX1_prepared")

print("=== STEP 40: STRUCTURE CONTENT VERIFICATION ===")

for ligand in ["EE", "LNG", "NGMN", "MPA"]:
    print(f"\n--- {ligand} ---")

    files = {
        "Original pose": base / f"GPX1_{ligand}_pose1.pdb",
        "Minimized complex": (
            base / "minimized_complexes" /
            f"GPX1_{ligand}_minimized.pdb"
        )
    }

    for label, path in files.items():
        if not path.exists():
            print(f"{label}: FILE NOT FOUND")
            continue

        lines = path.read_text().splitlines()

        protein_atoms = [
            line for line in lines
            if line.startswith("ATOM  ")
        ]

        ligand_atoms = [
            line for line in lines
            if line.startswith("HETATM")
        ]

```

### Source lines 3581–3785
```python
            line[12:16].strip()
        )
        atoms[key] = {
            "coord": np.array([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ]),
            "element": line[76:78].strip().upper()
        }
    return atoms

def align_kabsch(mobile, target):
    m_center = mobile.mean(axis=0)
    t_center = target.mean(axis=0)
    u, s, vt = np.linalg.svd(
        (mobile - m_center).T @ (target - t_center)
    )
    d = np.linalg.det(u @ vt)
    correction = np.diag([1, 1, np.sign(d)])
    rotation = u @ correction @ vt
    translation = t_center - m_center @ rotation
    return rotation, translation

def rmsd(a, b):
    return np.sqrt(np.mean(np.sum((a - b)**2, axis=1)))

original = read_pdb_atoms(original_protein)
results = []

for ligand, sdf_path in ligand_files.items():
    minimized_path = (
        base / "minimized_complexes" /
        f"GPX1_{ligand}_minimized.pdb"
    )
    minimized = read_pdb_atoms(minimized_path)

    ca_keys = sorted(
        key for key in original.keys() & minimized.keys()
        if key[3] == "CA"
    )

    assert len(ca_keys) == 186, (
        f"{ligand}: expected 186 matched Cα atoms"
    )

    orig_ca = np.array([original[k]["coord"] for k in ca_keys])
    min_ca = np.array([minimized[k]["coord"] for k in ca_keys])

    rotation, translation = align_kabsch(min_ca, orig_ca)
    aligned_ca = min_ca @ rotation + translation

    whole_rmsd = rmsd(orig_ca, aligned_ca)

    mol = Chem.SDMolSupplier(
        str(sdf_path), removeHs=False
    )[0]
    assert mol is not None, f"Cannot read {sdf_path}"

    original_ligand = np.array([
        list(mol.GetConformer().GetAtomPosition(a.GetIdx()))
        for a in mol.GetAtoms()
        if a.GetAtomicNum() > 1
    ])

    pocket_residues = set()
    for key, atom in original.items():
        if atom["element"] in ("H", "D"):
            continue
        if cdist(
            atom["coord"][None, :],
            original_ligand
        ).min() <= 5.0:
            pocket_residues.add(key[:3])

    pocket_ca_keys = [
        k for k in ca_keys if k[:3] in pocket_residues
    ]

    pocket_ca_rmsd = rmsd(
        np.array([original[k]["coord"] for k in pocket_ca_keys]),
        np.array([
            minimized[k]["coord"] @ rotation + translation
            for k in pocket_ca_keys
        ])
    )

    pocket_heavy_keys = [
        k for k in original.keys() & minimized.keys()
        if k[:3] in pocket_residues
        and original[k]["element"] not in ("H", "D")
    ]

    pocket_heavy_rmsd = rmsd(
        np.array([original[k]["coord"] for k in pocket_heavy_keys]),
        np.array([
            minimized[k]["coord"] @ rotation + translation
            for k in pocket_heavy_keys
        ])
    )

    min_ligand_atoms = [
        np.array([
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        ])
        for line in minimized_path.read_text().splitlines()
        if line.startswith("HETATM")
        and line[76:78].strip().upper() not in ("H", "D")
    ]

    assert len(min_ligand_atoms) == len(original_ligand), (
        f"{ligand}: ligand heavy-atom count mismatch"
    )

    aligned_min_ligand = (
        np.array(min_ligand_atoms) @ rotation + translation
    )

    centroid_shift = np.linalg.norm(
        original_ligand.mean(axis=0)
        - aligned_min_ligand.mean(axis=0)
    )

    results.append({
        "Ligand A": ligand,
        "Whole-receptor Cα RMSD (Å)": round(whole_rmsd, 4),
        "Pocket Cα RMSD (Å)": round(pocket_ca_rmsd, 4),
        "Pocket heavy-atom RMSD (Å)": round(pocket_heavy_rmsd, 4),
        "Ligand A centroid shift (Å)": round(centroid_shift, 4),
        "Pocket residues": len(pocket_residues)
    })

df = pd.DataFrame(results)

print("=== STEP 42: GPX1 TABLE 4 RESULTS ===")
print(df.to_string(index=False))

output = base / "GPX1_Table4_Structural_Displacement.csv"
df.to_csv(output, index=False)

print("\nSaved:", output)
print("=== STEP 42 COMPLETE ===")

from pathlib import Path

base = Path("/content/GPX1_prepared")

keywords = [
    "minimizeEnergy",
    "maxIterations",
    "tolerance",
    "minimization"
]

print("=== STEP 43: MINIMIZATION SETTINGS SEARCH ===")

found = False

for path in base.rglob("*"):
    if not path.is_file():
        continue

    if path.suffix.lower() not in (".py", ".txt", ".json", ".log"):
        continue

    try:
        lines = path.read_text(errors="ignore").splitlines()
    except Exception:
        continue

    matches = [
        (i, line.strip())
        for i, line in enumerate(lines, 1)
        if any(word.lower() in line.lower() for word in keywords)
    ]

    if matches:
        found = True
        print(f"\nFILE: {path}")
        for number, line in matches[:20]:
            print(f"Line {number}: {line[:200]}")

if not found:
    print("No saved minimization settings found in project text files.")
    print("The original Colab code may still contain the settings.")

print("\n=== STEP 43 COMPLETE ===")

import re
from IPython import get_ipython

print("=== STEP 44: ORIGINAL MINIMIZATION CODE SEARCH ===")

ip = get_ipython()
history = list(ip.history_manager.get_range(
    output=False,
    raw=True
))

found = 0

for session, cell_number, code in history:
    if re.search(
```

### Source lines 3869–4001
```python

        coords.append([
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        ])

    return np.array(coords)

def contact_residues(protein_atoms, ligand_coords):
    protein_coords = np.array([a[0] for a in protein_atoms])

    distances = cdist(protein_coords, ligand_coords)
    contacting = np.where((distances < cutoff).any(axis=1))[0]

    return {
        (protein_atoms[i][1],
         protein_atoms[i][2],
         protein_atoms[i][3])
        for i in contacting
    }

rows = []

for first, second in permutations(ligands, 2):

    complex_path = (
        base / "minimized_complexes" /
        f"GPX1_{first}_minimized.pdb"
    )

    docking_path = (
        results_dir /
        f"GPX1_{first}_then_{second}.pdbqt"
    )

    protein = read_heavy_atoms(complex_path, "ATOM")
    ligand_a = read_heavy_atoms(complex_path, "HETATM")
    ligand_b = read_docked_ligand(docking_path)

    coords_a = np.array([atom[0] for atom in ligand_a])

    contacts_a = contact_residues(protein, coords_a)
    contacts_b = contact_residues(protein, ligand_b)

    shared = contacts_a & contacts_b

    shared_names = ", ".join(
        f"{res}{num}"
        for res, chain, num in sorted(
            shared,
            key=lambda x: int(x[2])
        )
    )

    rows.append({
        "Sequence": f"{first} → {second}",
        "Ligand A GPX1 contacts": len(contacts_a),
        "Ligand B GPX1 contacts": len(contacts_b),
        "Shared contact residues": shared_names or "None"
    })

df = pd.DataFrame(rows)

print("=== STEP 45: GPX1 TABLE 6 RESULTS ===")
print(df.to_string(index=False))

output = results_dir / "GPX1_Table6_Contact_Analysis.csv"
df.to_csv(output, index=False)

print("\nSaved:", output)
print("=== STEP 45 COMPLETE ===")

from pathlib import Path

base = Path("/content/GPX1_prepared")
protein_path = base / "GPX1_protein_repaired.pdb"

residues = {}
for line in protein_path.read_text().splitlines():
    if not line.startswith("ATOM  "):
        continue

    chain = line[21].strip() or "-"
    number = int(line[22:26])
    name = line[17:20].strip()
    atom = line[12:16].strip()

    if atom == "CA":
        residues[(chain, number)] = name

print("=== STEP 46: GPX1 CATALYTIC-SITE REFERENCE ===")
print("Protein:", protein_path.name)
print("Number of residues:", len(residues))
print("Residue numbering range:",
      min(k[1] for k in residues),
      "to",
      max(k[1] for k in residues))

print("\nResidues 30–110:")
for (chain, number), name in sorted(residues.items()):
    if 30 <= number <= 110:
        print(f"Chain {chain} | {name}{number}")

print("\n=== STEP 46 COMPLETE ===")

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from itertools import permutations

base = Path("/content/GPX1_prepared")
results_dir = base / "sequential_docking/docking_results"

ligands = ["EE", "LNG", "NGMN", "MPA"]

# Candidate GPX1 catalytic-region residues in the
# prepared structure; verify their identities below.
candidate_numbers = [47, 74, 93, 94, 162]

protein_path = base / "GPX1_protein_repaired.pdb"

protein_atoms = []
for line in protein_path.read_text().splitlines():
    if not line.startswith("ATOM  "):
        continue
    if line[76:78].strip().upper() in ("H", "D"):
        continue

    protein_atoms.append({
        "resnum": int(line[22:26]),
        "resname": line[17:20].strip(),
```

### Source lines 4042–4130
```python
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line.split()[-1].upper() in ("H", "HD"):
            continue
        coords.append([
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        ])
    return np.array(coords)

catalytic_atoms = [
    a for a in protein_atoms
    if a["resnum"] in candidate_numbers
]

reference_coords = np.array([
    a["coord"] for a in catalytic_atoms
])

rows = []

for first, second in permutations(ligands, 2):
    a = ligand_a_coords(first)
    b = ligand_b_coords(first, second)

    distance_a = cdist(a, reference_coords).min()
    distance_b = cdist(b, reference_coords).min()

    rows.append({
        "Sequence": f"{first} → {second}",
        "Ligand A nearest candidate residue (Å)": round(distance_a, 3),
        "Ligand B nearest candidate residue (Å)": round(distance_b, 3)
    })

df = pd.DataFrame(rows)

print("\n=== STEP 47: CATALYTIC-REGION PROXIMITY ===")
print(df.to_string(index=False))
print("\n=== STEP 47 COMPLETE ===")

from pathlib import Path
from itertools import permutations
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

base = Path("/content/GPX1_prepared")
results_dir = base / "sequential_docking/docking_results"

ligands = ["EE", "LNG", "NGMN", "MPA"]
candidate_numbers = [47, 74, 93, 94, 162]

protein_path = base / "GPX1_protein_repaired.pdb"

reference_atoms = []

for line in protein_path.read_text().splitlines():
    if not line.startswith("ATOM  "):
        continue

    resnum = int(line[22:26])
    if resnum not in candidate_numbers:
        continue

    element = line[76:78].strip().upper()
    if element in ("H", "D"):
        continue

    reference_atoms.append({
        "residue": f"{line[17:20].strip()}{resnum}",
        "atom": line[12:16].strip(),
        "xyz": np.array([
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        ])
    })

def get_ligand_a(first):
    path = (
        base / "minimized_complexes" /
        f"GPX1_{first}_minimized.pdb"
    )
    coords = []

    for line in path.read_text().splitlines():
```

### Source lines 4157–4245
```python
            continue

        coords.append([
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        ])

    return np.array(coords)

ref_xyz = np.array([a["xyz"] for a in reference_atoms])

def nearest_reference(coords):
    distances = cdist(coords, ref_xyz)
    ligand_idx, ref_idx = np.unravel_index(
        distances.argmin(), distances.shape
    )
    return (
        reference_atoms[ref_idx]["residue"],
        round(float(distances[ligand_idx, ref_idx]), 3)
    )

rows = []

for first, second in permutations(ligands, 2):
    residue_a, distance_a = nearest_reference(
        get_ligand_a(first)
    )
    residue_b, distance_b = nearest_reference(
        get_ligand_b(first, second)
    )

    rows.append({
        "Sequence": f"{first} → {second}",
        "Ligand A nearest residue": residue_a,
        "Ligand A distance (Å)": distance_a,
        "Ligand B nearest residue": residue_b,
        "Ligand B distance (Å)": distance_b
    })

df = pd.DataFrame(rows)

print("=== STEP 48: NEAREST CANDIDATE CATALYTIC RESIDUES ===")
print(df.to_string(index=False))
print("\n=== STEP 48 COMPLETE ===")

from pathlib import Path
from itertools import permutations
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

base = Path("/content/GPX1_prepared")
results_dir = base / "sequential_docking/docking_results"

ligands = ["EE", "LNG", "NGMN", "MPA"]

protein_path = base / "GPX1_protein_repaired.pdb"

gly47_atoms = []

for line in protein_path.read_text().splitlines():
    if not line.startswith("ATOM  "):
        continue

    if int(line[22:26]) != 47:
        continue

    if line[76:78].strip().upper() in ("H", "D"):
        continue

    gly47_atoms.append([
        float(line[30:38]),
        float(line[38:46]),
        float(line[46:54])
    ])

gly47_coords = np.array(gly47_atoms)

assert len(gly47_coords) > 0, "GLY47 not found."

def get_ligand_a(first):
    path = (
        base / "minimized_complexes" /
        f"GPX1_{first}_minimized.pdb"
    )
    coords = []

    for line in path.read_text().splitlines():
```

### Source lines 4260–4301
```python
    path = (
        results_dir /
        f"GPX1_{first}_then_{second}.pdbqt"
    )
    coords = []

    for line in path.read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line.split()[-1].upper() in ("H", "HD"):
            continue

        coords.append([
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        ])

    return np.array(coords)

rows = []

for first, second in permutations(ligands, 2):
    a = get_ligand_a(first)
    b = get_ligand_b(first, second)

    distance_a = cdist(a, gly47_coords).min()
    distance_b = cdist(b, gly47_coords).min()

    rows.append({
        "Sequence": f"{first} → {second}",
        "Ligand A to GLY47 (Å)": round(float(distance_a), 3),
        "Ligand B to GLY47 (Å)": round(float(distance_b), 3)
    })

df = pd.DataFrame(rows)

print("=== STEP 49: DISTANCE TO GPX1 GLY47 ===")
print(df.to_string(index=False))
print("\n=== STEP 49 COMPLETE ===")
```
