---
name: mamba-environment
description: Create and configure a new Mamba/Conda environment with a requested program or set of packages, register it as a Jupyter kernel, and verify the install. Use whenever the user asks to create, build, or install a new Mamba/Conda environment, regardless of the software involved.
---

# Mamba environment creation

A general-purpose workflow for building a new Mamba environment. It is not tied to any
specific software, package, repository, or environment name — always take the
`<environment_name>`, `<requested_program>`, and `<requested_packages>` from the user's
request. Never hard-code an example name or package from this skill into a real run.

## 1. Always use this Mamba/Conda installation

Before creating an environment, initialize the correct Conda/Mamba installation in the
terminal:

```bash
source myconda_v1
```

All environments must be created under:

```
/data/kanferg/conda_v1/envs
```

Never create the environment with another Conda/Mamba installation, and never install into
the base environment.

## 2. Create the environment

Use the exact `<environment_name>` requested by the user. Pick a Python version compatible
with `<requested_program>`/`<requested_packages>`; if the user specifies a version, use it
unless it is incompatible with the requested software.

```bash
source myconda_v1
mamba create -y -n <environment_name> python=<compatible_python_version>
mamba activate <environment_name>
```

Confirm the environment lives at `/data/kanferg/conda_v1/envs/<environment_name>` (e.g. with
`mamba env list` or `which python`) before installing anything else.

## 3. Install the requested software

Install whatever program, package, GitHub repository, or set of packages the user asked for.
Choose the installation method appropriate for that software:

- `mamba install <requested_package>` or `mamba install -c conda-forge <requested_package>` /
  `-c bioconda` when a specific channel is required.
- `pip install <requested_package>` when pip is required or recommended (e.g. the package
  isn't on conda-forge/bioconda, or it's a GitHub repo/editable install).
- Follow any install instructions specific to `<requested_program>` (e.g. a documented
  channel order, extra system dependencies, or a `pip install git+...` form).

Prefer Mamba/conda-forge first, and fall back to pip only when necessary.

## 4. Always install Jupyter and ipykernel

Every environment created with this skill must include Jupyter, regardless of what else was
requested:

```bash
mamba install -y -c conda-forge jupyter ipykernel
```

## 5. Always register a Jupyter kernel

Register the environment as a kernel, using the environment name as the kernel name unless
the user asks otherwise:

```bash
python -m ipykernel install --user --name=<environment_name>
```

## 6. Verify the environment

Before reporting success, verify all of the following:

- The correct environment is active.
- It is located under `/data/kanferg/conda_v1/envs/<environment_name>`.
- Python runs: `which python` and `python --version`.
- `<requested_program>`/`<requested_packages>` are installed — check with the relevant
  `--version` flag, or `python -c "import <pkg>; print(<pkg>.__version__)"` when it's a
  Python package.
- Jupyter is installed: `jupyter --version`.
- The kernel was registered: `jupyter kernelspec list` should include `<environment_name>`.

If any check fails, diagnose and fix it (e.g. reinstall the failing package) before reporting
completion — do not report success on an unverified environment.

## 7. Final response

Report a short summary, then the exact commands used:

```
Environment: <environment_name>
Location: /data/kanferg/conda_v1/envs/<environment_name>
Python: <python_version>

Installed:
- <requested_program_or_package>
- <additional_requested_packages>
- Jupyter
- ipykernel

Jupyter kernel:
<environment_name>
```

```bash
source myconda_v1

mamba create -n <environment_name> python=<python_version>
mamba activate <environment_name>

<commands used to install the requested software>

mamba install -c conda-forge jupyter ipykernel

python -m ipykernel install --user --name=<environment_name>
```

Also provide the commands to use the environment later:

```bash
source myconda_v1
mamba activate <environment_name>

jupyter notebook --ip localhost --port $PORT1 --no-browser
```

## Workflow summary

1. User asks for a new environment.
2. Get `<environment_name>` and `<requested_program>`/`<requested_packages>` from the request.
3. `source myconda_v1`.
4. Create the environment under `/data/kanferg/conda_v1/envs`.
5. Install the requested software.
6. Install Jupyter + ipykernel.
7. Register the Jupyter kernel.
8. Verify every step in section 6.
9. Return the installation summary and exact commands.
