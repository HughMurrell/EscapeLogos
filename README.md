# EscapeLogos

This repository provides a `Julia` script to generate an escape logo 
from a curated amino acid alignment.

### Installing Julia (version 1.10.5)

We recommend you use the `juliaup` version manager to install julia.
from a terminal you can do this as follows:

```bash
curl -fsSL https://install.julialang.org | sh
```

This should install the Julia version manager, `juliaup` as well as
the latest version of Julia. To find out how to use the version manager 
to makesure you have version 1.10.5 as your default, go here:

[https://github.com/JuliaLang/juliaup]

Once Julia is installed, make sure you can enter the julia REPL from 
the command line and check the version number by logging out and in again 

```bash
exit
ssh root@.....
```

and then from your new terminal session:

```bash
juliaup status
julia --version
```

If the version number is not 1.10.5 then you need to use `juliaup` to install
that version and make it the default. 

```bash
juliaup add 1.10.5
juliaup default 1.10.5
```

for further details concerning `juliaup` go here:

[https://github.com/JuliaLang/juliaup?tab=readme-ov-file#using-juliaup]

### cloning the EscapeLogos repository

Now that `Julia` is setup we clone the `EscapeLogos` repository

```bash
git clone https://github.com/HughMurrell/EscapeLogos.git
```

### setting up the Julia package environment

then navigate to the `EscapeLogos` directory and start the 
Julia REPL. 

Enter the package manager using `]` and then enter

```julia
activate .
instantiate
precompile
```

This will activate, install, and precompile the `julia` environment
specified by the  `Project.toml` and `Manifest.toml` files. 
The `precompile` command above is not strictly needed but is useful 
if there are issues with installing the `julia` packages listed in
`Project.toml`

### running the `escapelogos` script

Place your alignments in the alignments directory and run the script

```julia
julia escapelogos.jl
```
and the resulting `svg` logos will appear in the `escape_logos` directory.
