I package everything with [Nix][1] and try to keep it clean.
You should be able to build and run each script with all its dependencies like this:

    cd jeff/scriptname
    nix-build
    ./result/bin/scriptname --help

Each one should have an `examples` directory with valid inputs and/or outputs you can try.

Finally there's a top-level [default.nix]() to build everything at once,
which is probably only useful to me.

You can also run things without Nix, but you'll still need to read
`default.nix` to figure out what's needed. Usually some Debian, R, and/or
Python packages.

[1]: https://nixos.org/nix
