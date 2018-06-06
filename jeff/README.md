I package everything with [Nix][1] and try to keep it clean.
You should be able to build each project like this:

    cd projectname
    nix-build
    ./result/bin/projectname --help

If not, email and/or file an issue! I'll try to fix it.

Each one should also have an `examples` directory with valid inputs you can try out.

Finally there's a top-level [default.nix]() to build everything at once,
which is probably only useful to me.

[1]: https://nixos.org/nix
