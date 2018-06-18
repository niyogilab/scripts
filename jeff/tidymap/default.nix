with import <nixpkgs> {};
let tidymap =

{ stdenv, makeWrapper, python27Packages }:

stdenv.mkDerivation {
  name = "tidymap-0.1";
  src = ./.;
  buildInputs = [ makeWrapper python27Packages.docopt ];
  installPhase = ''
    mkdir -p $out/bin
    export usageTxt="$(cat $src/usage.txt)"
    substituteAll $src/tidymap.py $out/bin/tidymap
    chmod +x $out/bin/tidymap
    wrapProgram $out/bin/tidymap \
      --prefix PYTHONPATH : $(toPythonPath ${python27Packages.docopt})
  '';
}

; in callPackage tidymap {}
