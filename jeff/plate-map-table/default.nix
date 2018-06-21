with import <nixpkgs> {};
let plate-map-table =

{ stdenv, makeWrapper, python27Packages }:

stdenv.mkDerivation {
  name = "plate-map-table-1.0";
  src = ./.;
  buildInputs = [ makeWrapper python27Packages.docopt ];
  installPhase = ''
    mkdir -p $out/bin
    export usageTxt="$(cat $src/usage.txt)"
    substituteAll $src/plate-map-table.py $out/bin/plate-map-table
    chmod +x $out/bin/plate-map-table
    wrapProgram $out/bin/plate-map-table \
      --prefix PYTHONPATH : $(toPythonPath ${python27Packages.docopt})
  '';
}

; in callPackage plate-map-table {}
