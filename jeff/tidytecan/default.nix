with import <nixpkgs> {};
let tidytecan =

{ python27Packages, xlsx2csv }:
with python27Packages;

buildPythonPackage {
  name = "tidytecan-0.4";
  namePrefix = "";
  src = ./.;
  propagatedBuildInputs = [ xlsx2csv docopt ];
  doCheck = false;
  dontStrip = true;
}

; in callPackage tidytecan {}
