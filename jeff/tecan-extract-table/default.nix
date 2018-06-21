with import <nixpkgs> {};
let tecan-extract-table =

{ python27Packages, xlsx2csv }:
with python27Packages;

buildPythonPackage {
  name = "tecan-extract-table-1.0";
  namePrefix = "";
  src = ./.;
  propagatedBuildInputs = [ xlsx2csv docopt ];
  doCheck = false;
  dontStrip = true;
}

; in callPackage tecan-extract-table {}
