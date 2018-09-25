{ pythonPackages }:
with pythonPackages;

buildPythonPackage {
  name = "cyano-cpf1-primers-0.1";
  namePrefix = "";
  src = ./.;
  propagatedBuildInputs = [
    docopt
    biopython
    (import ./primer3.nix)
  ];
  doCheck = false;
  dontStrip = true;
}
