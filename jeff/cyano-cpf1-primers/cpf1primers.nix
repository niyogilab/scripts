{ pythonPackages }:
with pythonPackages;

buildPythonPackage {
  name = "cyano-cpf1-primers-0.1";
  namePrefix = "";
  src = ./.;
  propagatedBuildInputs = [
    docopt
  ];
  doCheck = false;
  dontStrip = true;
}
