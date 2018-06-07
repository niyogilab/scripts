{ pythonPackages }:
with pythonPackages;

buildPythonPackage {
  name = "qrlabels-0.1";
  namePrefix = "";
  src = ./.;
  propagatedBuildInputs = [
    docopt
    reportlab
    shortuuid
  ];
  doCheck = false;
  dontStrip = true;
}
