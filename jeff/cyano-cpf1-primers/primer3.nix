with import <nixpkgs> {};
with pythonPackages;

buildPythonPackage rec {
  pname = "primer3-py";
  version = "0.5.7";
  src = fetchPypi {
    inherit pname version;
    sha256 = "4d128e48124349a107690905aa60cd50739982d9a9ae7e82b0cb3979b40c1fdd";
  };
  propagatedBuildInputs = [ libusb1 cython ];
  doCheck = false; # TODO do check!
}
