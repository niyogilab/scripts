with import <nixpkgs> {};
let
  mergetecan =

{ stdenv, rPackages, makeWrapper }:

stdenv.mkDerivation rec {
  name = "mergetecan-0.1";
  src = ./.;
  buildInputs = with rPackages; [ R dplyr docopt makeWrapper ];
  installPhase = ''
    mkdir -p $out/bin
    export usageTxt=$src/usage.txt
    substituteAll $src/mergetecan.R $out/bin/mergetecan
    chmod +x $out/bin/mergetecan
    wrapProgram $out/bin/mergetecan --set R_LIBS_SITE "$R_LIBS_SITE"
  '';
}

; in callPackage mergetecan {}
