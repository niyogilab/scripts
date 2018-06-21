with import <nixpkgs> {};
let
  tecan-merge-metadata =

{ stdenv, rPackages, makeWrapper }:

stdenv.mkDerivation rec {
  name = "tecan-merge-metadata-1.0";
  src = ./.;
  buildInputs = with rPackages; [ R dplyr docopt makeWrapper ];
  installPhase = ''
    mkdir -p $out/bin
    export usageTxt=$src/usage.txt
    substituteAll $src/tecan-merge-metadata.R $out/bin/tecan-merge-metadata
    chmod +x $out/bin/tecan-merge-metadata
    wrapProgram $out/bin/tecan-merge-metadata --set R_LIBS_SITE "$R_LIBS_SITE"
  '';
}

; in callPackage tecan-merge-metadata {}
