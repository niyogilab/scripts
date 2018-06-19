with import <nixpkgs> {};
let
  tnseq-growth-curves =

{ stdenv, rPackages, makeWrapper }:

stdenv.mkDerivation rec {
  name = "tnseq-growth-curves-1.0";
  src = ./.;
  buildInputs = with rPackages; [
    R
    cowplot
    docopt
    dplyr
    ggplot2
    gridBase # TODO gridExtra?
    makeWrapper
    tidyr
  ];
  installPhase = ''
    mkdir -p $out/bin
    export usageTxt=$src/usage.txt
    substituteAll $src/tnseq-growth-curves.R $out/bin/tnseq-growth-curves
    chmod +x $out/bin/tnseq-growth-curves
    wrapProgram $out/bin/tnseq-growth-curves --set R_LIBS_SITE "$R_LIBS_SITE"
  '';
}

; in callPackage tnseq-growth-curves {}
