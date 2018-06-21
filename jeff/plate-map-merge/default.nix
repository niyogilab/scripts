with import <nixpkgs> {};
let
  plate-map-table = import ../plate-map-table;
  plate-map-merge =

{ stdenv, rPackages, makeWrapper }:

stdenv.mkDerivation rec {
  name = "plate-map-merge-1.0";
  src = ./.;
  buildInputs = with rPackages; [
    R
    docopt
    dplyr
    makeWrapper
    plate-map-table
    reshape
  ];
  installPhase = ''
    mkdir -p $out/bin
    export usageTxt=$src/usage.txt
    substituteAll $src/plate-map-merge.R $out/bin/plate-map-merge
    chmod +x $out/bin/plate-map-merge
    wrapProgram $out/bin/plate-map-merge \
      --set R_LIBS_SITE "$R_LIBS_SITE" \
      --prefix PATH : "${plate-map-table}/bin"
  '';
}

; in callPackage plate-map-merge {}
