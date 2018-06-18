with import <nixpkgs> {};
let
  tidymap = import ../tidymap;
  mergemeta =

{ stdenv, rPackages, makeWrapper }:

stdenv.mkDerivation rec {
  name = "mergemeta-0.1";
  src = ./.;
  buildInputs = with rPackages; [ R dplyr docopt reshape makeWrapper tidymap ];
  installPhase = ''
    mkdir -p $out/bin
    export usageTxt=$src/usage.txt
    substituteAll $src/mergemeta.R $out/bin/mergemeta
    chmod +x $out/bin/mergemeta
    wrapProgram $out/bin/mergemeta \
      --set R_LIBS_SITE "$R_LIBS_SITE" \
      --prefix PATH : "${tidymap}/bin"
  '';
}

; in callPackage mergemeta {}
