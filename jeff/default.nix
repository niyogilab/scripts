with import <nixpkgs> {};

buildEnv {
  name = "knlab-jeff-env";
  paths = [
    # (import ./absorbance-scans)
    # (import ./led-array-programs)
    # (import ./plate-metadata-generator)
    # (import ./simple-growth-curves)
    # (import ./tecan-growth-curves)
    (import ./pcc7942-cpf1-primers)
    (import ./qrcode-labels)

    # tecan scripts
    (import ./tidymap)
    (import ./tidytecan)
    (import ./mergemeta)
    (import ./mergetecan)
  ];
}
