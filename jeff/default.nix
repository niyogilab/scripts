with import <nixpkgs> {};

buildEnv {
  name = "knlab-jeff-env";
  paths = [
    # (import ./absorbance-scans)
    # (import ./led-array-programs)
    # (import ./plate-metadata-generator)
    # (import ./tecan-growth-curves)

    (import ./pcc7942-cpf1-primers)
    (import ./plate-map-merge)
    (import ./plate-map-table)
    (import ./qrcode-labels)
    (import ./tecan-extract-table)
    (import ./tecan-merge-metadata)
    (import ./tnseq-growth-curves)
  ];
}
