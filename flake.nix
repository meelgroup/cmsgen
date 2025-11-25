{
  description = "A heuristic CNF sampler built on CryptoMiniSat";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
  };
  outputs =
    {
      self,
      nixpkgs,
    }:
    let
      inherit (nixpkgs) lib;
      systems = lib.intersectLists lib.systems.flakeExposed lib.platforms.linux;
      forAllSystems = lib.genAttrs systems;
      nixpkgsFor = forAllSystems (system: nixpkgs.legacyPackages.${system});
      fs = lib.fileset;

      cmsgen-package =
        {
          stdenv,
          fetchFromGitHub,
          cmake,
          pkg-config,
          zlib,
        }:
        stdenv.mkDerivation {
          name = "cmsgen";
          src = fs.toSource {
            root = ./.;
            fileset = fs.unions [
              ./src
              ./CMakeLists.txt
              ./cmake
              ./scripts
              ./cmsgenConfig.cmake.in
            ];
          };

          nativeBuildInputs = [
            cmake
            pkg-config
          ];
          buildInputs = [
            zlib
          ];
        };
    in
    {
      packages = forAllSystems (
        system:
        let
          cmsgen = nixpkgsFor.${system}.callPackage cmsgen-package { };
        in
        {
          inherit cmsgen;
          default = cmsgen;
        }
      );
    };
}
