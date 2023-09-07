#!/bin/sh

version=$1

wget https://github.com/goreleaser/nfpm/releases/download/v1.7.0/nfpm_1.7.0_Linux_x86_64.tar.gz
tar xzvf nfpm_1.7.0_Linux_x86_64.tar.gz

cat > krait.desktop <<EOF
[Desktop Entry]
Version=${version}
Name=Krait
Comment=microsatellite investigation and primer design
Exec=/usr/lib/Krait/Krait
Icon=/usr/share/icons/krait_logo.png
Terminal=false
Type=Application
Categories=Application
EOF

cat > nfpm.yaml <<EOF
name: Krait
arch: amd64
platform: linux
version: ${version}
maintainer: lmdu <adullb@qq.com>
description: Microsatellite investigation and primer design
vendor: Bioinformatics and Integrative Genomics
homepage: https://github.com/lmdu/krait
license: AGPLv3
files:
  ./Krait/**/*: /usr/lib/Krait
  ./krait.desktop: /usr/share/applications/krait.desktop
  ./krait_logo.png: /usr/share/icons/krait_logo.png
EOF

# copy logo file
cp ../src/icons/krait_logo.png .

./nfpm pkg -t Krait-v${version}-amd64.deb
