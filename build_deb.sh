#!/bin/sh

version=$1
appdir=krait_${version}

#make app folder
mkdir -p ${appdir}/usr/lib
mkdir -p ${appdir}/usr/share/icons
mkdir -p ${appdir}/usr/share/applications
mkdir -p ${appdir}/DEBIAN

#move app to folder
mv dist/Krait ${appdir}/usr/lib

#copy icons
cp src/icons/krait_logo.png ${appdir}/usr/share/icons

#write desktop file
desktop="[Desktop Entry]
Version=${version}
Name=Krait
Comment=microsatellite investigation and primer design
Exec=/usr/lib/Krait/Krait
Icon=/usr/share/icons/krait_logo.png
Terminal=false
Type=Application
Categories=Application
"
echo "$desktop" > ${appdir}/usr/share/applications/krait.desktop
chmod a+x ${appdir}/usr/share/applications/krait.desktop

#write control file
control="package: Krait
version: ${version}
architecture: amd64
maintainer: lmdu <adullb@qq.com>
description: krait installer on ubuntu and debian
"
echo "$control" > ${appdir}/DEBIAN/control

dpkg -b ${appdir} Krait-v${version}-amd64.deb