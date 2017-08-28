# -*- mode: python -*-

block_cipher = None


a = Analysis(['src/main.py'],
             pathex=[''],
             binaries=[],
             datas=[('src/primer3_config', 'perimer3_config'), ('src/template', 'template'), ('src/logo.icns', '.'), ('example', 'example')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=['gtk', 'PyQt4', 'PyQt5', 'Tkinter', 'wx'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='Krait',
          debug=False,
          strip=False,
          upx=True,
          console=False , icon='src/logo.icns')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='Krait')
app = BUNDLE(coll,
             name='Krait.app',
             icon='src/logo.icns',
             bundle_identifier=None)
