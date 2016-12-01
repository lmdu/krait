#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
from db import SSRDB
from ssr import SSRSearchEngine
from PySide.QtCore import *
from PySide.QtGui import *
from PySide.QtSql import *

#import ctypes
#myappid = 'Mencent.GWSSR.SSR.0.1' # arbitrary string
#ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

class MainWindow(QMainWindow):
	def __init__(self):
		super(MainWindow, self).__init__()

		self.setWindowTitle("Gmia v0.1")
		self.setWindowIcon(QIcon(QPixmap("logo.png")))

		self.tabs = QTabWidget()
		self.setCentralWidget(self.tabs)

		self.generalTab = GeneralTab()
		self.tabs.insertTab(0, self.generalTab, "General")
		#self.tabs.insertTab(9, QWidget(), "Statistics")
		#self.tabs.insertTab(1, QWidget(), "Search SSRs")

		self.createActions()
		self.createMenus()
		self.createToolBars()
		self.createStatusBar()

		self.readSettings()

		self.createSSRDb()

		self.show()

	def readSettings(self):
		self.settings = QSettings("config.ini", QSettings.IniFormat)
		self.resize(self.settings.value("size", QSize(900, 600)))

	def writeSettings(self):
		self.settings.setValue("size", self.size())

	def closeEvent(self, event):
		self.writeSettings()

	def createSSRDb(self):
		db = QSqlDatabase.addDatabase('QSQLITE')
		db.setDatabaseName("file:ssr?mode=memory&cache=shared")
		db.open()
		self.db = SSRDB("file:ssr?mode=memory&cache=shared")

	def createActions(self):
		#open a project action
		self.openProjectAct = QAction(self.tr("Open project"), self)
		self.openProjectAct.setShortcut(QKeySequence.Open)
		self.openProjectAct.triggered.connect(self.doOpenProject)

		#close a project action
		self.closeProjectAct = QAction(self.tr("Close project"), self)
		self.closeProjectAct.setShortcut(QKeySequence.Close)
		
		#save a project action
		self.saveProjectAct = QAction(self.tr("Save project"), self)
		self.saveProjectAct.setShortcut(QKeySequence.Save)
		
		#save as a project action
		self.saveAsProjectAct = QAction(self.tr("Save project as..."), self)
		self.saveAsProjectAct.setShortcut(QKeySequence.SaveAs)
		
		#load fasta file or genome action
		self.loadFastaAct = QAction(self.tr("Load fasta sequence"), self)
		self.loadFastaAct.triggered.connect(self.doLoadFasta)
		self.loadFastasAct = QAction(self.tr("Load Fastas in folder"), self)
		self.loadFastasAct.triggered.connect(self.doLoadFastas)
		
		#export the Results
		self.exportResAct = QAction(self.tr("Export Results"), self)
		
		#exit action
		self.exitAct = QAction(self.tr("Exit"), self)
		self.exitAct.setShortcut(QKeySequence.Quit)
		self.exitAct.triggered.connect(self.close)
		
		#copy action
		self.copyAct = QAction(self.tr("Copy"), self)
		self.copyAct.setShortcut(QKeySequence.Copy)
		self.copyAct.triggered.connect(self.doCopy)
		
		self.cutAct = QAction(self.tr("Cut"), self)
		self.cutAct.setShortcut(QKeySequence.Cut)
		self.cutAct.triggered.connect(self.doCut)
		
		self.pasteAct = QAction(self.tr("Paste"), self)
		self.pasteAct.setShortcut(QKeySequence.Paste)
		self.pasteAct.triggered.connect(self.doPaste)
	
		self.selectAllAct = QAction(self.tr("Select All"), self)
		self.selectAllAct.setShortcut(QKeySequence.SelectAll)
		self.selectAllAct.triggered.connect(self.doSelectAll)

		self.preferenceAct = QAction(self.tr("Preferences"), self)
		self.preferenceAct.setShortcut(QKeySequence.Preferences)
		self.preferenceAct.triggered.connect(self.doPreference)

		#toolbar actions
		#search ssrs tool button
		self.searchAct = QAction(QIcon("icons/tandem.png"), self.tr("Search SSRs"), self)
		self.searchAct.triggered.connect(self.doSearchSSRs)
		self.searchRunAct = QAction(QIcon("icons/tandem.png"), self.tr("Run search SSRs"), self)
		self.searchRunAct.triggered.connect(self.doSearchSSRs)
		self.searchParaAct = QAction(self.tr("Set minimum repeats"), self)
		self.searchParaAct.triggered.connect(self.doPreference)

		self.identifyAct = QAction(QIcon("icons/compound.png"), self.tr("Identify cSSRs"), self)
		self.identifyRunAct = QAction(QIcon("icons/compound.png"), self.tr("Run identify cSSRs"), self)
		self.identifyParaAct = QAction(self.tr("dMax setting"), self)
		self.removeIdentifyAct = QAction(self.tr("Remove cSSR results"), self)

	def createMenus(self):
		self.fileMenu = self.menuBar().addMenu("&File")
		self.editMenu = self.menuBar().addMenu("&Edit")
		self.helpMenu = self.menuBar().addMenu("&Help")
		
		self.fileMenu.addAction(self.openProjectAct)
		self.fileMenu.addAction(self.closeProjectAct)
		self.fileMenu.addSeparator()
		self.fileMenu.addAction(self.saveProjectAct)
		self.fileMenu.addAction(self.saveAsProjectAct)
		self.fileMenu.addSeparator()
		self.fileMenu.addAction(self.loadFastaAct)
		self.fileMenu.addAction(self.loadFastasAct)
		self.fileMenu.addSeparator()
		self.fileMenu.addAction(self.exportResAct)
		self.fileMenu.addSeparator()
		self.fileMenu.addAction(self.exitAct)
		
		self.editMenu.addAction(self.copyAct)
		self.editMenu.addAction(self.cutAct)
		self.editMenu.addAction(self.pasteAct)
		self.editMenu.addSeparator()
		self.editMenu.addAction(self.selectAllAct)
		self.editMenu.addSeparator()
		self.editMenu.addAction(self.preferenceAct)

		#tool bar menus
		#search ssrs tool button menu
		self.searchMenu = QMenu()
		self.searchMenu.addAction(self.searchRunAct)
		self.searchMenu.addSeparator()
		self.searchMenu.addAction(self.loadFastaAct)
		self.searchMenu.addAction(self.loadFastasAct)
		self.searchMenu.addSeparator()
		self.searchMenu.addAction(self.searchParaAct)

		self.identifyMenu = QMenu()
		self.identifyMenu.addAction(self.identifyRunAct)
		self.identifyMenu.addAction(self.removeIdentifyAct)
		self.identifyMenu.addSeparator()
		self.identifyMenu.addAction(self.identifyParaAct)
		

	def createToolBars(self):
		self.toolBar = self.addToolBar('')
		self.toolBar.setMovable(False)
		self.toolBar.setIconSize(QSize(36, 36))
		self.toolBar.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)

		#search ssr action and menus
		self.searchAct.setMenu(self.searchMenu)
		self.toolBar.addAction(self.searchAct)

		self.identifyAct.setMenu(self.identifyMenu)
		self.toolBar.addAction(self.identifyAct)

		self.annotToolBtn = QAction(QIcon("icons/annotation.png"), self.tr("Locate SSRs"), self)
		#self.annotToolBtn.setDisabled(True)
		self.annotToolBtnMenu = QMenu()
		self.annotToolBtnMenu.addAction("Settings")
		self.annotToolBtnMenu.addAction("Testings")
		self.annotToolBtn.setMenu(self.annotToolBtnMenu)
		self.toolBar.addAction(self.annotToolBtn)

		self.statToolBtn = QAction(QIcon("icons/statistics.png"), self.tr("Analysis"), self)
		#self.statToolBtn.setDisabled(True)
		self.statToolBtnMenu = QMenu()
		self.statToolBtnMenu.addAction("Settings")
		self.statToolBtnMenu.addAction("Testings")
		self.statToolBtn.setMenu(self.statToolBtnMenu)
		self.toolBar.addAction(self.statToolBtn)

		self.reportToolBtn = QAction(QIcon("icons/report.png"), self.tr("Statistics"), self)
		#self.reportToolBtn.setDisabled(True)
		self.reportToolBtnMenu = QMenu()
		self.reportToolBtnMenu.addAction("Settings")
		self.reportToolBtnMenu.addAction("Testings")
		self.reportToolBtn.setMenu(self.reportToolBtnMenu)
		self.toolBar.addAction(self.reportToolBtn)

	def createStatusBar(self):
		self.statusBar = self.statusBar()
		self.statusBar.showMessage("Genome-wide microsatellites analysis tool.")
		self.progressBar = QProgressBar(self)
		self.progressBar.setAlignment(Qt.AlignRight)
		self.progressBar.setMaximumSize(140, 14)
		self.statusBar.addPermanentWidget(self.progressBar)




	def doOpenProject(self):
		for k in self.settings.allKeys():
			print k,"\t", self.settings.value(k)

	def doLoadFasta(self):
		fasta, _ = QFileDialog.getOpenFileName(self, filter="Fasta (*.fa *.fas *.fasta);;All files (*.*)")
		if fasta:
			self.inputFasta = fasta
			self.pushStatusMessage("Load fasta file %s" % fasta)

	def doLoadFastas(self):
		folder = QFileDialog.getExistingDirectory(self)
		if folder:
			self.inputFasta = folder
			self.pushStatusMessage("Load all fastas in folder %s" % folder)

	def doCopy(self):
		focus = QApplication.focusWidget()
		if focus is 0: return
		QApplication.postEvent(focus, QKeyEvent(QEvent.KeyPress, Qt.Key_C, Qt.ControlModifier))
		QApplication.postEvent(focus, QKeyEvent(QEvent.KeyRelease, Qt.Key_C, Qt.ControlModifier))

	def doCut(self):
		focus = QApplication.focusWidget()
		if focus is 0: return
		QApplication.postEvent(focus, QKeyEvent(QEvent.KeyPress, Qt.Key_X, Qt.ControlModifier))
		QApplication.postEvent(focus, QKeyEvent(QEvent.KeyRelease, Qt.Key_X, Qt.ControlModifier))
	
	def doPaste(self):
		focus = QApplication.focusWidget()
		if focus is 0: return
		QApplication.postEvent(focus, QKeyEvent(QEvent.KeyPress, Qt.Key_V, Qt.ControlModifier))
		QApplication.postEvent(focus, QKeyEvent(QEvent.KeyRelease, Qt.Key_V, Qt.ControlModifier))
	
	def doSelectAll(self):
		focus = QApplication.focusWidget()
		if focus is 0: return
		QApplication.postEvent(focus, QKeyEvent(QEvent.KeyPress, Qt.Key_A, Qt.ControlModifier))
		QApplication.postEvent(focus, QKeyEvent(QEvent.KeyRelease, Qt.Key_A, Qt.ControlModifier))

	def doPreference(self):
		dialog = PreferenceDialog(self, self.settings)
		if dialog.exec_() == QDialog.Accepted:
			dialog.saveSettings()

	def doSearchSSRs(self):
		#if self.db.isTableExists('ssr'):
		#	status = QMessageBox.warning(self, 
		#		self.tr("Warning"), 
		#		self.tr("The SSRs have been searched.\nDo you want to remove result and search again?"), 
		#		QMessageBox.Ok | QMessageBox.Cancel
		#	)

		#	if status == QMessageBox.Cancel:
			
		minimunRepeats = self.getMinimunRepeats()
		self.db.createSSRTable()
		if os.path.isfile(self.inputFasta):
			self.pushStatusMessage("Build fasta index for %s" % self.inputFasta)
			genome = pyfaidx.Fasta(self.inputFasta)

			self.db

	def getMinimunRepeats(self):
		return {
			1: int(self.settings.value('mono', 12)),
			2: int(self.settings.value('di', 7)),
			3: int(self.settings.value('tri', 5)), 
			4: int(self.settings.value('tetra', 4)),
			5: int(self.settings.value('penta', 3)),
			6: int(self.settings.value('hexa', 3))
		}
		

	def pushStatusMessage(self, msg):
		self.generalTab.writeMessage(msg)


class PreferenceDialog(QDialog):
	def __init__(self, parent=None, settings=None):
		super(PreferenceDialog, self).__init__(parent)
		self.settings = settings
		self.setWindowTitle(self.tr("Preferences"))
		self.setMinimumWidth(400)

		repeatsGroup = QGroupBox(self.tr("Minimum repeats"))
		monoLabel = QLabel("Mono-nucleotide")
		self.monoValue = QSpinBox()
		diLabel = QLabel("Di-nucleotide")
		self.diValue = QSpinBox()
		triLabel = QLabel("Tri-nucleotide")
		self.triValue = QSpinBox()
		tetraLabel = QLabel("Tetra-nucleotide")
		self.tetraValue = QSpinBox()
		pentaLabel = QLabel("Penta-nucleotide")
		self.pentaValue = QSpinBox()
		hexaLabel = QLabel("Hexa-nucleotide")
		self.hexaValue = QSpinBox()
		repeatLayout = QGridLayout()
		repeatLayout.setVerticalSpacing(10)
		repeatLayout.setHorizontalSpacing(10)
		repeatLayout.setColumnStretch(1, 1)
		repeatLayout.setColumnStretch(3, 1)
		repeatLayout.addWidget(monoLabel, 0, 0)
		repeatLayout.addWidget(self.monoValue, 0, 1)
		repeatLayout.addWidget(diLabel, 0, 2)
		repeatLayout.addWidget(self.diValue, 0, 3)
		repeatLayout.addWidget(triLabel, 1, 0)
		repeatLayout.addWidget(self.triValue, 1, 1)
		repeatLayout.addWidget(tetraLabel, 1, 2)
		repeatLayout.addWidget(self.tetraValue, 1, 3)
		repeatLayout.addWidget(pentaLabel, 2, 0)
		repeatLayout.addWidget(self.pentaValue, 2, 1)
		repeatLayout.addWidget(hexaLabel, 2, 2)
		repeatLayout.addWidget(self.hexaValue, 2, 3)
		repeatsGroup.setLayout(repeatLayout)

		distanceGroup = QGroupBox(self.tr("Compound SSR maximal distance"))
		distanceLabel = QLabel("Maximal distance (dmax): ")
		self.distanceValue = QSpinBox()
		distanceLayout = QHBoxLayout()
		distanceLayout.addWidget(distanceLabel)
		distanceLayout.addWidget(self.distanceValue, 1)
		distanceGroup.setLayout(distanceLayout)

		flankGroup = QGroupBox(self.tr("Flanking sequence"))
		flankLabel = QLabel("Flanking sequence length: ")
		self.flankValue = QSpinBox()
		self.flankValue.setMaximum(1000)
		flankLayout = QHBoxLayout()
		flankLayout.addWidget(flankLabel)
		flankLayout.addWidget(self.flankValue, 1)
		flankGroup.setLayout(flankLayout)

		buttonBox = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
		buttonBox.accepted.connect(self.accept)
		buttonBox.rejected.connect(self.reject)

		spacerItem = QSpacerItem(10, 20, QSizePolicy.Minimum, QSizePolicy.Expanding)

		mainLayout = QVBoxLayout()
		mainLayout.addWidget(repeatsGroup)
		mainLayout.addWidget(distanceGroup)
		mainLayout.addWidget(flankGroup)
		mainLayout.addItem(spacerItem)
		mainLayout.addWidget(buttonBox)
		self.setLayout(mainLayout)

		self.getSettings()

	def getSettings(self):
		self.monoValue.setValue(int(self.settings.value('mono', 12)))
		self.diValue.setValue(int(self.settings.value('di', 7)))
		self.triValue.setValue(int(self.settings.value('tri', 5)))
		self.tetraValue.setValue(int(self.settings.value('tetra', 4)))
		self.pentaValue.setValue(int(self.settings.value('penta', 3)))
		self.hexaValue.setValue(int(self.settings.value('hexa', 3)))
		self.distanceValue.setValue(int(self.settings.value('dmax', 10)))
		self.flankValue.setValue(int(self.settings.value('flank', 100)))


	def saveSettings(self):
		self.settings.setValue('mono', self.monoValue.value())
		self.settings.setValue('di', self.diValue.value())
		self.settings.setValue('tri', self.triValue.value())
		self.settings.setValue('tetra', self.tetraValue.value())
		self.settings.setValue('penta', self.pentaValue.value())
		self.settings.setValue('hexa', self.hexaValue.value())
		self.settings.setValue('dmax', self.distanceValue.value())
		self.settings.setValue('flank', self.flankValue.value())
		

class GeneralTab(QWidget):
	def __init__(self):
		super(GeneralTab, self).__init__()
		self.runningStatusWdg = QTextBrowser(self)
		layout = QVBoxLayout()
		layout.addWidget(QLabel("Running message and status:"))
		layout.addWidget(self.runningStatusWdg)
		self.setLayout(layout)

	def writeMessage(self, msg):
		self.runningStatusWdg.append(msg)
			
class SSRTab(QWidget):
	def __init__(self):
		super(SSRTab, self).__init__()
		#self.SSRTableWdg =


if __name__ == '__main__':
	app = QApplication(sys.argv)
	app.setOrganizationName('Mencent')
	app.setOrganizationDomain('mencent.com')
	app.setApplicationName('Gmia')
	mainWin = MainWindow()
	mainWin.show()
	sys.exit(app.exec_())