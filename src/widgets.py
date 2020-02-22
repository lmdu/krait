#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import csv
import apsw
import time
import json

import shutil
import requests
import platform

from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtSql import *
from PySide2.QtWebEngineWidgets import *
from PySide2.QtWidgets import *
from PySide2.QtPrintSupport import *

#from PySide.QtCore import *
#from PySide.QtGui import *
#from PySide.QtWebKit import *

from db import *
from utils import *
from detail import *
from workers import *
from config import *

class SSRMainWindow(QMainWindow):
	def __init__(self):
		super(SSRMainWindow, self).__init__()

		self.setWindowTitle("Krait v%s" % VERSION)
		#self.setWindowIcon(QIcon('icons/logo.png'))
		self.setWindowIcon(QIcon(':/icons/logo.ico'))

		#stacked widget
		self.main_widget = QStackedWidget(self)
		self.setCentralWidget(self.main_widget)
		
		self.table = SSRTableView(self)
		self.createTableModel()

		self.browser = QWebEngineView(self)
		#self.browser = BrowserWidget(self)
		#QWebSettings.globalSettings().setAttribute(QWebSettings.DeveloperExtrasEnabled, True);

		
		self.main_widget.addWidget(self.table)
		self.main_widget.addWidget(self.browser)

		#search text input
		self.filter = SSRFilterInput(self)
		self.filter.returnPressed.connect(self.filterTable)

		#create fasta table
		#self.fasta_table = FastaTable()

		self.createActions()
		self.createMenus()
		self.createToolBars()
		self.createStatusBar()

		#connect to database
		self.db = Database()

		#opened project
		self.opened_project = None

		#annotation file
		self.annot_file = None

		#statistical results
		self.statis_result = None

		#changed rows in database
		self.changed_rows = 0

		#read settings
		self.readSettings()

		#read home page
		#self.homepage()

		#Enable dragging and dropping onto the main window
		self.setAcceptDrops(True)

		self.show()


	def swichMainWidget(self, widget):
		if widget == 'table':
			self.exportTableAct.setEnabled(True)
			self.exportStatsAct.setDisabled(True)
			self.main_widget.setCurrentIndex(0)
			if self.model.tableName() in ('ssr', 'issr', 'cssr', 'vntr'):
				self.exportFastaAct.setEnabled(True)
				self.exportGFFAct.setEnabled(True)
			else:
				self.exportFastaAct.setDisabled(True)
				self.exportGFFAct.setDisabled(True)

		else:
			self.exportTableAct.setDisabled(True)
			self.exportStatsAct.setEnabled(True)
			self.main_widget.setCurrentIndex(1)
			self.exportFastaAct.setDisabled(True)
			self.exportGFFAct.setDisabled(True)

	def homepage(self):
		#content = template_render('index.html')
		#self.browser.setHtml(content, QUrl("qrc:/"))
		self.browser.load(QUrl('https://github.com/lmdu/krait'))

	def createTableModel(self):
		self.model = TableModel()
		self.table.setModel(self.model)
		self.model.row_col.connect(self.changeRowColCount)
		self.model.sel_row.connect(self.changeSelectCount)

	def readSettings(self):
		self.settings = QSettings(CONFIG_FILE, QSettings.IniFormat)
		if len(self.settings.allKeys()) == 0:
			dialog = PreferenceDialog(self, self.settings)
			dialog.saveSettings()
		self.resize(self.settings.value("size", QSize(1000, 600)))

	def writeSettings(self):
		if not self.isMaximized():
			self.settings.setValue("size", self.size())

	def closeEvent(self, event):
		self.writeSettings()
		if self.changed_rows == self.db.changes():
			event.accept()
			return

		ret = QMessageBox.question(self, "Closing", 
			"Would you like to save results before exiting",
			QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
		)

		if ret == QMessageBox.Yes:
			self.saveProject()
		elif ret == QMessageBox.Cancel:
			event.ignore()
		elif ret == QMessageBox.No:
			event.accept()

	def dragEnterEvent(self, event):
		if event.mimeData().hasUrls():
			event.accept()
		else:
			event.ignore()

	def dragMoveEvent(self, event):
		if event.mimeData().hasUrls():
			event.accept()
		else:
			event.ignore()

	def dropEvent(self, event):
		if not event.mimeData().hasUrls():
			return

		urls = event.mimeData().urls()
		if len(urls) > 1:
			return

		drag_file = urls[0].toLocalFile()

		if drag_file.endswith('.kdb'):
			self.openProject(drag_file)
			return

		for suffix in ['.fa', '.fna', '.fas', '.fasta', '.fa.gz', '.fna.gz', '.fas.gz', '.fasta.gz']:
			if drag_file.endswith(suffix):
				self.importFasta(drag_file)
				return
		warn_msg = (
			"Please drag a project file with .kdb suffix to open project, or "
			"Drag a fasta formatted sequence file with .fa, .fna, .fas, .fasta, "
			".fa.gz, .fna.gz, .fas.gz, .fasta.gz suffix to search repeats"
		)
		QMessageBox.warning(self, "Input file format not right", warn_msg)

	def createActions(self):
		#open a project action
		self.openProjectAct = QAction(self.tr("Open Project"), self)
		self.openProjectAct.setShortcut(QKeySequence.Open)
		self.openProjectAct.triggered.connect(self.openProject)

		#close a project action
		self.closeProjectAct = QAction(self.tr("Close Project"), self)
		self.closeProjectAct.setShortcut(QKeySequence.Close)
		self.closeProjectAct.triggered.connect(self.closeProject)
		
		#save a project action
		self.saveProjectAct = QAction(self.tr("Save Project"), self)
		self.saveProjectAct.setShortcut(QKeySequence.Save)
		self.saveProjectAct.triggered.connect(self.saveProject)
		
		#save as a project action
		self.saveAsProjectAct = QAction(self.tr("Save Project As..."), self)
		self.saveAsProjectAct.setShortcut(QKeySequence.SaveAs)
		self.saveAsProjectAct.triggered.connect(self.saveProjectAs)
		
		#load fasta file or genome action
		self.loadFastaAct = QAction(self.tr("Import Fasta"), self)
		self.loadFastaAct.triggered.connect(self.importFasta)
		self.loadFastaAct.setShortcut(QKeySequence(Qt.CTRL+Qt.SHIFT+Qt.Key_O))
		self.loadFastasAct = QAction(self.tr("Import Fastas in Folder"), self)
		self.loadFastasAct.triggered.connect(self.importFastas)
		
		#export the Results
		self.exportWholeTableAct = QAction(self.tr("Export Whole Table"), self)
		self.exportWholeTableAct.triggered.connect(self.exportWholeTable)
		self.exportWholeTableAct.setDisabled(True)
		self.exportSelectedRowsAct = QAction(self.tr("Export Selected Rows"), self)
		self.exportSelectedRowsAct.triggered.connect(self.exportSelectedRows)
		self.exportSelectedRowsAct.setDisabled(True)

		self.exportTableAct = QAction(self.tr("Exprot Selected as Table"), self)
		self.exportTableAct.triggered.connect(self.exportTableRows)
		self.exportTableAct.setDisabled(True)
		self.exportTableAct.setShortcut(QKeySequence(Qt.CTRL+Qt.SHIFT+Qt.Key_T))
		self.exportFastaAct = QAction(self.tr("Export Selected as Fasta"), self)
		self.exportFastaAct.triggered.connect(self.exportTableFastas)
		self.exportFastaAct.setDisabled(True)
		self.exportGFFAct = QAction(self.tr("Export Selected as GFF"), self)
		self.exportGFFAct.triggered.connect(self.exportTableGFF)
		self.exportGFFAct.setDisabled(True)
		self.exportStatsAct = QAction(self.tr("Export Statistical Report"), self)
		self.exportStatsAct.setShortcut(QKeySequence(Qt.CTRL+Qt.Key_P))
		self.exportStatsAct.setDisabled(True)
		self.exportStatsAct.triggered.connect(self.exportStatisResult)
		
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
		self.preferenceAct.triggered.connect(self.setPreference)

		#toolbar actions
		#search perfect ssrs tool button
		self.SSRSearchAct = QAction(QIcon(":/icons/ssr.png"), self.tr("SSRs"), self)
		self.SSRSearchAct.setToolTip(self.tr("Search for Perfect SSRs"))
		self.SSRSearchAct.setStatusTip(self.tr("Search for Perfect Microsatellites"))
		self.SSRSearchAct.triggered.connect(self.searchOrShowSSR)
		self.SSRForceAct = QAction(self.tr("Redo Search for SSRs"), self)
		self.SSRForceAct.setShortcut(QKeySequence(Qt.CTRL+Qt.Key_1))
		self.SSRForceAct.triggered.connect(self.searchSSR)
		self.SSRShowAct = QAction(self.tr("Show Perfect SSRs"), self)
		self.SSRShowAct.setShortcut(QKeySequence(Qt.CTRL+Qt.SHIFT+Qt.Key_1))
		self.SSRShowAct.triggered.connect(self.showSSR)
		self.SSRRemoveAct = QAction(self.tr("Remove Perfect SSRs"), self)
		self.SSRRemoveAct.triggered.connect(self.removeSSR)
		self.SSRSetAct = QAction(self.tr("Specify Minimum Repeats"), self)
		self.SSRSetAct.triggered.connect(self.setPreference)
		
		#search compound ssrs tool button
		self.CSSRSearchAct = QAction(QIcon(":/icons/cssr.png"), self.tr("cSSRs"), self)
		self.CSSRSearchAct.setToolTip(self.tr("Search for Compound SSRs"))
		self.CSSRSearchAct.setStatusTip(self.tr("Search for Compound Microsatellites"))
		self.CSSRSearchAct.triggered.connect(self.searchOrShowCSSR)
		self.CSSRForceAct = QAction(self.tr("Redo Search for cSSRs"), self)
		self.CSSRForceAct.setShortcut(QKeySequence(Qt.CTRL+Qt.Key_2))
		self.CSSRForceAct.triggered.connect(self.searchCSSR)
		self.CSSRShowAct = QAction(self.tr("Show Compound SSRs"), self)
		self.CSSRShowAct.setShortcut(QKeySequence(Qt.CTRL+Qt.SHIFT+Qt.Key_2))
		self.CSSRShowAct.triggered.connect(self.showCSSR)
		self.CSSRRemoveAct = QAction(self.tr("Remove Compound SSRs"), self)
		self.CSSRRemoveAct.triggered.connect(self.removeCSSR)
		#self.bestDmaxAct = QAction(self.tr("Estimate best dMax"), self)
		#self.bestDmaxAct.triggered.connect(self.estimateBestMaxDistance)
		self.CSSRSetAct = QAction(self.tr("Specify Maximal Distance"), self)
		self.CSSRSetAct.triggered.connect(self.setPreference)

		#search VNTRs
		self.VNTRSearchAct = QAction(QIcon(":/icons/vntr.png"), self.tr("VNTRs"), self)
		self.VNTRSearchAct.setToolTip(self.tr("Search for Minisatellites or Macrosatellites"))
		self.VNTRSearchAct.setStatusTip(self.tr("Search for Minisatellites or Macrosatellites"))
		self.VNTRSearchAct.triggered.connect(self.searchOrShowVNTR)
		self.VNTRForceAct = QAction(self.tr("Redo Search for VNTRs"), self)
		self.VNTRForceAct.setShortcut(QKeySequence(Qt.CTRL+Qt.Key_4))
		self.VNTRForceAct.triggered.connect(self.searchVNTR)
		self.VNTRShowAct = QAction(self.tr("Show VNTRs"), self)
		self.VNTRShowAct.setShortcut(QKeySequence(Qt.CTRL+Qt.SHIFT+Qt.Key_4))
		self.VNTRShowAct.triggered.connect(self.showVNTR)
		self.VNTRRemoveAct = QAction(self.tr("Remove VNTRs"), self)
		self.VNTRRemoveAct.triggered.connect(self.removeVNTR)
		self.VNTRSetAct = QAction(self.tr("Specify Search Parameters"), self)
		self.VNTRSetAct.triggered.connect(self.setPreference)

		#search imperfect microsatellites
		self.ISSRSearchAct = QAction(QIcon(":/icons/issr.png"), self.tr("iSSRs"), self)
		self.ISSRSearchAct.setToolTip(self.tr("Search for Imperfect SSRs"))
		self.ISSRSearchAct.setStatusTip(self.tr("Search for Imperfect Microsatellites"))
		self.ISSRSearchAct.triggered.connect(self.searchOrShowISSR)
		self.ISSRForceAct = QAction(self.tr("Redo Search for iSSRs"), self)
		self.ISSRForceAct.setShortcut(QKeySequence(Qt.CTRL+Qt.Key_3))
		self.ISSRForceAct.triggered.connect(self.searchISSR)
		self.ISSRShowAct = QAction(self.tr("Show Imperfect SSRs"), self)
		self.ISSRShowAct.setShortcut(QKeySequence(Qt.CTRL+Qt.SHIFT+Qt.Key_3))
		self.ISSRShowAct.triggered.connect(self.showISSR)
		self.ISSRRemoveAct = QAction(self.tr("Remove Imperfect SSRs"), self)
		self.ISSRRemoveAct.triggered.connect(self.removeISSR)
		self.ISSRSetAct = QAction(self.tr("Specify Search Parameters"), self)
		self.ISSRSetAct.triggered.connect(self.setPreference)

		#locate ssrs
		self.locateAct = QAction(QIcon(":/icons/locate.png"), self.tr("Mapping"), self)
		self.locateAct.setToolTip(self.tr("Mapping tandem repeats to genic regions"))
		self.locateAct.setStatusTip(self.tr("Mapping tandem repeats to genic regions"))
		#self.locateAct.triggered.connect(self.LocateOrShowTandem)
		self.locateAct.triggered.connect(self.locateTandem)
		
		self.locateToolAct = QAction(self.tr("Redo Mapping"), self)
		self.locateToolAct.triggered.connect(self.locateTandem)
		#self.locateShowAct = QAction(self.tr("Show mapping reslult"), self)
		#self.locateShowAct.triggered.connect(self.showLocation)

		self.locateSetAct = QAction(self.tr("Import Annotation File"), self)
		self.locateSetAct.triggered.connect(self.provideAnnotation)
		self.removeLocateAct = QAction(self.tr("Remove mapping reslult"), self)
		self.removeLocateAct.triggered.connect(self.removeMarker)
		
		cds_icon = QPixmap(16, 16)
		cds_icon.fill(QColor(245, 183, 177))
		self.showCDSAct = QAction(QIcon(cds_icon), self.tr("Show repeats in CDS"), self)
		self.showCDSAct.triggered.connect(self.showCDSMarker)

		exon_icon = QPixmap(16, 16)
		exon_icon.fill(QColor(169, 223, 191))
		self.showExonAct = QAction(QIcon(exon_icon), self.tr("Show repeats in exon"), self)
		self.showExonAct.triggered.connect(self.showExonMarker)

		utr_icon = QPixmap(16, 16)
		utr_icon.fill(QColor(250, 215, 160))
		self.showUTRAct = QAction(QIcon(utr_icon), self.tr("Show repeats in UTR"), self)
		self.showUTRAct.triggered.connect(self.showUTRMarker)

		intron_icon = QPixmap(16, 16)
		intron_icon.fill(QColor(174, 214, 241))
		self.showIntronAct = QAction(QIcon(intron_icon), self.tr("Show repeats in intron"), self)
		self.showIntronAct.triggered.connect(self.showIntronMarker)


		#design primer
		self.primerDesignAct = QAction(QIcon(":/icons/primer.png"), self.tr("Primer"), self)
		self.primerDesignAct.setToolTip(self.tr("Design primers"))
		self.primerDesignAct.setStatusTip(self.tr("Start design primers"))
		self.primerDesignAct.triggered.connect(self.designOrShowPrimer)
		self.primerForceAct = QAction(self.tr("Redo Design Primers"), self)
		self.primerForceAct.triggered.connect(self.designPrimer)
		self.primerShowAct = QAction(self.tr("Show Designed Primer"), self)
		self.primerShowAct.triggered.connect(self.showPrimer)
		self.primerRemoveAct = QAction(self.tr("Remove Designed Primer"), self)
		self.primerRemoveAct.triggered.connect(self.removePrimer)
		self.primerSetAct = QAction(self.tr("Specify Primer3 Settings"), self)
		self.primerSetAct.triggered.connect(self.setPrimerSettings)

		#statistics report
		self.statisticsAct = QAction(QIcon(":/icons/report.png"), self.tr("Statistics"), self)
		self.statisticsAct.triggered.connect(self.doOrShowStatistics)
		self.statisticsAct.setStatusTip(self.tr("Generate Statistical Report"))
		self.statisticsForceAct = QAction(self.tr("Redo Statistical Analysis"), self)
		self.statisticsForceAct.triggered.connect(self.performStatistics)
		self.statisticsShowAct = QAction(self.tr("Show Statistical Result"), self)
		self.statisticsShowAct.triggered.connect(self.showStatistics)
		self.statisticsRemoveAct = QAction(self.tr("Remove Statistical Result"), self)
		self.statisticsRemoveAct.triggered.connect(self.removeStatistics)

		#tool action
		self.downloadNCBIAct = QAction(self.tr("Download sequence from NCBI"), self)
		self.downloadNCBIAct.triggered.connect(self.downloadFasta)

		#about action
		self.aboutAct = QAction(self.tr("About Krait"), self)
		self.aboutAct.triggered.connect(self.openAboutMessage)

		#documentation action
		self.documentAct = QAction(self.tr("Documentation"), self)
		self.documentAct.setShortcut(QKeySequence.HelpContents)
		self.documentAct.triggered.connect(self.openDocumentation)

		#report issue action
		self.issueAct = QAction(self.tr("Report issue..."), self)
		self.issueAct.triggered.connect(self.reportIssue)

		#view input action
		self.showInputAct = QAction(self.tr("Show Input Files"), self)
		self.showInputAct.triggered.connect(self.showInputFastas)
		

	def createMenus(self):
		self.fileMenu = self.menuBar().addMenu("&File")
		self.editMenu = self.menuBar().addMenu("&Edit")
		self.searchMenu = self.menuBar().addMenu("&Search")
		self.viewMenu = self.menuBar().addMenu("&View")
		self.toolMenu = self.menuBar().addMenu("&Tool")
		self.helpMenu = self.menuBar().addMenu("&Help")
		
		self.fileMenu.addAction(self.openProjectAct)
		self.fileMenu.addAction(self.closeProjectAct)
		self.fileMenu.addSeparator()
		self.fileMenu.addAction(self.saveProjectAct)
		self.fileMenu.addAction(self.saveAsProjectAct)
		self.fileMenu.addSeparator()
		self.fileMenu.addAction(self.loadFastaAct)
		self.fileMenu.addAction(self.loadFastasAct)
		self.fileMenu.addAction(self.locateSetAct)
		self.fileMenu.addSeparator()
		self.fileMenu.addAction(self.exportWholeTableAct)
		self.fileMenu.addAction(self.exportSelectedRowsAct)
		#self.fileMenu.addAction(self.exportTableAct)
		#self.fileMenu.addAction(self.exportFastaAct)
		#self.fileMenu.addAction(self.exportGFFAct)
		self.fileMenu.addAction(self.exportStatsAct)
		self.fileMenu.addSeparator()
		self.fileMenu.addAction(self.exitAct)
		
		self.editMenu.addAction(self.copyAct)
		self.editMenu.addAction(self.cutAct)
		self.editMenu.addAction(self.pasteAct)
		self.editMenu.addSeparator()
		self.editMenu.addAction(self.selectAllAct)
		self.editMenu.addSeparator()
		self.editMenu.addAction(self.preferenceAct)

		self.searchMenu.addAction(self.SSRForceAct)
		self.searchMenu.addAction(self.CSSRForceAct)
		self.searchMenu.addAction(self.ISSRForceAct)
		self.searchMenu.addAction(self.VNTRForceAct)

		self.viewMenu.addAction(self.showInputAct)
		self.viewMenu.addSeparator()
		self.viewMenu.addAction(self.SSRShowAct)
		self.viewMenu.addAction(self.CSSRShowAct)
		self.viewMenu.addAction(self.ISSRShowAct)
		self.viewMenu.addAction(self.VNTRShowAct)
		self.viewMenu.addSeparator()
		self.viewMenu.addAction(self.SSRRemoveAct)
		self.viewMenu.addAction(self.CSSRRemoveAct)
		self.viewMenu.addAction(self.ISSRRemoveAct)
		self.viewMenu.addAction(self.VNTRRemoveAct)

		#self.toolMenu.addAction(self.bestDmaxAct)
		self.toolMenu.addAction(self.primerForceAct)
		self.toolMenu.addAction(self.locateToolAct)
		self.toolMenu.addAction(self.statisticsForceAct)
		self.toolMenu.addSeparator()
		self.toolMenu.addAction(self.downloadNCBIAct)

		self.helpMenu.addAction(self.documentAct)
		self.helpMenu.addAction(self.issueAct)
		self.helpMenu.addSeparator()
		self.helpMenu.addAction(self.aboutAct)


		#tool bar menus
		#search ssrs tool button menu
		self.SSRMenu = QMenu()
		self.SSRMenu.addAction(self.SSRForceAct)
		self.SSRMenu.addAction(self.SSRShowAct)
		self.SSRMenu.addSeparator()
		self.SSRMenu.addAction(self.SSRRemoveAct)
		self.SSRMenu.addSeparator()
		self.SSRMenu.addAction(self.loadFastaAct)
		self.SSRMenu.addAction(self.loadFastasAct)
		self.SSRMenu.addSeparator()
		self.SSRMenu.addAction(self.SSRSetAct)

		self.CSSRMenu = QMenu()
		self.CSSRMenu.addAction(self.CSSRForceAct)
		self.CSSRMenu.addAction(self.CSSRShowAct)
		self.CSSRMenu.addSeparator()
		self.CSSRMenu.addAction(self.CSSRRemoveAct)
		self.CSSRMenu.addSeparator()
		#self.CSSRMenu.addAction(self.bestDmaxAct)
		self.CSSRMenu.addAction(self.CSSRSetAct)

		self.VNTRMenu = QMenu()
		self.VNTRMenu.addAction(self.VNTRForceAct)
		self.VNTRMenu.addAction(self.VNTRShowAct)
		self.VNTRMenu.addSeparator()
		self.VNTRMenu.addAction(self.VNTRRemoveAct)
		self.VNTRMenu.addSeparator()
		self.VNTRMenu.addAction(self.VNTRSetAct)

		self.ISSRMenu = QMenu()
		self.ISSRMenu.addAction(self.ISSRForceAct)
		self.ISSRMenu.addAction(self.ISSRShowAct)
		self.ISSRMenu.addSeparator()
		self.ISSRMenu.addAction(self.ISSRRemoveAct)
		self.ISSRMenu.addSeparator()
		self.ISSRMenu.addAction(self.ISSRSetAct)

		self.locateMenu = QMenu()
		self.locateMenu.addAction(self.locateSetAct)
		self.locateMenu.addAction(self.locateToolAct)
		#self.locateMenu.addAction(self.locateShowAct)
		self.locateMenu.addSeparator()
		self.locateMenu.addAction(self.removeLocateAct)
		self.locateMenu.addSeparator()
		self.locateMenu.addAction(self.showCDSAct)
		self.locateMenu.addAction(self.showExonAct)
		self.locateMenu.addAction(self.showUTRAct)
		self.locateMenu.addAction(self.showIntronAct)

		self.primerMenu = QMenu()
		self.primerMenu.addAction(self.primerForceAct)
		self.primerMenu.addAction(self.primerShowAct)
		self.primerMenu.addSeparator()
		self.primerMenu.addAction(self.primerRemoveAct)
		self.primerMenu.addSeparator()
		self.primerMenu.addAction(self.primerSetAct)

		self.statisticsMenu = QMenu()
		self.statisticsMenu.addAction(self.statisticsForceAct)
		self.statisticsMenu.addAction(self.statisticsShowAct)
		self.statisticsMenu.addSeparator()
		self.statisticsMenu.addAction(self.statisticsRemoveAct)
		

	def createToolBars(self):
		self.toolBar = self.addToolBar('')
		self.toolBar.setMovable(False)
		self.toolBar.setIconSize(QSize(36, 36))
		self.toolBar.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)

		#search ssr action and menus
		self.SSRSearchAct.setMenu(self.SSRMenu)
		self.toolBar.addAction(self.SSRSearchAct)

		self.CSSRSearchAct.setMenu(self.CSSRMenu)
		self.toolBar.addAction(self.CSSRSearchAct)

		self.ISSRSearchAct.setMenu(self.ISSRMenu)
		self.toolBar.addAction(self.ISSRSearchAct)

		self.VNTRSearchAct.setMenu(self.VNTRMenu)
		self.toolBar.addAction(self.VNTRSearchAct)

		self.locateAct.setMenu(self.locateMenu)
		self.toolBar.addAction(self.locateAct)

		self.primerDesignAct.setMenu(self.primerMenu)
		self.toolBar.addAction(self.primerDesignAct)

		self.statisticsAct.setMenu(self.statisticsMenu)
		self.toolBar.addAction(self.statisticsAct)

		#self.reportToolBtn.setDisabled(True)

		#search input
		#self.filter.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
		self.toolBar.addWidget(self.filter)

	def createStatusBar(self):
		self.statusBar = self.statusBar()
		self.statusBar.showMessage("Genome-wide microsatellites analysis tool.")
		
		#add row and column counts widget
		self.rowCounts = QLabel("Row: 0", self)
		self.rowCounts.setStyleSheet("margin-right:10px;")
		self.statusBar.addPermanentWidget(self.rowCounts)
		
		self.colCounts = QLabel("Column: 0", self)
		self.colCounts.setStyleSheet("margin-right:10px;")
		self.statusBar.addPermanentWidget(self.colCounts)

		self.selectCounts = QLabel("Select: 0", self)
		self.selectCounts.setStyleSheet("margin-right:10px;")
		self.statusBar.addPermanentWidget(self.selectCounts)
		
		#add progressing bar
		self.progressBar = QProgressBar(self)
		self.progressBar.setMaximum(100)
		self.progressBar.setMinimum(0)
		self.statusBar.addPermanentWidget(self.progressBar)

	def showInputFastas(self, event):
		inputs = []
		
		if self.opened_project:
			inputs.append("<b>Input project:</b><p>{}</p>".format(self.opened_project))

		if self.annot_file:
			inputs.append("<b>Input annotation:</b><p>{}</p>".format(self.annot_file))

		fastas = self.db.get_all("SELECT * FROM fasta")

		if fastas:
			inputs.append("<b>Input Fastas:</b>")

		for fasta in fastas:
			inputs.append("<p>{}</p>".format(fasta[1]))

		if not inputs:
			QMessageBox.information(self, "Input files", "No input file found, import fasta file for analysis")
		else:
			QMessageBox.information(self, "Input files", "".join(inputs))

		
	def openProject(self, dbfile=None):
		if not dbfile:
			dbfile, _ = QFileDialog.getOpenFileName(self, filter="Krait Database (*.kdb)")
		
		if not dbfile:
			return

		self.opened_project = dbfile

		self.db.drop_tables()
		self.db.open(dbfile)

		if not self.db.is_empty('ssr'):
			self.showSSR()

		elif not self.db.is_empty('cssr'):
			self.showCSSR()

		elif not self.db.is_empty('issr'):
			self.showISSR()

		elif not self.db.is_empty('vntr'):
			self.showVNTR()

		self.setStatusMessage("%s has been successfully opened" % dbfile)

	def saveProject(self):
		if self.opened_project and not self.db.changes():
			return

		if self.opened_project is None:
			dbfile, _ = QFileDialog.getSaveFileName(self, filter="Krait Database (*.kdb)")
			if not dbfile:
				return
			self.opened_project = dbfile
		
		self.db.save(self.opened_project)
		worker = SaveProjectWorker(self.opened_project)
		self.executeTask(worker, lambda: 1)
		
		self.changed_rows = self.db.changes()
		

	def saveProjectAs(self):
		dbfile, _ = QFileDialog.getSaveFileName(self, filter="Krait Database (*.kdb)")
		if not dbfile:
			return

		self.changed_rows = self.db.changes()
		worker = SaveProjectWorker(dbfile)
		self.executeTask(worker, lambda: 1)
		self.opened_project = dbfile

	def closeProject(self):
		if self.changed_rows != self.db.changes():
			ret = QMessageBox.question(self, "Closing", 
				"Would you like to save results before exiting",
				QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
			)

			if ret == QMessageBox.Cancel:
				return

			if ret == QMessageBox.Yes:
				self.saveProject()
		
		self.db.drop_tables()
		self.createTableModel()
		self.browser.setHtml('')
		self.changed_rows = self.db.changes()

	def importFasta(self, fasta=None):
		'''
		Import a fasta file from a directory
		'''
		if not self.db.is_empty('fasta'):
			ret = QMessageBox.warning(self, "Closing", 
				"Are you sure you want to import new fasta file?<br>This will remove the previously imported fasta file.",
				QMessageBox.Yes | QMessageBox.No
			)

			if ret == QMessageBox.No:
				return

		if not fasta:
			fasta, _ = QFileDialog.getOpenFileName(self, filter="Fasta (*.fa *.fna *.fas *.fasta *.fna.gz *.fa.gz *.fasta.gz);;All files (*.*)")
	
		if not fasta:
			return

		self.db.clear('fasta')

		self.db.get_cursor().execute('INSERT INTO fasta VALUES (?,?)', (None, fasta))
		self.setStatusMessage("Import fasta %s" % fasta)

	def importFastas(self):
		'''
		import all fasta files from a directory
		'''
		if not self.db.is_empty('fasta'):
			ret = QMessageBox.warning(self, "Closing", 
				"Are you sure you want to import new fasta file?<br>This will remove the previously imported fasta file.",
				QMessageBox.Yes | QMessageBox.No
			)

			if ret == QMessageBox.No:
				return

		directory = QFileDialog.getExistingDirectory(self)
		if not directory:
			return
		
		folder = QDir(directory)
		fastas = folder.entryList(QDir.Files)

		if not fastas:
			return
		else:
			self.db.clear('fasta')
		
		count = 0
		for fasta in fastas:
			count += 1
			self.db.get_cursor().execute("INSERT INTO fasta VALUES (?,?)", (None, folder.absoluteFilePath(fasta)))
		
		self.setStatusMessage("Import %s fastas in %s" % (count, directory))


	def downloadFasta(self):
		'''
		download fasta file from NCBI database
		'''
		dialog = DownloadDialog(self)
		status = dialog.exec_()
		if not status: return
		acc, out = dialog.get()
		if not acc or not out:
			return
		self.progressBar.setMaximum(0)
		worker = EutilWorker(acc, out)
		self.executeTask(worker, lambda: self.progressBar.setMaximum(100))

	def exportTable(self, selected):
		suffix = ['Tabular (*.tsv)', 'CSV (*.csv)']
		if self.main_widget.currentWidget() == self.table:
			if self.model.tableName() in ('ssr', 'issr', 'cssr', 'vntr'):
				suffix.append('Fasta (*.fa);;GFF3 (*.gff)')
		suffix.append('TXT (*.txt)')

		outfile, _ = QFileDialog.getSaveFileName(self, filter=";;".join(suffix))
		if not outfile:
			return

		if self.model.tableName() == 'primer':
			worker = ExportPrimerWorker(selected, self.model, outfile)

		elif self.model.tableName() == 'feature':
			worker = ExportFeatureWorker(selected, self.model, outfile)

		elif outfile.endswith('.fa'):
			flank = int(self.settings.value('ssr/flank', 100))
			worker = ExportFastaWorker(selected, self.model, flank, outfile)
		
		else:
			worker = ExportTableWorker(selected, self.model, outfile)


		#worker.process()
		self.executeTask(worker, 
			lambda: QMessageBox.information(self, "Export Successed", "Successfully exported to {}".format(outfile))
		)

	def exportWholeTable(self):
		self.exportTable('whole')

	def exportSelectedRows(self):
		if not self.model.selected:
			return QMessageBox.warning(self, 'Warning', "No rows in table have been selected to export.")

		self.exportTable('selecte')

	def exportTableRows(self):
		#selected = self.model.getSelectedRows()
		if not self.model.selected:
			return QMessageBox.warning(self, 'Warning', "Please select rows in table to export.")

		exp_file, _ = QFileDialog.getSaveFileName(self, filter="CSV (*.csv);;Tabular text (*.txt)")
		if not exp_file:
			return

		worker = ExportTableWorker(self.model, exp_file)
		self.executeTask(worker,
			lambda: QMessageBox.information(self, "Information", "Successfully exported to %s" % exp_file)
		)


	def exportTableGFF(self):
		if not self.model.selected:
			return QMessageBox.warning(self, 'Warning', "Please select rows in table to export.")

		exp_file, _ = QFileDialog.getSaveFileName(self, filter="GFF3 (*.gff)")
		if not exp_file:
			return

		worker = ExportTableWorker(self.model, exp_file)
		self.executeTask(worker, 
			lambda: QMessageBox.information(self, "Information", "Successfully exported to %s" % exp_file)
		)
		

	def exportTableFastas(self):
		if not self.model.selected:
			return QMessageBox.warning(self, 'Warning', "Please select rows in table to export.")

		#if table not in ('ssr', 'issr', 'cssr', 'vntr'):
		#	return QMessageBox.warning(self, 'Warning', "Your selected rows are not SSRs")

		exp_file, _ = QFileDialog.getSaveFileName(self, filter="Fasta (*.fa);;Fasta (*.fasta)")
		if not exp_file:
			return

		flank = int(self.settings.value('ssr/flank', 100))
		
		worker = ExportFastaWorker(self.model, flank, exp_file)
		self.executeTask(worker,
			lambda: QMessageBox.information(self, "Information", "Successfully exported to %s" % exp_file)
		)

	def exportStatisResult(self):
		pdfname, _ = QFileDialog.getSaveFileName(self, filter="HTML Report (*.html);; PDF Report (*.pdf)")
		
		if not pdfname: return

		if pdfname.endswith('.pdf'):
			#page_layout = QPageLayout(
			#)
			self.browser.page().printToPdf(pdfname, QPageLayout(QPageSize(QPageSize.A4), QPageLayout.Portrait, QMarginsF(32,32,32,32)))
			#printer = QPrinter(QPrinter.HighResolution)
			#printer.setPageSize(QPrinter.A4)
			#printer.setColorMode(QPrinter.Color)
			#printer.setOutputFormat(QPrinter.PdfFormat)
			#printer.setOutputFileName(pdfname)
			#self.browser.print_(printer)

		elif pdfname.endswith('.html'):
			content = self.browser.page().save(pdfname, QWebEngineDownloadItem.SingleHtmlSaveFormat)
			#content = self.browser.page().mainFrame().toHtml()
			#with open(pdfname, 'w') as fw:
			#	fw.write(content)

		self.setStatusMessage("Statistical report was successfully save to %s" % pdfname)


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

	def setPreference(self):
		dialog = PreferenceDialog(self, self.settings)
		if dialog.exec_() == QDialog.Accepted:
			dialog.saveSettings()

	def setPrimerSettings(self):
		dialog = PreferenceDialog(self, self.settings)
		dialog.gotoPrimer()
		if dialog.exec_() == QDialog.Accepted:
			dialog.saveSettings()


	def executeTask(self, worker, finish_callback):
		#check the running task
		if hasattr(self, 'work_thread') and self.work_thread.isRunning():
			return QMessageBox.warning(self, "Warning", "Task is running! Please wait until finished.")

		self.worker = worker
		self.worker.update_message.connect(self.setStatusMessage)
		self.worker.update_progress.connect(self.setProgress)
		self.worker.error_message.connect(lambda x: QMessageBox.critical(self, 'Error', x))
		self.worker.finished.connect(finish_callback)
		self.work_thread = QThread()
		self.work_thread.started.connect(self.worker.run)
		self.worker.finished.connect(self.work_thread.quit)
		self.worker.finished.connect(self.worker.deleteLater)
		self.worker.failed.connect(self.work_thread.quit)
		self.worker.failed.connect(self.worker.deleteLater)
		self.worker.moveToThread(self.work_thread)
		self.work_thread.start()

	def getInputFastas(self):
		'''
		get all fasta files with path provided by user
		'''
		if self.db.is_empty('fasta'):
			QMessageBox.critical(self, "Error occurred", "No fasta file imported. Please import fasta files.")
			return None
		else:
			return self.db.get_all("SELECT * FROM fasta")


	#handle perfect SSRs search
	def searchSSR(self):
		fastas = self.getInputFastas()
		if not fastas:
			return

		self.removeSSR()
		
		rules = [
			int(self.settings.value('ssr/mono', 12)),
			int(self.settings.value('ssr/di', 7)),
			int(self.settings.value('ssr/tri', 5)), 
			int(self.settings.value('ssr/tetra', 4)),
			int(self.settings.value('ssr/penta', 4)),
			int(self.settings.value('ssr/hexa', 4))
		]
		level = int(self.settings.value('ssr/level', 3))
		worker = SSRWorker(fastas, rules, level)
		self.executeTask(worker, self.showSSR)

	def searchOrShowSSR(self):
		if self.db.is_empty('ssr'):
			self.searchSSR()
		else:
			self.showSSR()
	
	def showSSR(self):
		self.model.setTable('ssr')
		self.model.select()
		self.swichMainWidget('table')

	def removeSSR(self):
		self.model.remove('ssr')
	
	#handle compound SSRs search
	def searchCSSR(self):
		if self.db.is_empty('ssr'):
			QMessageBox.warning(self, "Warning", "Please search perfect SSRs first, before search compound SSRs.")
			return

		self.removeCSSR()

		dmax = int(self.settings.value('ssr/dmax', 10))
		worker = CSSRWorker(dmax)
		self.executeTask(worker, self.showCSSR)

	def searchOrShowCSSR(self):
		if self.db.is_empty('cssr'):
			self.searchCSSR()
		else:
			self.showCSSR()

	def showCSSR(self):
		self.model.setTable('cssr')
		self.model.select()
		self.swichMainWidget('table')

	def removeCSSR(self):
		self.model.remove('cssr')

	#handle VNTRs search
	def searchVNTR(self):
		fastas = self.getInputFastas()
		if not fastas:
			return

		self.removeVNTR()

		min_motif = int(self.settings.value('ssr/vmin', 7))
		max_motif = int(self.settings.value('ssr/vmax', 30))
		min_repeat = int(self.settings.value('ssr/vrep', 2))
		worker = VNTRWorker(fastas, min_motif, max_motif, min_repeat)
		self.executeTask(worker, self.showVNTR)

	def searchOrShowVNTR(self):
		if self.db.is_empty('vntr'):
			self.searchVNTR()
		else:
			self.showVNTR()
		
	def showVNTR(self):
		self.model.setTable('vntr')
		self.model.select()
		self.swichMainWidget('table')

	def removeVNTR(self):
		self.model.remove('vntr')
	
	#handle imperfect SSRs search
	def searchISSR(self):
		fastas = self.getInputFastas()
		if not fastas:
			return

		self.removeISSR()

		seed_repeat = int(self.settings.value('ssr/srep', 3))
		seed_length = int(self.settings.value('ssr/slen', 8))
		max_eidts = int(self.settings.value('ssr/error', 2))
		mis_penalty = int(self.settings.value('ssr/mismatch', 1))
		gap_penalty = int(self.settings.value('ssr/gap', 2))
		score = int(self.settings.value('ssr/score', 12))
		level = int(self.settings.value('ssr/level', 3))
		worker = ISSRWorker(fastas, seed_repeat, seed_length, max_eidts, mis_penalty, gap_penalty, score, level)
		self.executeTask(worker, self.showISSR)

	def searchOrShowISSR(self):
		if self.db.is_empty('issr'):
			self.searchISSR()
		else:
			self.showISSR()

	def showISSR(self):
		self.model.setTable('issr')
		self.model.select()
		self.swichMainWidget('table')

	def removeISSR(self):
		self.model.remove('issr')

	def getPrimerSettings(self):
		p3_settings = dict(
			PRIMER_TASK = 'generic',
			PRIMER_PICK_LEFT_PRIMER = 1,
			PRIMER_PICK_INTERNAL_OLIGO = 0,
			PRIMER_PICK_RIGHT_PRIMER = 1
		)

		repeat_library = ['',
			os.path.join(PRIMER3_CONFIG, 'humrep_and_simple.txt'),
			os.path.join(PRIMER3_CONFIG, 'rodent_ref.txt'),
			os.path.join(PRIMER3_CONFIG, 'rodrep_and_simple.txt'),
			os.path.join(PRIMER3_CONFIG, 'drosophila.w.transposons.txt')
		]

		self.settings.beginGroup("primer")
		keys = self.settings.childKeys()
		for key in keys:
			if key == 'PRIMER_PRODUCT_SIZE_RANGE':
				sgs = self.settings.value(key)
				p3_settings[key] = [list(map(int, sg.split('-'))) for sg in sgs.split()]
			else:
				p3_settings[key] = int(self.settings.value(key))
		self.settings.endGroup()

		p3_settings['PRIMER_MISPRIMING_LIBRARY'] = repeat_library[p3_settings['PRIMER_MISPRIMING_LIBRARY']]

		return p3_settings

	#design primers for ssrs
	def designPrimer(self):
		if not self.model.selected:
			return QMessageBox.warning(self, "Warning", "Please select tandem repeat makers")
		
		flank = int(self.settings.value('ssr/flank'))
		primer3_settings = self.getPrimerSettings()
		worker = PrimerWorker(self.model, flank, primer3_settings)

		self.removePrimer()
		#worker.process()
		#self.showPrimer()
		self.executeTask(worker, self.showPrimer)

	def designOrShowPrimer(self):
		if self.db.is_empty('primer'):
			self.designPrimer()
		else:
			self.showPrimer()

	def showPrimer(self):
		self.model.setTable('primer')
		self.model.select()
		self.swichMainWidget('table')

	def removePrimer(self):
		self.model.remove('primer')

	def provideAnnotation(self):
		#dialog = AnnotationDialog(self, self.annot_file)
		#dialog.exec_()
		#self.annot_file = dialog.get()
		filters = "GFF or GTF (*.gtf *.gtf.gz *.gff *.gff3 *.gff.gz *.gff3.gz);;ALL (*.*)"
		annot_file, _ = QFileDialog.getOpenFileName(self, "GFF or GTF annotation file", filter=filters)
		if not annot_file:
			return
		self.annot_file = annot_file
		self.setStatusMessage("Import annotation file %s" % self.annot_file)


	#def LocateOrShowTandem(self):
	#	if self.db.is_empty('feature'):
	#		self.locateTandem()
	#	else:
	#		self.showLocation()
	
	def locateTandem(self):
		if not self.annot_file:
			return QMessageBox.warning(self, "Warning", "Please provide gtf or gff well formated annotation file")

		table = self.model.table
		
		if not table:
			return QMessageBox.warning(self, "Warning", "You have not yet detected SSRs")

		if table not in ['ssr', 'issr', 'cssr', 'vntr']:
			return QMessageBox.warning(self, "Warning", "No SSRs displayed")

		if self.db.is_empty(table):
			return QMessageBox.warning(self, "warning", "No SSRs in table")

		worker = LocateWorker(table, self.annot_file)
		#worker.process()
		#self.showLocation()
		self.executeTask(worker, self.showLocation)

	def showLocation(self):
	#	self.model.setTable('feature')
	#	self.model.select()
	#	self.swichMainWidget('table')
		self.model.beginResetModel()
		self.model.endResetModel()

	def removeMarker(self):
		self.model.remove('location')

	def showMarker(self, marker):
		repeat_types = {'ssr': 1, 'cssr': 2, 'issr': 3, 'vntr': 4}
		feature_maps = {'CDS': 1, 'exon': 2, 'UTR': 3, 'intron': 4}
		
		table = self.model.table

		if table not in repeat_types:
			return QMessageBox.warning(self, 'Warning', "No tandem repeats displayed")
		
		#if not table or table == 'primer':
		#	return

		if feature_maps[marker] == 3:
			sql = "SELECT target FROM location WHERE reptype={} AND (feature=3 OR feature=5)".format(repeat_types[table])
		else:
			sql = "SELECT target FROM location WHERE reptype={} AND feature={}".format(repeat_types[table], feature_maps[marker])

		#data = self.db.get_column(sql)
		#if not data:
		#	return QMessageBox.warning(self, "Warning", "No %ss located in %s region" % (table.upper(), marker))
		
		self.model.setFilter('id IN (%s)' % sql)

	def showCDSMarker(self):
		self.showMarker('CDS')

	def showExonMarker(self):
		self.showMarker('exon')

	def showUTRMarker(self):
		self.showMarker('UTR')

	def showIntronMarker(self):
		self.showMarker('intron')

	def filterTable(self):
		filters = str(self.filter.text())
		if 'table=' in filters:
			self.model.setTable(filters.split('=')[1])
			self.model.select()
			return
		sqlwhere = format_sql_where(filters)
		self.model.setFilter(sqlwhere)

	def performStatistics(self):
		if self.db.is_empty('fasta'):
			return QMessageBox.warning(self, "Warning", "No fasta file inputted")

		worker = StatisWorker()
		#worker.process()
		self.progressBar.setMaximum(0)
		self.executeTask(worker, self.showStatistics)

	def doOrShowStatistics(self):
		if self.statis_result:
			self.swichMainWidget('browser')
			self.browser.setHtml(self.statis_result, QUrl("qrc:/"))
		else:
			self.performStatistics()

	def showStatistics(self):
		seq_statis = None
		ssr_statis = None
		issr_statis = None
		cssr_statis = None
		vntr_statis = None

		opt = self.db.get_option('seq_statis')
		if opt is not None:
			seq_statis = json.loads(opt)

		opt = self.db.get_option('ssr_statis')
		if opt is not None:
			ssr_statis = json.loads(opt)

		opt = self.db.get_option('issr_statis')
		if opt is not None:
			issr_statis = json.loads(opt)

		opt = self.db.get_option('cssr_statis')
		if opt is not None:
			cssr_statis = json.loads(opt)

		opt = self.db.get_option('vntr_statis')
		if opt is not None:
			vntr_statis = json.loads(opt)

		self.statis_result = template_render('report.html', 
			seq = seq_statis, 
			ssr = ssr_statis, 
			issr = issr_statis, 
			cssr = cssr_statis, 
			vntr = vntr_statis
		)

		self.swichMainWidget('browser')
		self.browser.setHtml(self.statis_result, QUrl("qrc:/"))

		self.progressBar.setMaximum(100)

	def removeStatistics(self):
		if self.statis_result:
			self.statis_result = None
			self.browser.setHtml('')

	def showSSRSequence(self, index):
		'''
		The row in table double clicked, show the sequence of SSR
		'''
		table = self.model.tableName()
		if table not in ['ssr', 'cssr']:
			return
		
		flank = int(self.settings.value('flank', 50))
		record = self.model.record(index.row())
		ssr = Data({record.fieldName(i) : record.value(i) for i in range(record.count())})
		sql = "SELECT f.path FROM fasta AS f, sequence AS s WHERE f.fid=s.fid AND s.name='%s'"
		query = QSqlQuery(sql % ssr.sequence)
		query.next()
		seq_file = query.value(0)
		start = int(record.value('start'))
		stop = int(record.value('stop'))
		html = get_ssr_sequence(seq_file, record.value('sequence'), start, stop, flank)
		#print ssr_seq
		dialog = BrowserDialog(self, html)
		if dialog.exec_() == QDialog.Accepted:
			pass


	def changeRowColCount(self, count):
		#self.table.setColumnWidth(0, 30)
		self.table.resizeColumnToContents(0)
		#self.table.horizontalHeader().setSortIndicator(1, Qt.AscendingOrder)
		#self.table.resizeColumnsToContents()
		labels = {'ssr':'SSRs', 'vntr':'VNTRs', 'cssr':'cSSRs', 'issr':'iSSRs', 'primer':'Primers'}
		label = labels.get(count[0], 'Row')
		self.rowCounts.setText("%s: %s" % (label, count[1]))
		self.colCounts.setText("Column: %s" % count[2])

		if count[1] > 0:
			self.exportWholeTableAct.setDisabled(False)

	def changeSelectCount(self, count):
		self.selectCounts.setText("Select: {}".format(count))
		if count > 0:
			self.exportSelectedRowsAct.setDisabled(False)
		else:
			self.exportSelectedRowsAct.setDisabled(True)

	
	def setProgress(self, percent):
		self.progressBar.setValue(percent)

	def setStatusMessage(self, msg):
		self.statusBar.showMessage(msg)

	def openAboutMessage(self):
		#system_info = "%s%s %s" % (platform.system(), platform.release(), platform.architecture()[0])
		#python_info = sys.version.split()[0]
		about_message =	"""
			<p><b>Krait - Microsatellite Identification and Primer Design</b></p>
			<p>Version v{version} Build {build}<p>
			<p>Krait is a robust and ultrafast tool that provides a user-friendly GUI for no computationally
			skilled biologists to extract perfect, imperfect and compound microsatellites and VNTRs from fasta
			formatted DNA sequences; and design primers; and perform statistical analysis.</p>
			<p><a href="https://wiki.qt.io/Qt_for_Python">PySide</a> for GUI. 
			<a href="https://github.com/lmdu/pyfastx">pyfastx</a> for parsing fasta.
			<a href="https://github.com/libnano/primer3-py">primer3-py</a> and 
			<a href="http://primer3.sourceforge.net/">primer3</a> for primer design. 
			<a href="https://github.com/hunt-genes/ncls">NCLS</a> for mapping SSRs.
			<a href="http://www.chartjs.org/">Chartjs</a> for plotting.
			</p>
			<p><b>Citation:</b><br>Du L, Zhang C, Liu Q, Zhang X, Yue B (2018) Krait: an ultrafast tool for genome-wide survey of microsatellites and primer design. Bioinformatics. 34(4):681-683.</p>
			<p><b>Contact:</b><br>dulianming@cdu.edu.cn</p>
		""".format(version=VERSION, build=BUILD)

		QMessageBox.about(self, "About Krait", about_message)

	def openDocumentation(self):
		QDesktopServices.openUrl(QUrl("https://github.com/lmdu/krait"))

	def reportIssue(self):
		QDesktopServices.openUrl(QUrl("https://github.com/lmdu/krait/issues"))

class SSRFilterInput(QLineEdit):
	def __init__(self, parent=None):
		super(SSRFilterInput, self).__init__(parent)
		self.setPlaceholderText("Filter data in table e.g. motif=AT and repeat>10")


#class SSRWebView(QWebView):
#	def __init__(self, parent=None):
#		super(SSRWebView, self).__init__(parent)
#		self.load(QUrl('http://www.baidu.com'))
#		self.page().setLinkDelegationPolicy(QWebPage.DelegateAllLinks)
#		self.linkClicked.connect(self.openUrl)
#	def openUrl(self, url):
#		QDesktopServices.openUrl(url)

class SSRTableView(QTableView):
	def __init__(self, parent=None):
		super(SSRTableView, self).__init__(parent)
		self.parent = parent
		self.verticalHeader().hide()
		self.horizontalHeader().setHighlightSections(False)
		#self.horizontalHeader().setDefaultSectionSize(150)
		self.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.setSelectionBehavior(QAbstractItemView.SelectRows)
		self.setSelectionMode(QAbstractItemView.SingleSelection)
		#self.setSelectionMode(QAbstractItemView.MultiSelection)
		self.setSortingEnabled(True)

		self.checkbox = QCheckBox(self.horizontalHeader())
		self.checkbox.setGeometry(QRect(3,5,20,20))
		self.checkbox.stateChanged.connect(self.checkboxAction)

		self.clicked.connect(self.showCellContent)


	def showCellContent(self, index):
		msg = self.model().data(index)
		self.parent.setStatusMessage(str(msg))

	def contextMenuEvent(self, event):
		self.current_row = self.rowAt(event.pos().y())
		if self.current_row == -1:
			return

		self.menu = QMenu(self)
		select_action = QAction('Select', self)
		select_action.triggered.connect(self.selectCurrentRow)
		deselect_action = QAction("Deselect", self)
		deselect_action.triggered.connect(self.deselectCurrentRow)
		select_all_action = QAction("Select All", self)
		select_all_action.setShortcut(QKeySequence(Qt.CTRL+Qt.Key_A))
		select_all_action.triggered.connect(self.selectAll)
		deselect_all_action = QAction("Deselect All", self)
		deselect_all_action.setShortcut(QKeySequence(Qt.CTRL+Qt.SHIFT+Qt.Key_A))
		deselect_all_action.triggered.connect(self.deselectAll)

		#delete_action = QAction("Delete All", self)

		detail_action = QAction("View Sequence", self)
		detail_action.triggered.connect(self.viewDetail)
		
		self.menu.addAction(select_action)
		self.menu.addAction(deselect_action)
		self.menu.addAction(select_all_action)
		self.menu.addAction(deselect_all_action)
		#self.menu.addSeparator()
		#self.menu.addAction(delete_action)
		self.menu.addSeparator()
		self.menu.addAction(detail_action)
		self.menu.popup(QCursor.pos())

	def checkboxAction(self, state):
		if state > 0:
			self.selectAll()
		else:
			self.deselectAll()

	def selectCurrentRow(self):
		self.model().selectRow(self.current_row)

	def deselectCurrentRow(self):
		self.model().deselectRow(self.current_row)

	def selectAll(self): 
		self.model().selectAll()

	def deselectAll(self):
		self.model().deselectAll()

	def viewDetail(self):
		'''
		view sequence and detail information of tandem repeats
		'''
		table = self.model().table
		flank = int(self.parent.settings.value('ssr/flank', 50))
		seed_repeat = int(self.parent.settings.value('ssr/srep', 3))
		seed_minlen = int(self.parent.settings.value('ssr/slen', 8))
		max_error = int(self.parent.settings.value('ssr/error', 3))
		_id = self.model().getCellId(self.current_row)
		
		if table == 'primer':
			content = PrimerDetail(table, _id, flank).generateHtml()
		elif table == 'issr':
			content = ISSRSeqDetail(table, _id, flank, seed_repeat, seed_minlen, max_error).generateHtml()
		else:
			content = SequenceDetail(table, _id, flank).generateHtml()

		SSRDetailDialog(self.parent, "%s Sequence" % table.upper(), content)


class TableModel(QAbstractTableModel):
	row_col = Signal(tuple)
	sel_row = Signal(int)
	def __init__(self, parent=None):
		super(TableModel, self).__init__(parent)
		self.headers = []

		#store ids of selected row
		self.selected = set()

		#store ids of displayed row
		self.displayed = []

		self.total_row_counts = 0
		self.readed_row_counts = 0
		
		#connect to database
		self.db = Database()

		#["Select", "Where", "Order"]
		self.query = ['', '', '']

		#table name
		self.table = None

		#cache row (index, data)
		self.cached_row = [-1, None]

		#mapping
		self.repeat_types = {'ssr': 1, 'cssr': 2, 'issr': 3, 'vntr': 4}

	def getRowCounts(self):
		return self.total_row_counts

	def getColCounts(self):
		return len(self.headers)

	def getAllItems(self):
		return self.total_row_counts

	def selectRow(self, row):
		if row not in self.selected:
			self.beginResetModel()
			self.selected.add(self.getCellId(row))
			self.endResetModel()
			self.sel_row.emit(len(self.selected))

	def deselectRow(self, row):
		if row in self.selected:
			self.beginResetModel()
			self.selected.remove(self.getCellId(row))
			self.endResetModel()
			self.sel_row.emit(len(self.selected))

	def selectAll(self):
		self.beginResetModel()
		sql = self.query[0] % 'COUNT(1)' + ' ' + self.query[1]
		if self.query[1]:
			self.selected = set(self.db.get_column(self.sql % 'id'))
		else:
			self.selected = set(range(1, self.total_row_counts+1))

		#if self.db.get_one(sql) == self.total_row_counts:
			#self.selected = set(range(self.total_row_counts))
		#else:
			#self.selected = set(range(len(self.displayed)))

		self.endResetModel()
		self.sel_row.emit(len(self.selected))

	def deselectAll(self):
		self.beginResetModel()
		self.selected = set()
		self.endResetModel()
		self.sel_row.emit(0)

	def setTable(self, table):
		self.table = table
		self.repeat_type = self.repeat_types.get(self.table, 0)
		self.headers = self.db.get_fields(self.table) or []
		self.query = ["SELECT %s FROM {0}".format(self.table), '', '']

	def tableName(self):
		return self.table

	def columnNames(self):
		return self.headers

	def setFilter(self, condition=None):
		if condition:
			self.query[1] = "WHERE %s" % condition
			self.select()
		else:
			self.query[1] = ''
			self.select()

	def sort(self, column, order):
		if column == 0:
			return
		
		colname = self.headers[column-1]
		if order == Qt.SortOrder.DescendingOrder:
			self.query[2] = "ORDER BY %s DESC" % colname
			self.select()
		elif order == Qt.AscendingOrder:
			self.query[2] = 'ORDER BY %s' % colname
			self.select()
		else:
			self.query[2] = ''
			self.select()

	def select(self):
		self.sql = " ".join(self.query)
		self.selected = set()
		self.cached_row = [-1, None]
		self.total_row_counts = self.db.get_one(self.query[0] % "COUNT(*)" + " " + self.query[1] + " LIMIT 1")

		self.beginResetModel()
		self.displayed = self.db.get_column(self.sql % 'id' + " LIMIT 100")
		self.readed_row_counts = len(self.displayed)
		self.endResetModel()

		self.row_col.emit((self.table, self.total_row_counts, len(self.headers)))
		self.sel_row.emit(len(self.selected))

	def clear(self):
		self.beginResetModel()
		self.total_row_counts = 0
		self.readed_row_counts = 0
		self.selected = set()
		self.cached_row = [-1, None]
		self.headers = []
		self.displayed = []
		self.endResetModel()

		self.row_col.emit((self.table, self.total_row_counts, len(self.headers)))
		self.sel_row.emit(len(self.selected))
		
	def remove(self, table):
		if table == self.table:
			self.clear()

		self.db.clear(table)


	#def getSelectedRows(self):
	#	if self.total_row_counts == len(self.selected):
	#		return 'whole'
	#	else:
	#		return [str(self.displayed[i]) for i in sorted(self.selected)]

	def getCellId(self, row):
		#return self.displayed[row]
		return self.db.get_one("%s LIMIT %s,1" % (self.sql % 'id', row))

	def value(self, index):
		row = index.row()
		col = index.column() - 1
		
		if row == self.cached_row[0]:
			return self.cached_row[1][col]
		
		ID = self.displayed[row]
		self.cached_row[0] = row
		sql = self.query[0] % '*' + " WHERE id=%s LIMIT 1" % ID
		self.cached_row[1] = self.db.get_row(sql)

		return self.cached_row[1][col]

	def rowColor(self, index):
		ID = self.displayed[index.row()]

		if self.table == 'primer' or self.db.is_empty('location'):
			return QColor(255, 255, 255)

		colors = {
			1: QColor(245, 183, 177), 
			2: QColor(250, 215, 160), 
			3: QColor(169, 223, 191),
			4: QColor(174, 214, 241)
		}
		
		sql = "SELECT feature FROM location WHERE reptype=%s AND target=%s LIMIT 1"
		
		feature = self.db.get_one(sql % (self.repeat_type, ID))

		return colors.get(feature, QColor(255, 255, 255))

	def rowCount(self, parent=QModelIndex()):
		if parent.isValid():
			return 0

		return self.readed_row_counts

	def columnCount(self, parent=QModelIndex()):
		if parent.isValid():
			return 0

		if len(self.headers) > 0:
			return len(self.headers) + 1

		return len(self.headers)

	def data(self, index, role=Qt.DisplayRole):
		if not index.isValid():
			return None

		if not (0 <= index.row() < self.rowCount()):
			return None

		elif role == Qt.DisplayRole:
			if index.column() > 0:
				return self.value(index)
			else:
				return None

		elif role == Qt.CheckStateRole:
			if index.column() == 0:
				if self.displayed[index.row()] in self.selected:
					return Qt.Checked
				else:
					return Qt.Unchecked

		#elif role == Qt.BackgroundRole:
		#	if index.row() in self.selected:
		#		return QColor(215, 242, 222)
		elif role == Qt.BackgroundColorRole:
			return self.rowColor(index)
			#return QColor(215, 242, 222)

		return None


	def headerData(self, section, orientation, role=Qt.DisplayRole):
		#if role == Qt.SizeHintRole:
		#	if section == 0:
		#		return QSize(20, -1)

		if role != Qt.DisplayRole:
			return None

		if orientation == Qt.Horizontal:
			if section == 0:
				return None
			else:
				return self.headers[section-1]
		
		elif orientation == Qt.Vertical:
			#return self.value(section, 1)
			if role == Qt.CheckStateRole:
				if index.column() == 0:
					if self.displayed[index.row()] in self.selected:
						return Qt.ItemIsUserCheckable
					else:
						return Qt.ItemIsUserCheckable

		return None

	def setHeaderData(self, section, orientation, value, role=Qt.DisplayRole):
		if role == Qt.CheckStateRole:
			if orientation == Qt.Vertical:
				return 

	def setData(self, index, value, role):
		if not index.isValid():
			return False

		if index.column() != 0:
			return False

		_id = self.getCellId(index.row())

		if role == Qt.CheckStateRole:
			if value == Qt.Checked:
				self.selected.add(_id)
			else:
				if _id in self.selected:
					self.selected.remove(_id)
			
			self.dataChanged.emit(index, index)
			self.sel_row.emit(len(self.selected))
			return True

		return False

	def flags(self, index):
		if not index.isValid():
			#return QAbstractItemModel.flags(index)
			return Qt.ItemIsSelectable

		flag = Qt.ItemIsEnabled | Qt.ItemIsSelectable

		if index.column() == 0:
			return flag | Qt.ItemIsUserCheckable

		return flag

	def canFetchMore(self, parent):
		return not parent.isValid() and (self.readed_row_counts < self.total_row_counts)

	def fetchMore(self, parent):
		if parent.isValid():
			return
		remainder = self.total_row_counts - self.readed_row_counts
		fetch_row = min(100, remainder)
		sql = self.sql % 'id' + " LIMIT %s,%s" % (self.readed_row_counts-1, fetch_row)
		IDs = self.db.get_column(sql)
		self.beginInsertRows(QModelIndex(), self.readed_row_counts, self.readed_row_counts+fetch_row-1)
		self.displayed.extend(IDs)
		self.readed_row_counts += fetch_row
		self.endInsertRows()

class AnnotationDialog(QDialog):
	def __init__(self, parent=None, gff_file='', rm_file=''):
		super(AnnotationDialog, self).__init__(parent)
		self.resize(QSize(450, 80))
		self.setWindowTitle(self.tr("Provide annotation file"))
		buttonBox = QDialogButtonBox(QDialogButtonBox.Ok)
		buttonBox.accepted.connect(self.accept)
		
		self.gff_label = QLabel(self.tr("Select GTF, GFF or gz compressed genome annotation file"), self)
		self.gff_input = QLineEdit(self)
		self.gff_input.setText(gff_file)
		self.gff_input.setReadOnly(True)
		self.gff_btn = QPushButton(self.tr("Select"), self)
		self.gff_btn.clicked.connect(self.selectAnnotationFile)

		annotLayout = QGridLayout()
		annotLayout.setColumnStretch(0, 1)
		annotLayout.addWidget(self.gff_label, 0, 0)
		annotLayout.addWidget(self.gff_input, 1, 0)
		annotLayout.addWidget(self.gff_btn, 1, 1)
		annotLayout.addWidget(buttonBox, 2, 1)

		self.setLayout(annotLayout)

	def selectAnnotationFile(self):
		filters = "GTF (*.gtf *.gtf.gz);;GFF (*.gff *.gff.gz);;ALL (*.*)"
		annot_file, _ = QFileDialog.getOpenFileName(self, "GFF or GTF annotation file", filter=filters)
		if annot_file:
			self.gff_input.setText(annot_file)

	def get(self):
		return self.gff_input.text()


class PreferenceDialog(QDialog):
	def __init__(self, parent=None, settings=None):
		super(PreferenceDialog, self).__init__(parent)
		self.settings = settings
		self.setWindowTitle(self.tr("Preferences"))
		#self.setMinimumWidth(500)

		self.general_tab = GeneralTab(self.settings)
		self.primer_tab = PrimerTab(self.settings)

		self.tabWidget = QTabWidget()
		self.tabWidget.addTab(self.general_tab, 'SSR search')
		self.tabWidget.addTab(self.primer_tab, 'Primer design')

		buttonBox = QDialogButtonBox(QDialogButtonBox.RestoreDefaults | QDialogButtonBox.Save | QDialogButtonBox.Cancel)
		buttonBox.accepted.connect(self.accept)
		buttonBox.rejected.connect(self.reject)
		buttonBox.button(QDialogButtonBox.RestoreDefaults).clicked.connect(self.resetSettings)

		spacerItem = QSpacerItem(10, 20, QSizePolicy.Minimum, QSizePolicy.Expanding)

		mainLayout = QVBoxLayout()
		mainLayout.addWidget(self.tabWidget)
		mainLayout.addItem(spacerItem)
		mainLayout.addWidget(buttonBox)

		self.setLayout(mainLayout)

	def saveSettings(self):
		self.general_tab.saveSettings()
		self.primer_tab.saveSettings()

	def resetSettings(self):
		self.settings.clear()
		self.general_tab.getSettings()
		self.primer_tab.getSettings()
		self.saveSettings()

	def gotoPrimer(self):
		self.tabWidget.setCurrentIndex(1)


class GeneralTab(QWidget):
	def __init__(self, settings, parent=None):
		super(GeneralTab, self).__init__(parent)
		self.settings = settings

		repeatsGroup = QGroupBox(self.tr("Minimum repeats for each SSR type"))
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
		repeatLayout.setColumnStretch(5, 1)
		repeatLayout.addWidget(monoLabel, 0, 0)
		repeatLayout.addWidget(self.monoValue, 0, 1)
		repeatLayout.addWidget(diLabel, 0, 2)
		repeatLayout.addWidget(self.diValue, 0, 3)
		repeatLayout.addWidget(triLabel, 0, 4)
		repeatLayout.addWidget(self.triValue, 0, 5)
		repeatLayout.addWidget(tetraLabel, 1, 0)
		repeatLayout.addWidget(self.tetraValue, 1, 1)
		repeatLayout.addWidget(pentaLabel, 1, 2)
		repeatLayout.addWidget(self.pentaValue, 1, 3)
		repeatLayout.addWidget(hexaLabel, 1, 4)
		repeatLayout.addWidget(self.hexaValue, 1, 5)
		repeatsGroup.setLayout(repeatLayout)

		distanceGroup = QGroupBox(self.tr("Compound microsatellite, cSSR"))
		distanceLabel = QLabel("Max distance allowed between two SSRs (dMAX)")
		self.distanceValue = QSpinBox()
		self.distanceValue.setSuffix(' bp')
		distanceLayout = QHBoxLayout()
		distanceLayout.addWidget(distanceLabel)
		distanceLayout.addWidget(self.distanceValue, 1)
		distanceGroup.setLayout(distanceLayout)

		satelliteGroup = QGroupBox(self.tr("Minisatellite and Macrosatellite, VNTRs"))
		min_tandem_label = QLabel("Min motif length")
		self.min_tandem_motif = QSpinBox()
		self.min_tandem_motif.setMinimum(7)

		max_tandem_label = QLabel("Max motif length")
		self.max_tandem_motif = QSpinBox()
		self.max_tandem_motif.setMinimum(7)

		repeat_tandem_label = QLabel("Min repeats")
		self.min_tandem_repeat = QSpinBox()
		self.min_tandem_repeat.setMinimum(2)

		satelliteLayout = QGridLayout()
		satelliteLayout.addWidget(min_tandem_label, 0, 0)
		satelliteLayout.addWidget(self.min_tandem_motif, 0, 1)
		satelliteLayout.addWidget(max_tandem_label, 0, 2)
		satelliteLayout.addWidget(self.max_tandem_motif, 0, 3)
		satelliteLayout.addWidget(repeat_tandem_label, 0, 4)
		satelliteLayout.addWidget(self.min_tandem_repeat, 0, 5)
		satelliteGroup.setLayout(satelliteLayout)

		issrGroup = QGroupBox(self.tr("Imperfect microsatellite, iSSR"))
		seed_mrep_label = QLabel("Min seed repeats")
		self.seed_min_repeat = QSpinBox()
		self.seed_min_repeat.setMinimum(1)
		seed_mlen_label = QLabel("Min seed length")
		self.seed_min_length = QSpinBox()
		max_error_label = QLabel("Max continuous edits")
		self.max_error = QSpinBox()
		mis_penalty_label = QLabel("Mismatch penalty")
		self.mis_penalty = QSpinBox()
		gap_penalty_label = QLabel("Gap penalty")
		self.gap_penalty = QSpinBox()
		min_score_label = QLabel("Min required score")
		self.min_score = QSpinBox()
		issrLayout = QGridLayout()
		issrLayout.addWidget(seed_mrep_label, 0, 0)
		issrLayout.addWidget(self.seed_min_repeat, 0, 1)
		issrLayout.addWidget(seed_mlen_label, 0, 2)
		issrLayout.addWidget(self.seed_min_length, 0, 3)
		issrLayout.addWidget(max_error_label, 0, 4)
		issrLayout.addWidget(self.max_error, 0, 5)
		issrLayout.addWidget(mis_penalty_label, 1, 0)
		issrLayout.addWidget(self.mis_penalty, 1, 1)
		issrLayout.addWidget(gap_penalty_label, 1, 2)
		issrLayout.addWidget(self.gap_penalty, 1, 3)
		issrLayout.addWidget(min_score_label, 1, 4)
		issrLayout.addWidget(self.min_score, 1, 5)
		issrGroup.setLayout(issrLayout)

		level_group = QGroupBox(self.tr("Motif standardization level"))
		level_label = QLabel(self.tr("Standard level"))
		self.level_select = QComboBox()
		#self.level_select.currentIndexChanged.connect(self.showStandardLevelDetail)
		#self.level_detail = QLabel()
		standard_level = [
			"Level 0  No standard",
			"Level 1  Similar motifs",
			"Level 2  Reverse complementary motifs + Level 1",
			"Level 3  complementary motifs + Level 2",
			"Level 4  Reverse motifs + Level 3 (not recommend)"
		]
		self.level_select.addItems(standard_level)
		level_layout = QHBoxLayout()
		#level_layout.setColumnStretch(1, 1)
		level_layout.addWidget(level_label)
		level_layout.addWidget(self.level_select, 1)
		#level_layout.addWidget(self.level_detail, 1, 1)
		level_group.setLayout(level_layout)

		flankGroup = QGroupBox(self.tr("Flanking sequence for exporting and primer design"))
		flankLabel = QLabel("Flanking sequence length")
		self.flankValue = QSpinBox()
		self.flankValue.setSuffix(' bp')
		self.flankValue.setMaximum(1000)
		flankLayout = QHBoxLayout()
		flankLayout.addWidget(flankLabel)
		flankLayout.addWidget(self.flankValue, 1)
		flankGroup.setLayout(flankLayout)
		
		mainLayout = QVBoxLayout()
		mainLayout.addWidget(repeatsGroup)
		mainLayout.addWidget(distanceGroup)
		mainLayout.addWidget(satelliteGroup)
		mainLayout.addWidget(issrGroup)
		mainLayout.addWidget(level_group)
		mainLayout.addWidget(flankGroup)
		self.setLayout(mainLayout)
		self.getSettings()

	def getSettings(self):
		self.monoValue.setValue(int(self.settings.value('ssr/mono', 12)))
		self.diValue.setValue(int(self.settings.value('ssr/di', 7)))
		self.triValue.setValue(int(self.settings.value('ssr/tri', 5)))
		self.tetraValue.setValue(int(self.settings.value('ssr/tetra', 4)))
		self.pentaValue.setValue(int(self.settings.value('ssr/penta', 4)))
		self.hexaValue.setValue(int(self.settings.value('ssr/hexa', 4)))
		self.distanceValue.setValue(int(self.settings.value('ssr/dmax', 10)))
		self.flankValue.setValue(int(self.settings.value('ssr/flank', 100)))
		self.min_tandem_motif.setValue(int(self.settings.value('ssr/vmin', 7)))
		self.max_tandem_motif.setValue(int(self.settings.value('ssr/vmax', 30)))
		self.min_tandem_repeat.setValue(int(self.settings.value('ssr/vrep', 2)))
		self.seed_min_repeat.setValue(int(self.settings.value('ssr/srep', 3)))
		self.seed_min_length.setValue(int(self.settings.value('ssr/slen', 8)))
		self.max_error.setValue(int(self.settings.value('ssr/error', 3)))
		self.min_score.setValue(int(self.settings.value('ssr/score', 10)))
		self.mis_penalty.setValue(int(self.settings.value('ssr/mismatch', 1)))
		self.gap_penalty.setValue(int(self.settings.value('ssr/gap', 2)))
		self.level_select.setCurrentIndex(int(self.settings.value('ssr/level', 3)))


	def saveSettings(self):
		self.settings.setValue('ssr/mono', self.monoValue.value())
		self.settings.setValue('ssr/di', self.diValue.value())
		self.settings.setValue('ssr/tri', self.triValue.value())
		self.settings.setValue('ssr/tetra', self.tetraValue.value())
		self.settings.setValue('ssr/penta', self.pentaValue.value())
		self.settings.setValue('ssr/hexa', self.hexaValue.value())
		self.settings.setValue('ssr/dmax', self.distanceValue.value())
		self.settings.setValue('ssr/flank', self.flankValue.value())
		self.settings.setValue('ssr/vmin', self.min_tandem_motif.value())
		self.settings.setValue('ssr/vmax', self.max_tandem_motif.value())
		self.settings.setValue('ssr/vrep', self.min_tandem_repeat.value())
		self.settings.setValue('ssr/level', self.level_select.currentIndex())
		self.settings.setValue('ssr/srep', self.seed_min_repeat.value())
		self.settings.setValue('ssr/slen', self.seed_min_length.value())
		self.settings.setValue('ssr/error', self.max_error.value())
		self.settings.setValue('ssr/score', self.min_score.value())
		self.settings.setValue('ssr/mismatch', self.mis_penalty.value())
		self.settings.setValue('ssr/gap', self.gap_penalty.value())

	#def showStandardLevelDetail(self, idx):
	#	if idx == 0:
	#		detail = 'No standard'

	#	elif idx == 1:
	#		detail = 'standard similar motifs'

	#	elif idx == 2:
	#		detail = 'Level 1 + reverse complementary motifs'

	#	elif idx == 3:
	#		detail = 'Level 2 + complementary motifs'

	#	elif idx == 4:
	#		detail = 'Level 3 + reverse motifs (not recommend)'
		
	#	self.level_detail.setText(detail)


class PrimerTab(QWidget):
	def __init__(self, settings, parent=None):
		super(PrimerTab, self).__init__(parent)
		self.settings = settings
		
		self.product_size = QLineEdit()
		self.primer_num_return = QSpinBox()
		self.repeat_library = QComboBox()
		self.repeat_library.addItems(['None', 'Human', 'Rodent', 'Rodent and Simple', 'Drosophila'])

		product_size_group = QGroupBox(self.tr('General Settings'))
		product_size_layout = QGridLayout()

		product_size_layout.addWidget(PrimerTagLabel('PRIMER_PRODUCT_SIZE_RANGE'), 0, 0)
		product_size_layout.addWidget(self.product_size, 0, 1, 1, 3)
		product_size_layout.addWidget(PrimerTagLabel('PRIMER_MISPRIMING_LIBRARY'), 1, 0)
		product_size_layout.addWidget(self.repeat_library, 1, 1)
		product_size_layout.addWidget(PrimerTagLabel("PRIMER_NUM_RETURN"), 1, 2)
		product_size_layout.addWidget(self.primer_num_return, 1, 3)
		
		product_size_group.setLayout(product_size_layout)

		primer_size_group = QGroupBox(self.tr("Primer Size and GC content"))
		primer_size_layout = QGridLayout()
		self.primer_size_min = QSpinBox()
		self.primer_size_opt = QSpinBox()
		self.primer_size_max = QSpinBox()
		self.primer_gc_min = QSpinBox()
		self.primer_gc_max = QSpinBox()
		self.primer_gc_clamp = QSpinBox()
		primer_size_layout.addWidget(PrimerTagLabel("PRIMER_MIN_SIZE"), 0, 0)
		primer_size_layout.addWidget(self.primer_size_min, 0, 1)
		primer_size_layout.addWidget(PrimerTagLabel("PRIMER_OPT_SIZE"), 0, 2)
		primer_size_layout.addWidget(self.primer_size_opt, 0, 3)
		primer_size_layout.addWidget(PrimerTagLabel("PRIMER_MAX_SIZE"), 0, 4)
		primer_size_layout.addWidget(self.primer_size_max, 0, 5)
		primer_size_layout.addWidget(PrimerTagLabel("PRIMER_MIN_GC"), 1, 0)
		primer_size_layout.addWidget(self.primer_gc_min, 1, 1)
		primer_size_layout.addWidget(PrimerTagLabel("PRIMER_MAX_GC"), 1, 2)
		primer_size_layout.addWidget(self.primer_gc_max, 1, 3)
		primer_size_layout.addWidget(PrimerTagLabel("PRIMER_GC_CLAMP"), 1, 4)
		primer_size_layout.addWidget(self.primer_gc_clamp, 1, 5)
		primer_size_group.setLayout(primer_size_layout)

		primer_tm_group = QGroupBox(self.tr("Primer Melting Temperature"))
		primer_tm_layout = QGridLayout()
		self.primer_tm_min = QSpinBox()
		self.primer_tm_opt = QSpinBox()
		self.primer_tm_max = QSpinBox()
		self.primer_tm_pair = QSpinBox()
		primer_tm_layout.addWidget(PrimerTagLabel("PRIMER_MIN_TM"), 0, 0)
		primer_tm_layout.addWidget(self.primer_tm_min, 0, 1)
		primer_tm_layout.addWidget(PrimerTagLabel("PRIMER_OPT_TM"), 0, 2)
		primer_tm_layout.addWidget(self.primer_tm_opt, 0, 3)
		primer_tm_layout.addWidget(PrimerTagLabel("PRIMER_MAX_TM"), 0, 4)
		primer_tm_layout.addWidget(self.primer_tm_max, 0, 5)
		primer_tm_group.setLayout(primer_tm_layout)

		primer_bind_group = QGroupBox(self.tr("Self-binding (primer-dimer and hairpins)"))
		primer_bind_layout = QGridLayout()
		self.primer_max_self_any = QSpinBox()
		self.primer_pair_max_compl_any = QSpinBox()
		self.primer_max_self_end = QSpinBox()
		self.primer_pair_max_compl_end = QSpinBox()
		self.primer_max_hairpin = QSpinBox()
		primer_bind_layout.addWidget(PrimerTagLabel("PRIMER_MAX_SELF_ANY_TH"), 0, 0)
		primer_bind_layout.addWidget(self.primer_max_self_any, 0, 1)
		primer_bind_layout.addWidget(PrimerTagLabel("PRIMER_PAIR_MAX_COMPL_ANY_TH"), 0, 2)
		primer_bind_layout.addWidget(self.primer_pair_max_compl_any, 0, 3)
		primer_bind_layout.addWidget(PrimerTagLabel("PRIMER_MAX_SELF_END_TH"), 1, 0)
		primer_bind_layout.addWidget(self.primer_max_self_end, 1, 1)
		primer_bind_layout.addWidget(PrimerTagLabel("PRIMER_PAIR_MAX_COMPL_END_TH"), 1, 2)
		primer_bind_layout.addWidget(self.primer_pair_max_compl_end, 1, 3)
		primer_bind_layout.addWidget(PrimerTagLabel("PRIMER_MAX_HAIRPIN_TH"), 2, 0)
		primer_bind_layout.addWidget(self.primer_max_hairpin, 2, 1)
		primer_bind_group.setLayout(primer_bind_layout)

		primer_other_group = QGroupBox(self.tr("PolyX and Other"))
		primer_other_layout = QGridLayout()
		self.primer_max_end_stability = QSpinBox()
		self.primer_max_ns_accepted = QSpinBox()
		self.primer_max_poly_x = QSpinBox()
		primer_other_layout.addWidget(PrimerTagLabel("PRIMER_MAX_END_STABILITY"), 0, 0)
		primer_other_layout.addWidget(self.primer_max_end_stability, 0, 1)
		primer_other_layout.addWidget(PrimerTagLabel("PRIMER_MAX_POLY_X"), 0, 2)
		primer_other_layout.addWidget(self.primer_max_poly_x, 0, 3)
		primer_other_layout.addWidget(PrimerTagLabel("PRIMER_MAX_NS_ACCEPTED"), 1, 0)
		primer_other_layout.addWidget(self.primer_max_ns_accepted, 1, 1)
		primer_other_layout.addWidget(PrimerTagLabel("PRIMER_PAIR_MAX_DIFF_TM"), 1, 2)
		primer_other_layout.addWidget(self.primer_tm_pair, 1, 3)
		primer_other_group.setLayout(primer_other_layout)

		mainLayout = QVBoxLayout()
		mainLayout.addWidget(product_size_group)
		mainLayout.addWidget(primer_size_group)
		mainLayout.addWidget(primer_tm_group)
		mainLayout.addWidget(primer_bind_group)
		mainLayout.addWidget(primer_other_group)

		self.setLayout(mainLayout)
		self.getSettings()

	def getSettings(self):
		self.product_size.setText(self.settings.value('primer/PRIMER_PRODUCT_SIZE_RANGE', '100-300'))
		self.primer_size_min.setValue(int(self.settings.value('primer/PRIMER_MIN_SIZE', 18)))
		self.primer_size_opt.setValue(int(self.settings.value('primer/PRIMER_OPT_SIZE', 20)))
		self.primer_size_max.setValue(int(self.settings.value('primer/PRIMER_MAX_SIZE', 27)))
		self.primer_tm_min.setValue(int(self.settings.value('primer/PRIMER_MIN_TM', 58)))
		self.primer_tm_opt.setValue(int(self.settings.value('primer/PRIMER_OPT_TM', 60)))
		self.primer_tm_max.setValue(int(self.settings.value('primer/PRIMER_MAX_TM', 65)))
		self.primer_gc_min.setValue(int(self.settings.value('primer/PRIMER_MIN_GC', 30)))
		self.primer_gc_max.setValue(int(self.settings.value('primer/PRIMER_MAX_GC', 80)))
		self.primer_gc_clamp.setValue(int(self.settings.value('primer/PRIMER_GC_CLAMP', 2)))
		self.primer_tm_pair.setValue(int(self.settings.value('primer/PRIMER_PAIR_MAX_DIFF_TM', 2)))
		self.primer_max_self_any.setValue(int(self.settings.value('primer/PRIMER_MAX_SELF_ANY_TH', 47)))
		self.primer_pair_max_compl_any.setValue(int(self.settings.value('primer/PRIMER_PAIR_MAX_COMPL_ANY_TH', 47)))
		self.primer_max_self_end.setValue(int(self.settings.value('primer/PRIMER_MAX_SELF_END_TH', 47)))
		self.primer_pair_max_compl_end.setValue(int(self.settings.value('primer/PRIMER_PAIR_MAX_COMPL_END_TH', 47)))
		self.primer_max_hairpin.setValue(int(self.settings.value('primer/PRIMER_MAX_HAIRPIN_TH', 47)))
		self.primer_max_end_stability.setValue(int(self.settings.value('primer/PRIMER_MAX_END_STABILITY', 100)))
		self.primer_max_ns_accepted.setValue(int(self.settings.value('primer/PRIMER_MAX_POLY_X', 5)))
		self.primer_max_poly_x.setValue(int(self.settings.value('primer/PRIMER_MAX_NS_ACCEPTED', 0)))
		self.primer_num_return.setValue(int(self.settings.value('primer/PRIMER_NUM_RETURN', 1)))
		self.repeat_library.setCurrentIndex(int(self.settings.value('primer/PRIMER_MISPRIMING_LIBRARY', 0)))


	def saveSettings(self):
		self.settings.setValue('primer/PRIMER_PRODUCT_SIZE_RANGE', self.product_size.text())
		self.settings.setValue('primer/PRIMER_MIN_SIZE', self.primer_size_min.value())
		self.settings.setValue('primer/PRIMER_OPT_SIZE', self.primer_size_opt.value())
		self.settings.setValue('primer/PRIMER_MAX_SIZE', self.primer_size_max.value())
		self.settings.setValue('primer/PRIMER_MIN_TM', self.primer_tm_min.value())
		self.settings.setValue('primer/PRIMER_OPT_TM', self.primer_tm_opt.value())
		self.settings.setValue('primer/PRIMER_MAX_TM', self.primer_tm_max.value())
		self.settings.setValue('primer/PRIMER_MIN_GC', self.primer_gc_min.value())
		self.settings.setValue('primer/PRIMER_MAX_GC', self.primer_gc_max.value())
		self.settings.setValue('primer/PRIMER_GC_CLAMP', self.primer_gc_clamp.value())
		self.settings.setValue('primer/PRIMER_PAIR_MAX_DIFF_TM', self.primer_tm_pair.value())
		self.settings.setValue('primer/PRIMER_MAX_SELF_ANY_TH', self.primer_max_self_any.value())
		self.settings.setValue('primer/PRIMER_PAIR_MAX_COMPL_ANY_TH', self.primer_pair_max_compl_any.value())
		self.settings.setValue('primer/PRIMER_MAX_SELF_END_TH', self.primer_max_self_end.value())
		self.settings.setValue('primer/PRIMER_PAIR_MAX_COMPL_END_TH', self.primer_pair_max_compl_end.value())
		self.settings.setValue('primer/PRIMER_MAX_HAIRPIN_TH', self.primer_max_hairpin.value())
		self.settings.setValue('primer/PRIMER_MAX_END_STABILITY', self.primer_max_end_stability.value())
		self.settings.setValue('primer/PRIMER_MAX_POLY_X', self.primer_max_ns_accepted.value())
		self.settings.setValue('primer/PRIMER_MAX_NS_ACCEPTED', self.primer_max_poly_x.value())
		self.settings.setValue('primer/PRIMER_NUM_RETURN', self.primer_num_return.value())
		self.settings.setValue('primer/PRIMER_MISPRIMING_LIBRARY', self.repeat_library.currentIndex())


class PrimerTagLabel(QLabel):
	base_url = "http://primer3.sourceforge.net/primer3_manual.htm#%s"
	def __init__(self, tag, parent=None):
		super(PrimerTagLabel, self).__init__(parent)
		#self.setOpenExternalLinks(True)
		self.tag = tag
		self.setText('<a href="#%s">%s</a>' % (self.tag, self.tag))
		self.linkActivated.connect(self.openLink)

	def openLink(self, link):
		url = QUrl.fromLocalFile(self.base_url % self.tag)
		QDesktopServices.openUrl(url)


class SSRDetailDialog(QDialog):
	def __init__(self, parent=None, title=None, content=None):
		super(SSRDetailDialog, self).__init__(parent)
		self.setWindowTitle(title)
		self.viewer = QWebEngineView(self)
		#self.viewer = QWebView(self)
		self.viewer.setHtml(content, QUrl("qrc:/"))

		#buttonBox = QDialogButtonBox(QDialogButtonBox.Ok)
		#buttonBox.accepted.connect(self.accept)
		
		mainLayout = QVBoxLayout()
		#mainLayout.setSpacing(0)
		#mainLayout.setMargin(0)
		mainLayout.addWidget(self.viewer)
		#mainLayout.addWidget(buttonBox)
		self.resize(700, 500)

		self.setLayout(mainLayout)
		self.open()

class DownloadDialog(QDialog):
	def __init__(self, parent=None):
		super(DownloadDialog, self).__init__(parent)
		self.setMinimumWidth(500)
		self.setWindowTitle("Download sequence from NCBI")
		self.acc_label = QLabel(self.tr("NCBI Accession Number:"), self)
		self.acc_input = QLineEdit(self)
		self.out_label = QLabel(self.tr("Select output path:"), self)
		self.out_input = QLineEdit(self)
		self.browser_btn = QPushButton(self.tr("Browser"), self)
		self.browser_btn.clicked.connect(self.select)
		self.out_input.setReadOnly(True)

		layout = QGridLayout()
		layout.setColumnStretch(0, 1)
		layout.addWidget(self.acc_label, 0, 0)
		layout.addWidget(self.acc_input, 1, 0, 1, 2)
		layout.addWidget(self.out_label, 2, 0)
		layout.addWidget(self.out_input, 3, 0)
		layout.addWidget(self.browser_btn, 3, 1)

		buttonBox = QDialogButtonBox(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
		buttonBox.accepted.connect(self.accept)
		buttonBox.rejected.connect(self.reject)

		layout.addWidget(buttonBox, 4, 0, Qt.AlignLeft)

		self.setLayout(layout)

	def select(self):
		exp_file, _ = QFileDialog.getSaveFileName(self, filter="Fasta (*.fa);;Fasta (*.fasta)")
		if exp_file:
			self.out_input.setText(exp_file)

	def get(self):
		acc = self.acc_input.text().strip()
		out = self.out_input.text()
		return acc, out


class BrowserWidget(QWebEngineView):
#class BrowserWidget(QWebView):
	def __init__(self, parent):
		super(BrowserWidget, self).__init__(parent)
		self.parent = parent
		self.db = Database()

		#self.pageAction(QWebEnginePage.DownloadImageToDisk).triggered.connect(self.saveImage)
		#self.pageAction(QWebEnginePage.OpenImageInNewWindow).triggered.connect(self.openImage)
		#self.page().setLinkDelegationPolicy(QWebEnginePage.DelegateAllLinks)
		#self.page().linkClicked.connect(self.saveTable)

	#def openImage(self):
	#	url = self.page().mainFrame().hitTestContent(QCursor.pos()).imageUrl()
	#	QDesktopServices.openUrl(url)

	#def saveImage(self):
	#	pm = self.page().mainFrame().hitTestContent(QCursor.pos()).pixmap()
	#	filepath, _ = QFileDialog.getSaveFileName(self, "Save image", 
	#		filter="TIFF File (*.tiff);;JPEG File (*.jpg);;PNG File (*.png)"
	#	)
		
	#	if not filepath: return
	#	pm.save(filepath)

	#	self.parent.setStatusMessage('Image %s has been successfully saved' % filepath)

	#def saveTable(self, url):
	#	url = url.toString()
	#	table, name = url.split('/')[-1].split('-')
	#	stats_str = self.db.get_option('%s_statis' % table)
	#	stats_obj = json.loads(stats_str)
	#	outfile, _ = QFileDialog.getSaveFileName(self, filter="CSV (*.csv)")
	#	if not outfile: return

	#	with open(outfile, 'wb') as fh:
	#		writer = csv.writer(fh)
	#		writer.write(stats_obj[name][0])
	#		write_to_csv(writer, stats_obj[name][1:])
		
	#	self.parent.setStatusMessage("Table %s has been successfully saved" % outfile)


#class SSRTableModel(QSqlTableModel):
#	def __init__(self):
#		super(SSRTableModel, self).__init__()
#		self.setTable('ssr')
#		self.select()

	#def data(self, index, role=Qt.DisplayRole):
	#	if role == Qt.TextAlignmentRole:
	#		return Qt.AlignCenter

	#	return QSqlTableModel.data(self, index, role)