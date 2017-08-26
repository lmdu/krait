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

from PySide.QtCore import *
from PySide.QtGui import *
from PySide.QtSql import *
from PySide.QtWebKit import *

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
		self.setWindowIcon(QIcon('logo.ico'))

		#stacked widget
		self.main_widget = QStackedWidget(self)
		self.setCentralWidget(self.main_widget)
		
		self.table = SSRTableView(self)
		self.createTableModel()

		self.browser = QWebView(self)
		self.browser.linkClicked.connect(self.saveStatTableFigure)

		self.main_widget.addWidget(self.table)
		self.main_widget.addWidget(self.browser)

		#self.setCentralWidget(self.browser)

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
		self.annot_file = ''

		#statistical results
		self.statis_result = None

		#read settings
		self.readSettings()

		#read home page
		self.homepage()

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
		content = template_render('index.html')
		self.browser.setHtml(content, QUrl.fromLocalFile(CACHE_PATH))

	def createTableModel(self):
		self.model = TableModel()
		self.table.setModel(self.model)
		self.model.row_col.connect(self.changeRowColCount)
		self.model.sel_row.connect(self.changeSelectCount)

	def saveStatTableFigure(self, href):
		href = href.toString()
		media, name = href.split(':')

		if media == 'table':
			table, name = name.split('-')
			stats_str = self.db.get_option('%s_statis' % table)
			stats_obj = json.loads(stats_str)
			outfile, _ = QFileDialog.getSaveFileName(self, filter="CSV (*.csv)")
			if not outfile: return
			write_to_csv(outfile, stats_obj[name][0], stats_obj[name][1:])

		elif media == 'figure':
			src = os.path.join(CACHE_PATH, "%s.png" % name)
			dst, fmat = QFileDialog.getSaveFileName(self, dir=name, filter="PNG (*.png);;JPG (*.jpg);;TIFF (*.tif)")
			if not dst: return
			img = QImage(src)
			img.save(dst, fmat.split()[0])

		else:
			pass

	def readSettings(self):
		self.settings = QSettings("config.ini", QSettings.IniFormat)
		self.resize(self.settings.value("size", QSize(900, 600)))

	def writeSettings(self):
		self.settings.setValue("size", self.size())

	def closeEvent(self, event):
		self.writeSettings()
		if not self.db.changes():
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
		#self.selectAllAct.setShortcut(QKeySequence.SelectAll)
		self.selectAllAct.triggered.connect(self.doSelectAll)

		self.preferenceAct = QAction(self.tr("Preferences"), self)
		self.preferenceAct.setShortcut(QKeySequence.Preferences)
		self.preferenceAct.triggered.connect(self.setPreference)

		#toolbar actions
		#search perfect ssrs tool button
		self.SSRSearchAct = QAction(QIcon(":/icons/ssr.png"), self.tr("SSRs"), self)
		self.SSRSearchAct.setToolTip(self.tr("Search for Perfect SSRs"))
		self.SSRSearchAct.triggered.connect(self.searchOrShowSSR)
		self.SSRForceAct = QAction(self.tr("Search for SSRs"), self)
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
		self.CSSRSearchAct.triggered.connect(self.searchOrShowCSSR)
		self.CSSRForceAct = QAction(self.tr("Search for cSSRs"), self)
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
		self.VNTRSearchAct.triggered.connect(self.searchOrShowVNTR)
		self.VNTRForceAct = QAction(self.tr("Search for VNTRs"), self)
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
		self.ISSRSearchAct.triggered.connect(self.searchOrShowISSR)
		self.ISSRForceAct = QAction(self.tr("Search for iSSRs"), self)
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
		self.locateAct = QAction(QIcon(":/icons/locate.png"), self.tr("Locate"), self)
		self.locateAct.setToolTip(self.tr("Locate SSRs in genes"))
		self.locateAct.triggered.connect(self.locateTandem)
		self.locateToolAct = QAction(self.tr("Locate SSRs in genes"), self)
		self.locateToolAct.triggered.connect(self.locateTandem)

		self.locateSetAct = QAction(self.tr("Import Annotation File"), self)
		self.locateSetAct.triggered.connect(self.provideAnnotation)
		self.removeLocateAct = QAction(self.tr("Remove locations"), self)
		self.removeLocateAct.triggered.connect(self.removeMarker)
		
		cds_icon = QPixmap(16, 16)
		cds_icon.fill(QColor(245, 183, 177))
		self.showCDSAct = QAction(QIcon(cds_icon), self.tr("Show SSRs in CDS"), self)
		self.showCDSAct.triggered.connect(self.showCDSMarker)

		exon_icon = QPixmap(16, 16)
		exon_icon.fill(QColor(169, 223, 191))
		self.showExonAct = QAction(QIcon(exon_icon), self.tr("Show SSRs in Exon"), self)
		self.showExonAct.triggered.connect(self.showExonMarker)

		utr_icon = QPixmap(16, 16)
		utr_icon.fill(QColor(250, 215, 160))
		self.showUTRAct = QAction(QIcon(utr_icon), self.tr("Show SSRs in UTR"), self)
		self.showUTRAct.triggered.connect(self.showUTRMarker)

		intron_icon = QPixmap(16, 16)
		intron_icon.fill(QColor(174, 214, 241))
		self.showIntronAct = QAction(QIcon(intron_icon), self.tr("Show SSRs in Intron"), self)
		self.showIntronAct.triggered.connect(self.showIntronMarker)


		#design primer
		self.primerDesignAct = QAction(QIcon(":/icons/primer.png"), self.tr("Primer"), self)
		self.primerDesignAct.setToolTip(self.tr("Design primers"))
		self.primerDesignAct.triggered.connect(self.designOrShowPrimer)
		self.primerForceAct = QAction(self.tr("Design Primers"), self)
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
		self.statisticsForceAct = QAction(self.tr("Statistical Analysis"), self)
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
		self.fileMenu.addSeparator()
		self.fileMenu.addAction(self.exportTableAct)
		self.fileMenu.addAction(self.exportFastaAct)
		self.fileMenu.addAction(self.exportGFFAct)
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
		

	def openProject(self):
		dbfile, _ = QFileDialog.getOpenFileName(self, filter="Database (*.db)")
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

	def saveProject(self):
		if self.opened_project is None:
			dbfile, _ = QFileDialog.getSaveFileName(self, filter="Database (*.db)")
			if not dbfile:
				return
			self.opened_project = dbfile

		os.remove(self.opened_project)

		self.db.save(self.opened_project)

	def saveProjectAs(self):
		dbfile, _ = QFileDialog.getSaveFileName(self, filter="Database (*.db)")
		if not dbfile:
			return

		self.db.save(dbfile)

	def closeProject(self):
		self.db.drop_tables()
		self.createTableModel()
	
	def importFasta(self):
		'''
		Import a fasta file from a directory
		'''
		fasta, _ = QFileDialog.getOpenFileName(self, 
			filter="Fasta (*.fa *.fna *.fas *.fasta *.fna.gz *.fa.gz *.fasta.gz);;All files (*.*)")
		if not fasta: return
		self.db.get_cursor().execute('INSERT INTO fasta VALUES (?,?)', (None, fasta))
		#self.fasta_table.insert(Data(fid=None, path=fasta))
		self.setStatusMessage("Import fasta %s" % fasta)

	def importFastas(self):
		'''
		import all fasta files from a directory
		'''
		directory = QFileDialog.getExistingDirectory(self)
		if not directory: return
		folder = QDir(directory)
		count = 0
		for fasta in  folder.entryList(QDir.Files):
			self.db.get_cursor().execute("INSERT INTO fasta VALUES (?,?)", (None, folder.absoluteFilePath(fasta)))
			count += 1
		self.setStatusMessage("Import %s fastas in %s" % (count, directory))

	def downloadFasta(self):
		'''
		download fasta file from NCBI database
		'''
		dialog = DownloadDialog(self)
		status = dialog.exec_()
		if not status: return
		acc, out = dialog.get()
		if not acc or not out: return
		self.progressBar.setMaximum(0)
		worker = EutilWorker(acc, out)
		self.executeTask(worker, lambda: self.progressBar.setMaximum(100))

	def exportTableRows(self):
		selected = self.model.getSelectedRows()
		if not selected:
			return QMessageBox.warning(self, 'Warning', "Please select rows in table to export.")

		table = self.model.tableName()
		headers = self.model.columnNames()

		exp_file, _ = QFileDialog.getSaveFileName(self, filter="CSV (*.csv);;GFF (*.gff);;GTF (*.gtf);;Tabular text (*.txt)")
		if not exp_file: return

		if len(selected) == self.db.get_one("SELECT COUNT(1) FROM %s" % table):
			sql = "SELECT * FROM %s" % table
		else:
			sql = "SELECT * FROM %s WHERE id IN (%s)" % (table, ",".join(map(str, selected.values())))

		rows = self.db.query(sql)

		if exp_file.endswith('.csv'):
			write_to_csv(exp_file, headers, rows)
		else:
			write_to_tab(exp_file, headers, rows)

		QMessageBox.information(self, "Information", "Successfully exported to %s" % exp_file)

	def exportTableGFF(self):
		selected = self.model.getSelectedRows()
		if not selected:
			return QMessageBox.warning(self, 'Warning', "Please select rows in table to export.")

		table = self.model.tableName()
		headers = self.model.columnNames()

		exp_file, _ = QFileDialog.getSaveFileName(self, filter="GFF3 (*.gff)")
		if not exp_file: return

		if len(selected) == self.db.get_one("SELECT COUNT(1) FROM %s" % table):
			sql = "SELECT * FROM %s" % table
		else:
			sql = "SELECT * FROM %s WHERE id IN (%s)" % (table, ",".join(map(str, selected.values())))
		
		rows = self.db.query(sql)

		write_to_gff(exp_file, table.upper(), rows)
		
		QMessageBox.information(self, "Information", "Successfully exported to %s" % exp_file)

	def exportTableFastas(self):
		selected = self.model.getSelectedRows()
		if not selected:
			return QMessageBox.warning(self, 'Warning', "Please select rows in table to export.")

		table = self.model.tableName()

		#if table not in ('ssr', 'issr', 'cssr', 'vntr'):
		#	return QMessageBox.warning(self, 'Warning', "Your selected rows are not SSRs")

		exp_file, _ = QFileDialog.getSaveFileName(self, filter="Fasta (*.fa);;Fasta (*.fasta)")
		if not exp_file: return

		flank = int(self.settings.value('ssr/flank'))
		worker = ExportFastaWorker(table, selected, flank, exp_file)
		self.executeTask(worker, lambda: QMessageBox.information(self, "Information", "Successfully exported to %s" % exp_file))

	def exportStatisResult(self):
		pdfname, _ = QFileDialog.getSaveFileName(self, filter="PDF (*.pdf)")
		if not file: return
		printer = QPrinter(QPrinter.HighResolution)
		printer.setPageSize(QPrinter.A4)
		printer.setColorMode(QPrinter.Color)
		printer.setOutputFormat(QPrinter.PdfFormat)
		printer.setOutputFileName(pdfname)
		self.browser.print_(printer)


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
		self.worker.finished.connect(finish_callback)
		self.work_thread = QThread()
		self.work_thread.started.connect(self.worker.process)
		self.worker.finished.connect(self.work_thread.quit)
		self.worker.finished.connect(self.worker.deleteLater)
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
			int(self.settings.value('ssr/mono')),
			int(self.settings.value('ssr/di')),
			int(self.settings.value('ssr/tri')), 
			int(self.settings.value('ssr/tetra')),
			int(self.settings.value('ssr/penta')),
			int(self.settings.value('ssr/hexa'))
		]
		level = int(self.settings.value('ssr/level'))
		worker = SSRWorker(fastas, rules, level)
		self.executeTask(worker, self.showSSR)
		#proc = SSRTask(fastas, rules, level)

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
		self.db.clear('ssr')

	
	#handle compound SSRs search
	def searchCSSR(self):
		if self.db.is_empty('ssr'):
			QMessageBox.warning(self, "Warning", "Please search perfect SSRs first, before search compound SSRs.")
			return

		self.removeCSSR()

		dmax = int(self.settings.value('ssr/dmax'))
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
		self.db.clear('cssr')


	#handle VNTRs search
	def searchVNTR(self):
		fastas = self.getInputFastas()
		if not fastas:
			return

		self.removeVNTR()

		min_motif = int(self.settings.value('ssr/vmin'))
		max_motif = int(self.settings.value('ssr/vmax'))
		min_repeat = int(self.settings.value('ssr/vrep'))
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
		self.db.clear('vntr')

	
	#handle imperfect SSRs search
	def searchISSR(self):
		fastas = self.getInputFastas()
		if not fastas:
			return

		self.removeISSR()

		seed_repeat = int(self.settings.value('ssr/srep'))
		seed_length = int(self.settings.value('ssr/slen'))
		max_eidts = int(self.settings.value('ssr/error'))
		mis_penalty = int(self.settings.value('ssr/mismatch'))
		gap_penalty = int(self.settings.value('ssr/gap'))
		score = int(self.settings.value('ssr/score'))
		level = int(self.settings.value('ssr/level'))
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
		self.db.clear('issr')


	def getPrimerSettings(self):
		p3_settings = dict(
			PRIMER_TASK = 'generic',
			PRIMER_PICK_LEFT_PRIMER = 1,
			PRIMER_PICK_INTERNAL_OLIGO = 0,
			PRIMER_PICK_RIGHT_PRIMER = 1
		)
		self.settings.beginGroup("primer")
		keys = self.settings.childKeys()
		for key in keys:
			if key == 'PRIMER_PRODUCT_SIZE_RANGE':
				sgs = self.settings.value(key)
				p3_settings[key] = [map(int, sg.split('-')) for sg in sgs.split()]
			else:
				p3_settings[key] = int(self.settings.value(key))
		self.settings.endGroup()
		return p3_settings

	#design primers for ssrs
	def designPrimer(self):
		rows = self.model.getSelectedRows()
		if not rows:
			return QMessageBox.warning(self, "Warning", "Please select SSR, iSSRs, cSSRs or VNTRs")	
		
		table = self.model.table
		flank = int(self.settings.value('ssr/flank'))
		primer3_settings = self.getPrimerSettings()
		worker = PrimerWorker(table, rows, flank, primer3_settings)

		self.removePrimer()

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
		self.db.clear('primer')

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


	def locateTandem(self):
		if not self.annot_file:
			return QMessageBox.warning(self, "Warning", "Please provide gtf or gff well formated annotation file")

		table = self.model.table
		worker = LocateWorker(table, self.annot_file)
		self.executeTask(worker, lambda : 1)

	def removeMarker(self):
		self.db.clear('location')

	def showMarker(self, marker):
		table = self.model.table
		sql = "SELECT target FROM location WHERE category='%s' AND feature='%s'" % (table, marker)
		data = self.db.get_column(sql)
		if not data:
			return QMessageBox.warning(self, "Warning", "No %ss located in %s region" % (table.upper(), marker))
		self.model.setFilter('id IN (%s)' % ",".join(map(str, data)))

	def showCDSMarker(self):
		self.showMarker('CDS')

	def showExonMarker(self):
		self.showMarker('EXON')

	def showUTRMarker(self):
		table = self.model.table
		sql = "SELECT target FROM location WHERE category='%s' AND feature IN ('3UTR', '5UTR')" % table
		data = self.db.get_column(sql)
		if not data:
			return QMessageBox.warning(self, "Warning", "No %ss located in UTR region" % table.upper())
		self.model.setFilter('id IN (%s)' % ",".join(map(str, data)))

	def showIntronMarker(self):
		self.showMarker('INTRON')

	def estimateBestMaxDistance(self):
		pass

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
		self.executeTask(worker, self.showStatistics)

	def doOrShowStatistics(self):
		if self.statis_result:
			self.swichMainWidget('browser')
			self.browser.setHtml(self.statis_result, QUrl.fromLocalFile(CACHE_PATH))
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
		self.browser.setHtml(self.statis_result, QUrl.fromLocalFile(CACHE_PATH))

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
		#self.table.resizeColumnsToContents()
		labels = {'ssr':'SSRs', 'vntr':'VNTRs', 'cssr':'cSSRs', 'issr':'iSSRs', 'primer':'Primers'}
		label = labels.get(count[0], 'Row')
		self.rowCounts.setText("%s: %s" % (label, count[1]))
		self.colCounts.setText("Column: %s" % count[2])

	def changeSelectCount(self, count):
		self.selectCounts.setText("Select: %s" % count)
	
	def setProgress(self, percent):
		self.progressBar.setValue(percent)

	def setStatusMessage(self, msg):
		self.statusBar.showMessage(msg)

	def openAboutMessage(self):
		#system_info = "%s%s %s" % (platform.system(), platform.release(), platform.architecture()[0])
		#python_info = sys.version.split()[0]
		about_message =	"""
			<p><b>Krait for microsatellite investigation</b></p>
			<p>Version v{version} Build {build}<p>
			<p>Krait is a robust and ultrafast tool and provides a user-friendly GUI for no computationally
			skilled biologists to extract perfect, imperfect and compound microsatellites and VNTRs with any length
			of motif from DNA fasta sequences and design PCR primers and do statistics analysis.</p>
			<p><a href="https://pypi.python.org/pypi/PySide/1.2.4">PySide</a> for GUI. 
			<a href="http://lh3lh3.users.sourceforge.net/kseq.shtml">Kseq.h</a> for parsing fasta.
			<a href="https://github.com/libnano/primer3-py">primer3-py</a> and 
			<a href="http://primer3.sourceforge.net/">primer3</a> for primer design. 
			Intersection from <a href="https://github.com/bxlab/bx-python">bx-python</a> for locating SSRs
			</p>
		""".format(version=VERSION, build=BUILD)

		QMessageBox.about(self, "About Krait", about_message)

	def openDocumentation(self):
		QDesktopServices.openUrl(QUrl("https://github.com/lmdu/niblet"))

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
		self.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.setSelectionBehavior(QAbstractItemView.SelectRows)
		self.setSelectionMode(QAbstractItemView.SingleSelection)
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

		delete_action = QAction("Delete All", self)

		detail_action = QAction("View Detail", self)
		detail_action.triggered.connect(self.viewDetail)
		
		self.menu.addAction(select_action)
		self.menu.addAction(deselect_action)
		self.menu.addAction(select_all_action)
		self.menu.addAction(deselect_all_action)
		self.menu.addSeparator()
		self.menu.addAction(delete_action)
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
		_id = self.model().getCellId(self.current_row)
		
		if table == 'primer':
			content = PrimerDetail(table, _id, flank).generateHtml()
		else:
			content = SequenceDetail(table, _id, flank).generateHtml()

		SSRDetailDialog(self.parent, "%s detail" % table.upper(), content)


class TableModel(QAbstractTableModel):
	row_col = Signal(tuple)
	sel_row = Signal(int)
	def __init__(self, parent=None):
		super(TableModel, self).__init__(parent)
		self.headers = []
		self.total = 0
		self.selected = {}
		self.read_row = 0
		self.db = Database()

		#["Select", "Where", "Order"]
		self.query = ['', '', '']

		#table name
		self.table = None

		self.cache_row_index = -1
		self.cache_row = None

	def getRowCounts(self):
		return self.total

	def getAllItems(self):
		return self.total

	def selectRow(self, row):
		if row not in self.selected:
			self.beginResetModel()
			self.selected[row] = self.getCellId(row)
			self.endResetModel()
			self.sel_row.emit(len(self.selected))

	def deselectRow(self, row):
		if row in self.selected:
			self.beginResetModel()
			self.selected.remove(row)
			self.endResetModel()
			self.sel_row.emit(len(self.selected))

	def selectAll(self):
		self.beginResetModel()
		if self.db.get_one(self.query[0] % 'COUNT(1)') == self.total:
			self.selected = {i:None for i in xrange(self.total)}
		else:
			ids = self.db.get_column(self.sql % 'id')
			self.selected = {i:j for i,j in enumerate(ids)}

		self.endResetModel()
		self.sel_row.emit(len(self.selected))

	def deselectAll(self):
		self.beginResetModel()
		self.selected = {}
		self.endResetModel()
		self.sel_row.emit(0)

	def setTable(self, table):
		self.table = table
		self.headers = self.db.get_fields(self.table)
		self.query = ["SELECT %s FROM {0}".format(self.table), '', '']

	def tableName(self):
		return self.table

	def columnNames(self):
		return self.headers

	def setFilter(self, condition=''):
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
		self.beginResetModel()
		self.read_row = 0
		self.selected = {}
		self.cache_row_index = -1
		self.cache_row = None
		self.total = self.db.get_one(self.sql % "COUNT(1)")
		self.endResetModel()
		self.row_col.emit((self.table, self.total, len(self.headers)))
		self.sel_row.emit(len(self.selected))

	def clear(self):
		self.total = 0
		self.headers = []
		self.selected = {}

	def getSelectedRows(self):
		return self.selected

	def getCellId(self, row):
		return self.db.get_one("%s LIMIT %s,1" % (self.sql % 'id', row))

	def value(self, index):
		#ID = self.dataset[index.row()]
		#col = self.headers[index.column()-1]
		col = index.column() - 1
		row = index.row()
		if row == self.cache_row_index:
			return self.cache_row[col]

		self.cache_row_index = row
		self.cache_row = self.db.get_row("%s LIMIT %s,1" % (self.sql % '*', row))
		return self.cache_row[col]

	def rowColor(self, index):
		#ID = self.dataset[index.row()]
		ID = self.db.get_one("%s LIMIT %s,1" % (self.sql % 'id', index.row()))
		sql = "SELECT feature FROM location WHERE target=%s AND category='%s' LIMIT 1" % (ID, self.table)
		feature = self.db.get_one(sql)
		if not feature:
			return QColor(255, 255, 255)
		if feature == 'CDS':
			return QColor(245, 183, 177)
		elif feature == 'EXON':
			return QColor(169, 223, 191)
		elif feature == '3UTR' or feature == '5UTR':
			return QColor(250, 215, 160) 
		elif feature == 'INTRON':
			return QColor(174, 214, 241)
		else:
			return QColor(255, 255, 255)

	def rowCount(self, parent=QModelIndex()):
		if parent.isValid():
			return 0

		return self.read_row

	def columnCount(self, parent=QModelIndex()):
		if parent.isValid():
			return 0

		if len(self.headers) > 0:
			return len(self.headers) + 1

		return len(self.headers)

	def data(self, index, role=Qt.DisplayRole):
		if not index.isValid():
			return None

		if not 0 <= index.row() < self.rowCount():
			return None

		elif role == Qt.DisplayRole:
			if index.column() > 0:
				return self.value(index)
			else:
				return None

		elif role == Qt.CheckStateRole:
			if index.column() == 0:
				if index.row() in self.selected:
					return Qt.Checked
				else:
					return Qt.Unchecked

		elif role == Qt.BackgroundColorRole:
			return self.rowColor(index)

		return None


	def headerData(self, section, orientation, role=Qt.DisplayRole):
		if role == Qt.SizeHintRole:
			if section == 0:
				return QSize(20, -1)

		if role != Qt.DisplayRole:
			return None

		if orientation == Qt.Horizontal:
			if section == 0:
				return None
			else:
				return self.headers[section-1]

		return None

	def setData(self, index, value, role):
		if not index.isValid():
			return False

		if index.column() != 0:
			return False

		if role == Qt.CheckStateRole:
			if value == Qt.Checked:
				self.selected[index.row()] = self.getCellId(index.row())
			else:
				if index.row() in self.selected:
					self.selected.remove(index.row())
			
			self.dataChanged.emit(index, index)
			self.sel_row.emit(len(self.selected))
			return True

		return False

	def flags(self, index):
		if not index.isValid():
			return QAbstractItemModel.flags(index)

		flag = Qt.ItemIsEnabled | Qt.ItemIsSelectable

		if index.column() == 0:
			return flag | Qt.ItemIsUserCheckable

		return flag

	def canFetchMore(self, parent):
		return not parent.isValid() and (self.read_row < self.total)

	def fetchMore(self, parent):
		if parent.isValid():
			return
		remainder = self.total - self.read_row
		fetch_row = min(100, remainder)
		self.beginInsertRows(QModelIndex(), self.read_row, self.read_row+fetch_row-1)
		self.read_row += fetch_row
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

		repeatsGroup = QGroupBox(self.tr("Minimal repeats for each SSR type"))
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
		distanceLabel = QLabel("Max distance (dMAX) allowed between two SSRs")
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
		max_error_label = QLabel("Max consecutive edits")
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
		self.level_select.currentIndexChanged.connect(self.showStandardLevelDetail)
		self.level_detail = QLabel()
		standard_level = [
			"Level 0  No standard",
			"Level 1  Similar motifs",
			"Level 2  Reverse complementary motifs",
			"Level 3  complementary motifs",
			"Level 4  Reverse motifs"
		]
		self.level_select.addItems(standard_level)
		level_layout = QGridLayout()
		level_layout.setColumnStretch(1, 1)
		level_layout.addWidget(level_label, 0, 0)
		level_layout.addWidget(self.level_select, 0, 1)
		level_layout.addWidget(self.level_detail, 1, 1)
		level_group.setLayout(level_layout)

		flankGroup = QGroupBox(self.tr("Flanking sequence"))
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
		self.max_error.setValue(int(self.settings.value('ssr/error', 2)))
		self.min_score.setValue(int(self.settings.value('ssr/score', 12)))
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

	def showStandardLevelDetail(self, idx):
		if idx == 0:
			detail = 'No standard'

		elif idx == 1:
			detail = 'standard similar motifs'

		elif idx == 2:
			detail = 'Level 1 + reverse complementary motifs'

		elif idx == 3:
			detail = 'Level 2 + complementary motifs'

		elif idx == 4:
			detail = 'Level 3 + reverse motifs (not recommend)'
		
		self.level_detail.setText(detail)


class PrimerTab(QWidget):
	def __init__(self, settings, parent=None):
		super(PrimerTab, self).__init__(parent)
		self.settings = settings
		
		product_size_label = PrimerTagLabel('PRIMER_PRODUCT_SIZE_RANGE')
		self.product_size = QLineEdit()
		product_size_group = QGroupBox(self.tr('Primer product size'))
		product_size_layout = QHBoxLayout()
		product_size_layout.addWidget(product_size_label)
		product_size_layout.addWidget(self.product_size, 1)
		product_size_group.setLayout(product_size_layout)

		primer_size_group = QGroupBox(self.tr("Primer size and melting temperature"))
		primer_size_layout = QGridLayout()
		self.primer_size_min = QSpinBox()
		self.primer_size_min.setSuffix(' bp')
		self.primer_size_opt = QSpinBox()
		self.primer_size_opt.setSuffix(' bp')
		self.primer_size_max = QSpinBox()
		self.primer_size_max.setSuffix(' bp')
		self.primer_tm_min = QSpinBox()
		self.primer_tm_min.setSuffix(' %sC' % chr(0260))
		self.primer_tm_opt = QSpinBox()
		self.primer_tm_opt.setSuffix(' %sC' % chr(0260))
		self.primer_tm_max = QSpinBox()
		self.primer_tm_max.setSuffix(' %sC' % chr(0260))
		primer_size_layout.addWidget(PrimerTagLabel("PRIMER_MIN_SIZE"), 0, 0)
		primer_size_layout.addWidget(self.primer_size_min, 0, 1)
		primer_size_layout.addWidget(PrimerTagLabel("PRIMER_OPT_SIZE"), 0, 2)
		primer_size_layout.addWidget(self.primer_size_opt, 0, 3)
		primer_size_layout.addWidget(PrimerTagLabel("PRIMER_MAX_SIZE"), 0, 4)
		primer_size_layout.addWidget(self.primer_size_max, 0, 5)
		primer_size_layout.addWidget(PrimerTagLabel("PRIMER_MIN_TM"), 1, 0)
		primer_size_layout.addWidget(self.primer_tm_min, 1, 1)
		primer_size_layout.addWidget(PrimerTagLabel("PRIMER_OPT_TM"), 1, 2)
		primer_size_layout.addWidget(self.primer_tm_opt, 1, 3)
		primer_size_layout.addWidget(PrimerTagLabel("PRIMER_MAX_TM"), 1, 4)
		primer_size_layout.addWidget(self.primer_tm_max, 1, 5)
		primer_size_group.setLayout(primer_size_layout)

		primer_gc_group = QGroupBox(self.tr("Primer GC content"))
		primer_gc_layout = QGridLayout()
		self.primer_gc_min = QSpinBox()
		self.primer_gc_min.setSuffix(' %')
		self.primer_gc_max = QSpinBox()
		self.primer_gc_max.setSuffix(' %')
		self.primer_gc_clamp = QSpinBox()
		self.primer_gc_end = QSpinBox()
		primer_gc_layout.addWidget(PrimerTagLabel("PRIMER_MIN_GC"), 0, 0)
		primer_gc_layout.addWidget(self.primer_gc_min, 0, 1)
		primer_gc_layout.addWidget(PrimerTagLabel("PRIMER_GC_CLAMP"), 0, 2)
		primer_gc_layout.addWidget(self.primer_gc_clamp, 0, 3)
		primer_gc_layout.addWidget(PrimerTagLabel("PRIMER_MAX_GC"), 1, 0)
		primer_gc_layout.addWidget(self.primer_gc_max, 1, 1)
		primer_gc_layout.addWidget(PrimerTagLabel("PRIMER_PAIR_MAX_DIFF_TM"), 1, 2)
		primer_gc_layout.addWidget(self.primer_gc_end, 1, 3)
		primer_gc_group.setLayout(primer_gc_layout)

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
		self.primer_num_return = QSpinBox()
		primer_other_layout.addWidget(PrimerTagLabel("PRIMER_MAX_END_STABILITY"), 0, 0)
		primer_other_layout.addWidget(self.primer_max_end_stability, 0, 1)
		primer_other_layout.addWidget(PrimerTagLabel("PRIMER_MAX_POLY_X"), 0, 2)
		primer_other_layout.addWidget(self.primer_max_poly_x, 0, 3)
		primer_other_layout.addWidget(PrimerTagLabel("PRIMER_MAX_NS_ACCEPTED"), 1, 0)
		primer_other_layout.addWidget(self.primer_max_ns_accepted, 1, 1)
		primer_other_layout.addWidget(PrimerTagLabel("PRIMER_NUM_RETURN"), 1, 2)
		primer_other_layout.addWidget(self.primer_num_return, 1, 3)
		primer_other_group.setLayout(primer_other_layout)

		mainLayout = QVBoxLayout()
		mainLayout.addWidget(product_size_group)
		mainLayout.addWidget(primer_size_group)
		mainLayout.addWidget(primer_gc_group)
		mainLayout.addWidget(primer_bind_group)
		mainLayout.addWidget(primer_other_group)

		self.setLayout(mainLayout)
		self.getSettings()

	def getSettings(self):
		self.product_size.setText(self.settings.value('primer/PRIMER_PRODUCT_SIZE_RANGE', '100-300'))
		self.primer_size_min.setValue(int(self.settings.value('primer/PRIMER_MIN_SIZE', 18)))
		self.primer_size_opt.setValue(int(self.settings.value('primer/PRIMER_OPT_SIZE', 20)))
		self.primer_size_max.setValue(int(self.settings.value('primer/PRIMER_MAX_SIZE', 27)))
		self.primer_tm_min.setValue(int(self.settings.value('primer/PRIMER_MIN_TM', 57)))
		self.primer_tm_opt.setValue(int(self.settings.value('primer/PRIMER_OPT_TM', 60)))
		self.primer_tm_max.setValue(int(self.settings.value('primer/PRIMER_MAX_TM', 63)))
		self.primer_gc_min.setValue(int(self.settings.value('primer/PRIMER_MIN_GC', 20)))
		self.primer_gc_max.setValue(int(self.settings.value('primer/PRIMER_MAX_GC', 80)))
		self.primer_gc_clamp.setValue(int(self.settings.value('primer/PRIMER_GC_CLAMP', 0)))
		self.primer_gc_end.setValue(int(self.settings.value('primer/PRIMER_PAIR_MAX_DIFF_TM', 100)))
		self.primer_max_self_any.setValue(int(self.settings.value('primer/PRIMER_MAX_SELF_ANY_TH', 47)))
		self.primer_pair_max_compl_any.setValue(int(self.settings.value('primer/PRIMER_PAIR_MAX_COMPL_ANY_TH', 47)))
		self.primer_max_self_end.setValue(int(self.settings.value('primer/PRIMER_MAX_SELF_END_TH', 47)))
		self.primer_pair_max_compl_end.setValue(int(self.settings.value('primer/PRIMER_PAIR_MAX_COMPL_END_TH', 47)))
		self.primer_max_hairpin.setValue(int(self.settings.value('primer/PRIMER_MAX_HAIRPIN_TH', 47)))
		self.primer_max_end_stability.setValue(int(self.settings.value('primer/PRIMER_MAX_END_STABILITY', 100)))
		self.primer_max_ns_accepted.setValue(int(self.settings.value('primer/PRIMER_MAX_POLY_X', 5)))
		self.primer_max_poly_x.setValue(int(self.settings.value('primer/PRIMER_MAX_NS_ACCEPTED', 0)))
		self.primer_num_return.setValue(int(self.settings.value('primer/PRIMER_NUM_RETURN', 5)))


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
		self.settings.setValue('primer/PRIMER_PAIR_MAX_DIFF_TM', self.primer_gc_end.value())
		self.settings.setValue('primer/PRIMER_MAX_SELF_ANY_TH', self.primer_max_self_any.value())
		self.settings.setValue('primer/PRIMER_PAIR_MAX_COMPL_ANY_TH', self.primer_pair_max_compl_any.value())
		self.settings.setValue('primer/PRIMER_MAX_SELF_END_TH', self.primer_max_self_end.value())
		self.settings.setValue('primer/PRIMER_PAIR_MAX_COMPL_END_TH', self.primer_pair_max_compl_end.value())
		self.settings.setValue('primer/PRIMER_MAX_HAIRPIN_TH', self.primer_max_hairpin.value())
		self.settings.setValue('primer/PRIMER_MAX_END_STABILITY', self.primer_max_end_stability.value())
		self.settings.setValue('primer/PRIMER_MAX_POLY_X', self.primer_max_ns_accepted.value())
		self.settings.setValue('primer/PRIMER_MAX_NS_ACCEPTED', self.primer_max_poly_x.value())
		self.settings.setValue('primer/PRIMER_NUM_RETURN', self.primer_num_return.value())


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
		self.viewer = QWebView(self)
		self.viewer.setHtml(content, QUrl.fromLocalFile(CACHE_PATH))

		buttonBox = QDialogButtonBox(QDialogButtonBox.Ok)
		buttonBox.accepted.connect(self.accept)
		
		mainLayout = QVBoxLayout()
		mainLayout.addWidget(self.viewer)
		mainLayout.addWidget(buttonBox)
		self.resize(600, 400)

		self.setLayout(mainLayout)
		self.exec_()

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

		layout.addWidget(buttonBox, 4, 0, 1)

		self.setLayout(layout)

	def select(self):
		exp_file, _ = QFileDialog.getSaveFileName(self, filter="Fasta (*.fa);;Fasta (*.fasta)")
		if exp_file:
			self.out_input.setText(exp_file)

	def get(self):
		acc = self.acc_input.text().strip()
		out = self.out_input.text()
		return acc, out


#class SSRTableModel(QSqlTableModel):
#	def __init__(self):
#		super(SSRTableModel, self).__init__()
#		self.setTable('ssr')
#		self.select()

	#def data(self, index, role=Qt.DisplayRole):
	#	if role == Qt.TextAlignmentRole:
	#		return Qt.AlignCenter

	#	return QSqlTableModel.data(self, index, role)