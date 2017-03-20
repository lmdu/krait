#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import apsw
import platform

from PySide.QtCore import *
from PySide.QtGui import *
from PySide.QtSql import *
from PySide.QtWebKit import *

from db import *
from fasta import *
from utils import *
from workers import *
from config import *
from primer import *

class SSRMainWindow(QMainWindow):
	def __init__(self):
		super(SSRMainWindow, self).__init__()

		self.setWindowTitle("Niblet v0.0.1")
		self.setWindowIcon(QIcon('logo.ico'))

		#self.browser = SSRWebView()
		
		self.table = SSRTableView()
		#self.table.verticalHeader().hide()
		#self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)s
		##self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
		#self.table.setSortingEnabled(True)
		#self.table.setContextMenuPolicy(Qt.CustomContextMenu)
		#self.table.doubleClicked.connect(self.showSSRSequence)
		#self.model = SSRTableModel()
		self.model = TableModel()
		self.table.setModel(self.model)
		self.model.row_col.connect(self.changeRowColCount)
		self.model.sel_row.connect(self.changeSelectCount)
		self.setCentralWidget(self.table)
		#self.setCentralWidget(self.browser)

		self.reportor = QTextBrowser()
		self.reportor.setSearchPaths([CACHE_PATH])
		self.reportor.setOpenLinks(False)
		self.reportor.setOpenExternalLinks(False)
		self.reportor.anchorClicked.connect(self.saveStatTableFigure)

		#search text input
		self.filter = SSRFilterInput(self)
		self.filter.returnPressed.connect(self.filterTable)

		#create fasta table
		self.fasta_table = FastaTable()

		self.createActions()
		self.createMenus()
		self.createToolBars()
		self.createStatusBar()

		#connect to database
		self.db = Database()

		#read settings
		self.readSettings()

		self.show()

	def saveStatTableFigure(self, url):
		action = url.toString()
		name = QFileDialog.getSaveFileName(self, directory=action, filter="SVG image (*.svg)")
		if not name: return
		print name

	def setDatabase(self, db_name):
		try:
			db = open_database(db_name)
		except Exception, e:
			QMessageBox.critical(self,
				self.tr('Error occurred'),
				self.tr(str(e))
			)
			return False
		return True

	def readSettings(self):
		self.settings = QSettings("config.ini", QSettings.IniFormat)
		self.resize(self.settings.value("size", QSize(900, 600)))

	def writeSettings(self):
		self.settings.setValue("size", self.size())

	def closeEvent(self, event):
		self.writeSettings()

	def createActions(self):
		#open a project action
		self.openProjectAct = QAction(self.tr("Open project"), self)
		self.openProjectAct.setShortcut(QKeySequence.Open)
		self.openProjectAct.triggered.connect(self.openProject)

		#close a project action
		self.closeProjectAct = QAction(self.tr("Close project"), self)
		self.closeProjectAct.setShortcut(QKeySequence.Close)
		self.closeProjectAct.triggered.connect(self.closeProject)
		
		#save a project action
		self.saveProjectAct = QAction(self.tr("Save project"), self)
		self.saveProjectAct.setShortcut(QKeySequence.Save)
		self.saveProjectAct.triggered.connect(self.saveProject)
		
		#save as a project action
		self.saveAsProjectAct = QAction(self.tr("Save project as..."), self)
		self.saveAsProjectAct.setShortcut(QKeySequence.SaveAs)
		self.saveAsProjectAct.triggered.connect(self.saveProjectAs)
		
		#load fasta file or genome action
		self.loadFastaAct = QAction(self.tr("Import fasta sequence"), self)
		self.loadFastaAct.triggered.connect(self.importFasta)
		self.loadFastasAct = QAction(self.tr("Import Fastas in folder"), self)
		self.loadFastasAct.triggered.connect(self.importFastas)
		
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
		self.preferenceAct.triggered.connect(self.setPreference)

		#toolbar actions
		#search perfect ssrs tool button
		self.perfectAct = QAction(QIcon("icons/ssr.png"), self.tr("SSRs"), self)
		self.perfectAct.setToolTip(self.tr("Search perfect microsatellites"))
		self.perfectAct.triggered.connect(self.searchMicrosatellites)
		self.perfectMenuAct = QAction(self.tr("Perform SSR search"), self)
		self.perfectMenuAct.triggered.connect(self.searchMicrosatellites)
		self.perfectResultAct = QAction(self.tr("Show perfect SSRs"), self)
		self.perfectResultAct.triggered.connect(self.showMicrosatellites)
		self.perfectRemoveAct = QAction(self.tr("Remove perfect SSRs"), self)
		self.perfectRemoveAct.triggered.connect(self.removePerfectSSRs)
		self.minRepeatAct = QAction(self.tr("Set minimum repeats"), self)
		self.minRepeatAct.triggered.connect(self.setPreference)
		
		#search compound ssrs tool button
		self.compoundAct = QAction(QIcon("icons/cssr.png"), self.tr("cSSRs"), self)
		self.compoundAct.setToolTip(self.tr("Identify compound microsatellites using dMax"))
		self.compoundAct.triggered.connect(self.searchCompoundSSRs)
		self.compoundMenuAct = QAction(self.tr("Perform cSSRs search"), self)
		self.compoundMenuAct.triggered.connect(self.searchCompoundSSRs)
		self.compoundResultAct = QAction(self.tr("Show compound SSRs"), self)
		self.compoundResultAct.triggered.connect(self.showCompoundSSRs)
		self.compoundRemoveAct = QAction(self.tr("Remove cSSR results"), self)
		self.compoundRemoveAct.triggered.connect(self.removeCompoundSSRs)
		self.bestDmaxAct = QAction(self.tr("Estimate best dMax"), self)
		self.bestDmaxAct.triggered.connect(self.estimateBestMaxDistance)
		self.maxDistanceAct = QAction(self.tr("Set Maximum distance dMAX"), self)
		self.maxDistanceAct.triggered.connect(self.setPreference)

		#search satellite dna
		self.satelliteAct = QAction(QIcon("icons/satellite.png"), self.tr("VNTRs"), self)
		self.satelliteAct.setToolTip(self.tr("Detect satellites"))
		self.satelliteAct.triggered.connect(self.detectSatellites)
		self.satelliteMenuAct = QAction(self.tr("Detect satellites"), self)
		self.satelliteMenuAct.triggered.connect(self.detectSatellites)
		self.satelliteResultAct = QAction(self.tr("Show satellites"), self)
		self.satelliteResultAct.triggered.connect(self.showSatellites)
		self.satelliteRemoveAct = QAction(self.tr("Remove satellites"), self)
		self.satelliteRemoveAct.triggered.connect(self.removeSatellites)
		self.satelliteRuleAct = QAction(self.tr("Set search parameter"), self)
		self.satelliteRuleAct.triggered.connect(self.setPreference)

		#search imperfect microsatellites
		self.imperfectAct = QAction(QIcon("icons/issr.png"), self.tr("iSSRs"), self)
		self.imperfectAct.setToolTip(self.tr("Find imperfect microsatellites"))

		#design primer
		self.primerAct = QAction(QIcon("icons/primer.png"), self.tr("Primer"), self)
		self.primerAct.setToolTip(self.tr("Design primers"))
		self.primerAct.triggered.connect(self.designPrimers)
		self.primerAllAct = QAction(self.tr("Design primer for all rows"), self)
		self.primerAllAct.setToolTip(self.tr("Design primer for all rows in table"))
		self.primerAllAct.triggered.connect(self.designPrimerForTable)
		self.primerSelectAct = QAction(self.tr("Design primer for selected rows"), self)
		self.primerSelectAct.setToolTip(self.tr("Design primer for selected rows in table"))
		self.primerSelectAct.triggered.connect(self.designPrimerForSelected)
		self.primerShowAct = QAction(self.tr("Show designed primers"), self)
		self.primerRemoveAct = QAction(self.tr("Remove designed primers"), self)

		#statistics report
		self.statisticsAct = QAction(QIcon("icons/report.png"), self.tr("Statistics"), self)
		self.statisticsAct.triggered.connect(self.generateStatisticsReport)
		self.statisticsMenuAct = QAction(self.tr("Perform statistics"), self)
		self.statisticsMenuAct.triggered.connect(self.generateStatisticsReport)

		#about action
		self.aboutAct = QAction(self.tr("About"), self)
		self.aboutAct.triggered.connect(self.openAboutMessage)

		#documentation action
		self.documentAct = QAction(self.tr("Documentation"), self)
		self.documentAct.triggered.connect(self.openDocumentation)
		

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

		self.searchMenu.addAction(self.perfectMenuAct)
		self.searchMenu.addAction(self.compoundMenuAct)
		self.searchMenu.addAction(self.satelliteMenuAct)

		self.viewMenu.addAction(self.perfectResultAct)
		self.viewMenu.addAction(self.compoundResultAct)
		self.viewMenu.addAction(self.satelliteResultAct)
		self.viewMenu.addSeparator()
		self.viewMenu.addAction(self.perfectRemoveAct)
		self.viewMenu.addAction(self.compoundRemoveAct)
		self.viewMenu.addAction(self.satelliteRemoveAct)

		self.toolMenu.addAction(self.bestDmaxAct)
		self.toolMenu.addAction(self.statisticsAct)

		self.helpMenu.addAction(self.documentAct)
		self.helpMenu.addSeparator()
		self.helpMenu.addAction(self.aboutAct)


		#tool bar menus
		#search ssrs tool button menu
		self.perfectMenu = QMenu()
		self.perfectMenu.addAction(self.perfectMenuAct)
		self.perfectMenu.addAction(self.perfectResultAct)
		self.perfectMenu.addAction(self.perfectRemoveAct)
		self.perfectMenu.addSeparator()
		self.perfectMenu.addAction(self.loadFastaAct)
		self.perfectMenu.addAction(self.loadFastasAct)
		self.perfectMenu.addSeparator()
		self.perfectMenu.addAction(self.minRepeatAct)

		self.compoundMenu = QMenu()
		self.compoundMenu.addAction(self.compoundMenuAct)
		self.compoundMenu.addAction(self.compoundResultAct)
		self.compoundMenu.addAction(self.compoundRemoveAct)
		self.compoundMenu.addSeparator()
		self.compoundMenu.addAction(self.bestDmaxAct)
		self.compoundMenu.addAction(self.maxDistanceAct)

		self.satelliteMenu = QMenu()
		self.satelliteMenu.addAction(self.satelliteMenuAct)
		self.satelliteMenu.addAction(self.satelliteResultAct)
		self.satelliteMenu.addAction(self.satelliteRemoveAct)
		self.satelliteMenu.addSeparator()
		self.satelliteMenu.addAction(self.satelliteRuleAct)

		self.imperfectMenu = QMenu()

		self.primerMenu = QMenu()
		self.primerMenu.addAction(self.primerAllAct)
		self.primerMenu.addAction(self.primerSelectAct)
		self.primerMenu.addAction(self.primerShowAct)
		self.primerMenu.addAction(self.primerRemoveAct)

		self.statisticsMenu = QMenu()
		self.statisticsMenu.addAction(self.statisticsMenuAct)
		

	def createToolBars(self):
		self.toolBar = self.addToolBar('')
		self.toolBar.setMovable(False)
		self.toolBar.setIconSize(QSize(36, 36))
		self.toolBar.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)

		#search ssr action and menus
		self.perfectAct.setMenu(self.perfectMenu)
		self.toolBar.addAction(self.perfectAct)

		self.compoundAct.setMenu(self.compoundMenu)
		self.toolBar.addAction(self.compoundAct)

		self.satelliteAct.setMenu(self.satelliteMenu)
		self.toolBar.addAction(self.satelliteAct)

		self.imperfectAct.setMenu(self.imperfectMenu)
		self.toolBar.addAction(self.imperfectAct)

		self.primerAct.setMenu(self.primerMenu)
		self.toolBar.addAction(self.primerAct)

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
		self.statusBar.addPermanentWidget(self.progressBar)
		

	def openProject(self):
		dbfile, _ = QFileDialog.getOpenFileName(self, filter="Database (*.db)")
		if not dbfile:
			return

		self.db.drop_tables()
		self.db.open(dbfile)

		self.showMicrosatellites()

	def saveProject(self):
		dbfile, _ = QFileDialog.getSaveFileName(self, filter="Database (*.db)")
		if not dbfile:
			return

		self.db.save(dbfile)

	def saveProjectAs(self):
		dbfile, _ = QFileDialog.getSaveFileName(self, filter="Database (*.db)")
		if not dbfile:
			return

		self.db.save(dbfile)

	def closeProject(self):
		self.db.drop_tables()
		self.model.select()
	
	def importFasta(self):
		'''
		Import a fasta file from a directory
		'''
		fasta, _ = QFileDialog.getOpenFileName(self, filter="Fasta (*.fa *.fna *.fas *.fasta);;All files (*.*)")
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

	def execute(self, worker):
		self.thread = QThread()
		self.thread.started.connect(worker.process)
		worker.finished.connect(self.thread.quit)
		worker.finished.connect(worker.deleteLater)
		self.thread.finished.connect(self.thread.deleteLater)
		worker.moveToThread(self.thread)
		self.thread.start()

	def searchMicrosatellites(self):
		#if self.db.isTableExists('ssr'):
		#	status = QMessageBox.warning(self, 
		#		self.tr("Warning"), 
		#		self.tr("The SSRs have been searched.\nDo you want to remove result and search again?"), 
		#		QMessageBox.Ok | QMessageBox.Cancel
		#	)

		#	if status == QMessageBox.Cancel:
			
		rules = [
			int(self.settings.value('ssr/mono')),
			int(self.settings.value('ssr/di')),
			int(self.settings.value('ssr/tri')), 
			int(self.settings.value('ssr/tetra')),
			int(self.settings.value('ssr/penta')),
			int(self.settings.value('ssr/hexa'))
		]
		level = int(self.settings.value('ssr/level'))
		fastas = self.db.get_all("SELECT * FROM fasta")
		self.worker = SSRWorker(fastas, rules, level)
		self.worker.update_message.connect(self.setStatusMessage)
		self.worker.update_progress.connect(self.setProgress)
		self.worker.finished.connect(self.showMicrosatellites)
		self.execute(self.worker)
	
	def showMicrosatellites(self):
		self.model.setTable('ssr')
		#self.table.horizontalHeader().setResizeMode(QHeaderView.Stretch)
		#self.model.refresh()
		self.model.select()

	def removePerfectSSRs(self):
		pass

	def searchCompoundSSRs(self):
		dmax = int(self.settings.value('ssr/dmax', 10))
		worker = CompoundWorker(self, dmax)
		worker.update_message.connect(self.setStatusMessage)
		worker.update_progress.connect(self.setProgress)
		worker.finished.connect(self.showCompoundSSRs)
		worker.start()

	def showCompoundSSRs(self):
		self.model.setTable('cssr')
		self.table.horizontalHeader().setResizeMode(QHeaderView.Interactive)
		self.model.refresh()

	def removeCompoundSSRs(self):
		pass

	def detectSatellites(self):
		min_motif = int(self.settings.value('ssr/smin'))
		max_motif = int(self.settings.value('ssr/smax'))
		min_repeat = int(self.settings.value('ssr/srep'))

		fastas = [fasta for fasta in self.fasta_table.fetchAll()]
		worker = SatelliteWorker(self, fastas, min_motif, max_motif, min_repeat)
		worker.update_message.connect(self.setStatusMessage)
		worker.update_progress.connect(self.setProgress)
		worker.finished.connect(self.showSatellites)
		worker.start()

	def showSatellites(self):
		self.model.setTable('satellite')
		self.model.refresh()

	def removeSatellites(self):
		pass

	def getPrimer3Settings(self):
		p3_settings = dict(
			PRIMER_TASK = 'generic',
			PRIMER_PICK_LEFT_PRIMER = 1,
			PRIMER_PICK_INTERNAL_OLIGO = 0,
			PRIMER_PICK_RIGHT_PRIMER = 1
		)
		self.settings.beginGroup("primer")
		keys = self.settings.childKeys()
		for key in keys:
			p3_settings[key] = self.settings.value(key)
		return p3_settings

	def designPrimers(self):
		rows = self.model.dataset
		table = self.model.table
		flank = min_motif = int(self.settings.value('ssr/flank'))
		primer3_settings = self.getPrimer3Settings()
		self.worker = PrimerWorker(table, rows, flank, primer3_settings)
		self.execute(self.worker)

	def designPrimerForTable(self):
		rows = self.model.dataset
		table = self.model.table


	def designPrimerForSelected(self):
		rows = self.model.getSelectedRows()
		table = self.model.table




	def showPrimers(self):
		pass

	def estimateBestMaxDistance(self):
		pass

	def filterTable(self):
		filters = str(self.filter.text())

		if filters.startswith('db'):
			self.model.setTable(filters.split('=')[1])
			return
		self.model.setFilter(filters)

	def generateStatisticsReport(self):
		worker = StatisticsWorker(self.reportor, self)
		worker.update_message.connect(self.showStatisticsReport)
		worker.start()

	def showStatisticsReport(self, html):
		self.reportor.setHtml(html)
		self.setCentralWidget(self.reportor)

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
		self.rowCounts.setText("Row: %s" % count[0])
		self.colCounts.setText("Column: %s" % count[1])

	def changeSelectCount(self, count):
		self.selectCounts.setText("Select: %s" % count)
	
	def setProgress(self, percent):
		self.progressBar.setValue(percent)

	def setStatusMessage(self, msg):
		self.statusBar.showMessage(msg)

	def openAboutMessage(self):
		system_info = "%s%s %s" % (platform.system(), platform.release(), platform.architecture()[0])
		python_info = sys.version.split()[0]
		about_message =	"""
			<p><b>Niblet for finding tandem repeats</b></p>
			<p>Version v0.1.0 Build 20170104<p>
			<p>System {system} Python {python}</p>
			<p>Niblet is a user-friendly tool that allows user to extract perfect microsatellites,
			compound microsatellites, imperfect microsatellites and tandem repeats with any length
			of motif from DNA fasta sequences and batch design PCR primers and statistics analysis.</p>
			<p>GUI was written by <a href="https://pypi.python.org/pypi/PySide/1.2.4">PySide</a>. 
			Fasta sequences are extracted by using <a href="https://github.com/mdshw5/pyfaidx">pyfaidx</a>.
			Primer design are performed by using <a href="https://github.com/libnano/primer3-py">primer3-py</a>.
			</p>
			<p><b>Cite:</b> Du L. Liu Q. and Yue B. Niblet: a flexible tool for genome-wide survey
			of tandem repeats.2017</p>
		""".format(system=system_info, python=python_info)

		QMessageBox.about(self, "About niblet", about_message)

	def openDocumentation(self):
		QDesktopServices.openUrl(QUrl("https://github.com/lmdu/niblet"))

class SSRFilterInput(QLineEdit):
	def __init__(self, parent=None):
		super(SSRFilterInput, self).__init__(parent)
		self.setPlaceholderText("Filter data in table e.g. motif=AT and repeat>10")


class SSRWebView(QWebView):
	def __init__(self, parent=None):
		super(SSRWebView, self).__init__(parent)
		self.load(QUrl('http://www.baidu.com'))
		self.page().setLinkDelegationPolicy(QWebPage.DelegateAllLinks)
		self.linkClicked.connect(self.openUrl)

	def openUrl(self, url):
		QDesktopServices.openUrl(url)


class SSRTableModel(QSqlTableModel):
	refreshed = Signal()
	def __init__(self):
		super(SSRTableModel, self).__init__()

	def refresh(self):
		self.select()
		self.refreshed.emit()


class SSRTableView(QTableView):
	def __init__(self, parent=None):
		super(SSRTableView, self).__init__(parent)
		self.verticalHeader().hide()
		self.horizontalHeader().setHighlightSections(False)
		self.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.setSelectionBehavior(QAbstractItemView.SelectRows)
		self.setSelectionMode(QAbstractItemView.SingleSelection)
		self.setSortingEnabled(True)

class TableModel(QAbstractTableModel):
	row_col = Signal(tuple)
	sel_row = Signal(int)
	def __init__(self, parent=None):
		super(TableModel, self).__init__(parent)
		self.headers = []
		self.dataset = []
		self.selected = set()
		self.read_row = 0
		self.db = Database()
		self.query = ['', '', '']

	def setTable(self, table):
		self.table = table
		self.headers = self.db.get_columns(self.table)
		self.query[0] = "SELECT id FROM %s" % self.table

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
		
		col = self.headers[column-1]
		if order == Qt.SortOrder.DescendingOrder:
			self.query[2] = "ORDER BY %s DESC" % (col)
			self.select()
		else:
			self.query[2] = ''
			self.select()

	def select(self):
		sql = " ".join(self.query)
		self.beginResetModel()
		self.dataset = self.db.get_column(sql)
		self.read_row = 0
		self.selected = set()
		self.endResetModel()
		self.row_col.emit((len(self.dataset), len(self.headers)))

	def clear(self):
		self.dataset = []
		self.headers = []
		self.selected = set()

	def getSelectedRows(self):
		if self.selected:
			return [self.dataset[idx] for idx in self.selected]

	def value(self, index):
		ID = self.dataset[index.row()]
		col = self.headers[index.column()-1]
		return self.db.get_one("SELECT %s FROM %s WHERE id=%s" % (col, self.table, ID))

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

		return None


	def headerData(self, section, orientation, role=Qt.DisplayRole):
		if role == Qt.SizeHintRole:
			if section == 0:
				return QSize(20, -1)

		if role != Qt.DisplayRole:
			return None

		if orientation == Qt.Horizontal:
			if section == 0:
				return ''
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
				self.selected.add(index.row())
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
		return self.read_row < len(self.dataset)

	def fetchMore(self, parent):
		remainder = len(self.dataset) - self.read_row
		fetch_row = min(100, remainder)
		self.beginInsertRows(QModelIndex(), self.read_row, self.read_row+fetch_row)
		self.read_row += fetch_row
		self.endInsertRows()

class PreferenceDialog(QDialog):
	def __init__(self, parent=None, settings=None):
		super(PreferenceDialog, self).__init__(parent)
		self.settings = settings
		self.setWindowTitle(self.tr("Preferences"))
		#self.setMinimumWidth(500)

		self.general_tab = GeneralTab(self.settings)
		self.primer_tab = PrimerTab(self.settings)

		tabWidget = QTabWidget()
		tabWidget.addTab(self.general_tab, 'SSR search')
		tabWidget.addTab(self.primer_tab, 'Primer design')

		buttonBox = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
		buttonBox.accepted.connect(self.accept)
		buttonBox.rejected.connect(self.reject)

		spacerItem = QSpacerItem(10, 20, QSizePolicy.Minimum, QSizePolicy.Expanding)

		mainLayout = QVBoxLayout()
		mainLayout.addWidget(tabWidget)
		mainLayout.addItem(spacerItem)
		mainLayout.addWidget(buttonBox)

		self.setLayout(mainLayout)

	def saveSettings(self):
		self.general_tab.saveSettings()
		self.primer_tab.saveSettings()


class GeneralTab(QWidget):
	def __init__(self, settings, parent=None):
		super(GeneralTab, self).__init__(parent)
		self.settings = settings

		repeatsGroup = QGroupBox(self.tr("Microsatellite minimum repeats"))
		monoLabel = QLabel("Mono-nucleotide")
		self.monoValue = QSpinBox()
		self.monoValue.setSuffix(' bp')
		diLabel = QLabel("Di-nucleotide")
		self.diValue = QSpinBox()
		self.diValue.setSuffix(' bp')
		triLabel = QLabel("Tri-nucleotide")
		self.triValue = QSpinBox()
		self.triValue.setSuffix(' bp')
		tetraLabel = QLabel("Tetra-nucleotide")
		self.tetraValue = QSpinBox()
		self.tetraValue.setSuffix(' bp')
		pentaLabel = QLabel("Penta-nucleotide")
		self.pentaValue = QSpinBox()
		self.pentaValue.setSuffix(' bp')
		hexaLabel = QLabel("Hexa-nucleotide")
		self.hexaValue = QSpinBox()
		self.hexaValue.setSuffix(' bp')
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

		distanceGroup = QGroupBox(self.tr("Compound microsatellite"))
		distanceLabel = QLabel("Maximal allowed distance between two SSRs (dMAX): ")
		self.distanceValue = QSpinBox()
		self.distanceValue.setSuffix(' bp')
		distanceLayout = QHBoxLayout()
		distanceLayout.addWidget(distanceLabel)
		distanceLayout.addWidget(self.distanceValue, 1)
		distanceGroup.setLayout(distanceLayout)

		satelliteGroup = QGroupBox(self.tr("Satellite"))
		min_tandem_label = QLabel("Minimum length of repeat unit ")
		self.min_tandem_motif = QSpinBox()
		self.min_tandem_motif.setMinimum(7)
		self.min_tandem_motif.setSuffix(' bp')

		max_tandem_label = QLabel("Maximum ")
		self.max_tandem_motif = QSpinBox()
		self.max_tandem_motif.setMinimum(7)
		self.max_tandem_motif.setSuffix(' bp')

		repeat_tandem_label = QLabel("Minimum allowed repeats ")
		self.min_tandem_repeat = QSpinBox()
		self.min_tandem_repeat.setMinimum(2)

		satelliteLayout = QGridLayout()
		satelliteLayout.addWidget(min_tandem_label, 0, 0)
		satelliteLayout.addWidget(self.min_tandem_motif, 0, 1)
		satelliteLayout.addWidget(max_tandem_label, 0, 2)
		satelliteLayout.addWidget(self.max_tandem_motif, 0, 3)
		satelliteLayout.addWidget(repeat_tandem_label, 1, 0)
		satelliteLayout.addWidget(self.min_tandem_repeat, 1, 1)
		satelliteGroup.setLayout(satelliteLayout)

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
		flankLabel = QLabel("Flanking sequence length: ")
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
		self.min_tandem_motif.setValue(int(self.settings.value('ssr/smin', 7)))
		self.max_tandem_motif.setValue(int(self.settings.value('ssr/smax', 30)))
		self.min_tandem_repeat.setValue(int(self.settings.value('ssr/srep', 2)))
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
		self.settings.setValue('ssr/smin', self.min_tandem_motif.value())
		self.settings.setValue('ssr/smax', self.max_tandem_motif.value())
		self.settings.setValue('ssr/srep', self.min_tandem_repeat.value())
		self.settings.setValue('ssr/level', self.level_select.currentIndex())


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



class BrowserDialog(QDialog):
	def __init__(self, parent=None, html=None):
		super(BrowserDialog, self).__init__(parent)
		self.webview = QTextEdit(self)
		self.webview.setHtml(html)

		buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		buttonBox.accepted.connect(self.accept)
		buttonBox.rejected.connect(self.reject)
		
		mainLayout = QVBoxLayout()
		mainLayout.addWidget(self.webview)
		mainLayout.addWidget(buttonBox)
		self.resize(600, 400)

		self.setLayout(mainLayout)



#class SSRTableModel(QSqlTableModel):
#	def __init__(self):
#		super(SSRTableModel, self).__init__()
#		self.setTable('ssr')
#		self.select()

	#def data(self, index, role=Qt.DisplayRole):
	#	if role == Qt.TextAlignmentRole:
	#		return Qt.AlignCenter

	#	return QSqlTableModel.data(self, index, role)