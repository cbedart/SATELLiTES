####################################################################################################
# SATELLiTES - version 1.0.10
# Copyright (C) 2024 Corentin BEDART
####################################################################################################

import sys
import os
import datetime

from PyQt6.QtWidgets import *
from PyQt6.QtCore import Qt, QTimer, QSize, QFile, QUrl
from PyQt6.QtGui import QFont, QPixmap, QImage, QTextCursor, QIcon, QDesktopServices

from rdkit import Chem
from rdkit.Chem import Draw, Descriptors

import multiprocessing as mp

import SATELLiTES
from SATELLiTES.SATELLiTES_gui import *
from SATELLiTES.SATELLiTES_tabs_2reagents import *
from SATELLiTES.SATELLiTES_tabs_3reagents import *

installdir_SATELLiTES = os.path.dirname(os.path.abspath(SATELLiTES.__file__)) + "/"

##################################################################################################################################

def launch_function(reagentsA, reagentsB, reagentsC, reaction, output_path, nb_cpu, modeA, modeB, modeC, inputA, inputB, inputC, step, step_path):
    # # Separate function that uses the information from the text fields
    # print("Launching function with the following inputs:")
    # print("Reagents A path =", reagentsA)
    # print("Reagents B path =", reagentsB)
    # print("Reagents C path =", reagentsC)

    # print("ModeA", modeA, "with additional input = ", inputA)
    # print("ModeB", modeB, "with additional input = ", inputB)
    # print("ModeC", modeC, "with additional input = ", inputC)

    # print("Step #", step)
    # print("SMARTS reaction =", reaction)
    # print("Output path =", output_path)
    # print("Number of threads =", nb_cpu)

    # if step == 2:
    #     print("Selected compounds from step #1 =", step_path)
    # if step == 3:
    #     print("Selected compounds from step #2 =", step_path)

    # Main function
    os.chdir(output_path)

    # job_name setup
    timestamp = datetime.datetime.now()
    
    # return ""
    if modeC == "":
        CB_run_2reagents(reagentsA, reagentsB, reaction, output_path, nb_cpu, modeA, modeB, inputA, inputB, window, app, job_name = timestamp.strftime("SATELLiTES_%d%b%Y_%H%M%S"))
    else:
        CB_run_3reagents(reagentsA, reagentsB, reagentsC, reaction, output_path, nb_cpu, modeA, modeB, modeC, inputA, inputB, inputC, step, step_path, window, app, job_name = timestamp.strftime("SATELLiTES_%d%b%Y_%H%M%S"))

##################################################################################################################################

class Window(QMainWindow):
    def update_test(self):
        self.launch_button.hide()
        # self.layout.addWidget(self.progress_bar)
        self.layout.addWidget(self.progress_text)
        
    def update_text(self, value, start = 0):
        if start == 1:
            self.progress_text.setHtml(value)
        else:
            temp = self.progress_text.toHtml() + value
            self.progress_text.setHtml(temp)
        self.progress_text.moveCursor(QTextCursor.MoveOperation.End)
    
    def update_progress(self, value):
        # self.progress_bar.setValue(value)
        pass

    def __init__(self):
        super().__init__()
        
        global installdir_SATELLiTES

        app_icon = QIcon()
        app_icon.addFile(installdir_SATELLiTES + 'SATELLiTES_favicon.png', QSize(256,256))
        self.setWindowIcon(app_icon)
        self.setWindowTitle("SATELLiTES - Version 1.0.10")
        # self.setGeometry(100, 100, 800, 300)
        self.setGeometry(0,0,1000, 300)

        self.lineheight = 30
        self.main_widget = QWidget(self)
        self.setCentralWidget(self.main_widget)

        self.layout = QVBoxLayout(self.main_widget)

        #########################

        self.layout_header = QHBoxLayout()

        #####

        self.left_column = QVBoxLayout()
        self.left_column_widget = QWidget()  # Create a widget to hold the layout
        self.left_column_widget.setLayout(self.left_column)
        self.left_column_widget.setFixedWidth(200)  # Fixed width for left column
        self.layout_header.addWidget(self.left_column_widget)

        #####

        self.label_logo = QLabel(self)
        self.label_logo.setPixmap(QPixmap(installdir_SATELLiTES + "SATELLiTES_header.png"))
        self.layout_header.addWidget(self.label_logo, alignment=Qt.AlignmentFlag.AlignCenter)
        
        #####

        self.right_column = QVBoxLayout()
        self.button_help = QPushButton('Wiki help')
        self.button_github = QPushButton('Github')
        self.button_cite = QPushButton('Cite us')
        self.right_column.addWidget(self.button_help)
        self.right_column.addWidget(self.button_github)
        self.right_column.addWidget(self.button_cite)
        self.button_help.clicked.connect(self.open_help)
        self.button_github.clicked.connect(self.open_github)
        self.button_cite.clicked.connect(self.open_cite)
        self.right_column_widget = QWidget()  # Create a widget to hold the layout
        self.right_column_widget.setLayout(self.right_column)
        self.right_column_widget.setFixedWidth(200)  # Fixed width for right column
        self.layout_header.addWidget(self.right_column_widget)

        #####

        self.layout.addLayout(self.layout_header)

        #########################

        self.layout_top = QVBoxLayout()
        self.layout.addLayout(self.layout_top)

        #####

        row3_layout = QHBoxLayout()
        self.layout_top.addLayout(row3_layout)

        self.label4 = QLabel("SMARTS chemical reaction")
        self.label4.setFixedHeight(self.lineheight)
        row3_layout.addWidget(self.label4)

        self.smartsReaction = QLineEdit(self)
        self.smartsReaction.setPlaceholderText("Input SMARTS chemical reactions in the form A.B>>P")
        row3_layout.addWidget(self.smartsReaction)

        #####

        row4_layout = QHBoxLayout()
        self.layout_top.addLayout(row4_layout)

        self.outputdir_label = QLabel("Output directory")
        self.outputdir_label.setFixedHeight(self.lineheight)
        row4_layout.addWidget(self.outputdir_label)

        self.outputdir_textbox = QLineEdit(self)
        self.outputdir_textbox.setPlaceholderText("Selected path for the output")
        row4_layout.addWidget(self.outputdir_textbox)

        self.outputdir_button = QPushButton("Browse", self)
        self.outputdir_button.clicked.connect(self.select_outputdir)
        row4_layout.addWidget(self.outputdir_button)

        #####

        row5_layout = QHBoxLayout()
        self.layout_top.addLayout(row5_layout)

        self.cpu_label = QLabel()
        self.cpu_label.setFixedHeight(self.lineheight)
        row5_layout.addWidget(self.cpu_label)

        self.cpu_slider = QSlider()
        self.cpu_slider.setMinimum(1)
        self.cpu_slider.setMaximum(mp.cpu_count())
        self.cpu_slider.setValue(round(mp.cpu_count() / 2))
        self.cpu_label.setText(f"Logical processors to allocate: {self.cpu_slider.value()}")
        self.cpu_slider.setOrientation(Qt.Orientation.Horizontal)
        self.cpu_slider.valueChanged.connect(self.update_cpu_label)
        row5_layout.addWidget(self.cpu_slider)
        
        #####
        # TYPE OF CHEMICAL REACTION

        self.reaction_layout = QHBoxLayout()
        self.layout.addLayout(self.reaction_layout, stretch=1)
        self.setStyleSheet('''
        QTabWidget::tab-bar {
            alignment: center;
        }''')
        self.reaction_tabs = QTabWidget()
        self.reaction_layout.addWidget(self.reaction_tabs)

        self.reaction_2reagents = QWidget()
        self.reaction_3reagents = QWidget()

        self.reaction_tabs.setStyleSheet("QTabBar::tab { color:black; height: 60px; width: 400px;}")

        self.reaction_tabs.addTab(self.reaction_2reagents, "2-reagents reaction")
        self.reaction_tabs.addTab(self.reaction_3reagents, "3-reagents reaction")

        #####
        # TABS - 2 REAGENTS
        tabs_2reagents(self)

        # TABS - 3 REAGENTS
        tabs_3reagents(self)

        #####
        # LAUNCH
        self.launch_button = QPushButton("Launch", self)
        self.launch_button.clicked.connect(self.launch_function_handler)
        self.layout.addWidget(self.launch_button)

        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)

        self.progress_text = QTextEdit()
        self.progress_text.setReadOnly(True)
        self.progress_text.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self.progress_text.setFixedHeight(10 * self.lineheight)

    #####
    # FUNCTIONS
    def open_help(self):
        url = QUrl('https://github.com/cbedart/SATELLiTES/issues')
        QDesktopServices.openUrl(url)

    def open_github(self):
        url = QUrl('https://github.com/cbedart/SATELLiTES')
        QDesktopServices.openUrl(url)

    def open_cite(self):
        url = QUrl('https://github.com/cbedart/SATELLiTES')
        QDesktopServices.openUrl(url)

    def getVariable(self):
        return

    ##### 2-REAGENTS #####

    def change_tabA_2R(self, index):
        if index == 3: # If FEL, change the other tab to FEL too
            self.tabsA_2R.setCurrentIndex(index)
            self.tabsB_2R.setCurrentIndex(index)
            self.dropdownB_2R.setCurrentIndex(index)
        else:
            if self.dropdownB_2R.currentIndex() == 3: # If it was FEL before, change the other tab automatically
                self.tabsB_2R.setCurrentIndex(index)
                self.dropdownB_2R.setCurrentIndex(index)
            self.tabsA_2R.setCurrentIndex(index)

    def change_tabB_2R(self, index):
        if index == 3: # If FEL, change the other tab to FEL too
            self.tabsB_2R.setCurrentIndex(index)
            self.tabsA_2R.setCurrentIndex(index)
            self.dropdownA_2R.setCurrentIndex(index)
        else:
            if self.dropdownA_2R.currentIndex() == 3: # If it was FEL before, change the other tab automatically
                self.tabsA_2R.setCurrentIndex(index)
                self.dropdownA_2R.setCurrentIndex(index)
            self.tabsB_2R.setCurrentIndex(index)

    #####

    def open_file_reagentA_2R(self):
        file_dialog = QFileDialog(self, "Select reagents A SMILES file")
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        file_dialog.setNameFilter("SMI Files (*.smi)")

        if file_dialog.exec():
            file_path_reagentA = file_dialog.selectedFiles()[0]
            self.path_textbox_A_2R.setText(file_path_reagentA)

    def open_file_reagentB_2R(self):
        file_dialog = QFileDialog(self, "Select reagents B SMILES file")
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        file_dialog.setNameFilter("SMI Files (*.smi)")

        if file_dialog.exec():
            file_path_reagentB = file_dialog.selectedFiles()[0]
            self.path_textbox_B_2R.setText(file_path_reagentB)



    ##### 3-REAGENTS #####
    def open_file_reagentA_3R(self):
        file_dialog = QFileDialog(self, "Select reagents A SMILES file")
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        file_dialog.setNameFilter("SMI Files (*.smi)")

        if file_dialog.exec():
            file_path_reagentA = file_dialog.selectedFiles()[0]
            self.path_textbox_A_3R.setText(file_path_reagentA)

    def open_file_reagentB_3R(self):
        file_dialog = QFileDialog(self, "Select reagents B SMILES file")
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        file_dialog.setNameFilter("SMI Files (*.smi)")

        if file_dialog.exec():
            file_path_reagentB = file_dialog.selectedFiles()[0]
            self.path_textbox_B_3R.setText(file_path_reagentB)

    def open_file_reagentC_3R(self):
        file_dialog = QFileDialog(self, "Select reagents C SMILES file")
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        file_dialog.setNameFilter("SMI Files (*.smi)")

        if file_dialog.exec():
            file_path_reagentC = file_dialog.selectedFiles()[0]
            self.path_textbox_C_3R.setText(file_path_reagentC)

    def update_cpu_label(self, value):
        self.cpu_label.setText(f"Logical processors to allocate: {value}")

    #########################

    def error_management(self, element, var_to_check):
        if var_to_check == "":
            element.setStyleSheet("QLineEdit { border: 2px solid red; }")
            return 1
        else:
            element.setStyleSheet("QLineEdit { }")
            return 0
        
    #########################



    def launch_function_handler(self):
        reaction = self.smartsReaction.text()
        output_path = self.outputdir_textbox.text()
        nb_cpu = self.cpu_slider.value()

        # 0 = 2-reagents // 1 = 3-reagents // +2 to match the number of reagents
        reagents_nb = self.reaction_tabs.currentIndex() + 2
        
        ##### 2-reagents #####
        if reagents_nb == 2:
            reagentsA = self.path_textbox_A_2R.text()
            reagentsB = self.path_textbox_B_2R.text()

            reagentsC = ""
            modeC = ""
            inputC = ""
            step_path = ""

            if not self.frame_step1_2R.isVisible() and not self.frame_step2_2R.isVisible():
                print("Nothing selected")
                return ""
            # STEP #1
            if self.frame_step1_2R.isVisible():
                step = 1

                # Input - Reagents A representative
                if self.frame1_A2R.isVisible():
                    modeA = 0
                    inputA = ""
                if self.frame2_A2R.isVisible():
                    modeA = 1
                    inputA = self.representative_A2R.text()
                if self.frame3_A2R.isVisible():
                    modeA = 2
                    inputA = self.smiles_A2R.text()

                # Input - Reagents B representative
                if self.frame1_B2R.isVisible():
                    modeB = 0
                    inputB = ""
                if self.frame2_B2R.isVisible():
                    modeB = 1
                    inputB = self.representative_B2R.text()
                if self.frame3_B2R.isVisible():
                    modeB = 2
                    inputB = self.smiles_B2R.text()

            # STEP #1
            if self.frame_step2_2R.isVisible():
                step = 2
                modeA = 3
                inputA = self.selected_step2_2R.text()
                modeB = 3
                inputB = self.selected_step2_2R.text()
                step_path = self.selected_step2_2R.text()

        ##### 3-reagents #####
        elif reagents_nb == 3:
            reagentsA = self.path_textbox_A_3R.text()
            reagentsB = self.path_textbox_B_3R.text()
            reagentsC = self.path_textbox_C_3R.text()

            if not self.frame_step1_3R.isVisible() and not self.frame_step2_3R.isVisible() and not self.frame_step3_3R.isVisible():
                return "ERROR - Impossible scenario"
            
            # STEP #1
            if self.frame_step1_3R.isVisible():
                step = 1
                step_path = ""
                # Input - Reagents A representative
                if self.frame1_step1_A3R.isVisible():
                    modeA = 0
                    inputA = ""
                if self.frame2_step1_A3R.isVisible():
                    modeA = 1
                    inputA = self.representative_step1_A3R.text()
                if self.frame3_step1_A3R.isVisible():
                    modeA = 2
                    inputA = self.smiles_step1_A3R.text()

                # Input - Reagents B representative
                if self.frame1_step1_B3R.isVisible():
                    modeB = 0
                    inputB = ""
                if self.frame2_step1_B3R.isVisible():
                    modeB = 1
                    inputB = self.representative_step1_B3R.text()
                if self.frame3_step1_B3R.isVisible():
                    modeB = 2
                    inputB = self.smiles_step1_B3R.text()
                
                # Input - Reagents C representative
                if self.frame1_step1_C3R.isVisible():
                    modeC = 0
                    inputC = ""
                if self.frame2_step1_C3R.isVisible():
                    modeC = 1
                    inputC = self.representative_step1_C3R.text()
                if self.frame3_step1_C3R.isVisible():
                    modeC = 2
                    inputC = self.smiles_step1_C3R.text()

            # STEP #2
            if self.frame_step2_3R.isVisible():
                step = 2
                step_path = self.selected_step2_3R.text()

                # Input - Reagents A representative
                if self.frame1_step2_A3R.isVisible():
                    modeA = 0
                    inputA = ""
                if self.frame2_step2_A3R.isVisible():
                    modeA = 1
                    inputA = self.representative_step2_A3R.text()
                if self.frame3_step2_A3R.isVisible():
                    modeA = 2
                    inputA = self.smiles_step2_A3R.text()

                # Input - Reagents B representative
                if self.frame1_step2_B3R.isVisible():
                    modeB = 0
                    inputB = ""
                if self.frame2_step2_B3R.isVisible():
                    modeB = 1
                    inputB = self.representative_step2_B3R.text()
                if self.frame3_step2_B3R.isVisible():
                    modeB = 2
                    inputB = self.smiles_step2_B3R.text()
                
                # Input - Reagents C representative
                if self.frame1_step2_C3R.isVisible():
                    modeC = 0
                    inputC = ""
                if self.frame2_step2_C3R.isVisible():
                    modeC = 1
                    inputC = self.representative_step2_C3R.text()
                if self.frame3_step2_C3R.isVisible():
                    modeC = 2
                    inputC = self.smiles_step2_C3R.text()

            # STEP #3
            if self.frame_step3_3R.isVisible():
                step = 3
                step_path = self.selected_step3_3R.text()
                modeA = 3 ; modeB = 3 ; modeC = 3
                inputA = "" ; inputB = "" ; inputC = ""
                
        else:
            return "ERROR - Impossible scenario"

        # Errors handler
        errors = 0

        errors += self.error_management(self.smartsReaction, reaction)
        errors += self.error_management(self.outputdir_textbox, output_path)
        if not os.path.exists(output_path):
            errors += self.error_management(self.outputdir_textbox, "")

        # Errors - 2 reagents
        if reagents_nb == 2:
            errors += self.error_management(self.path_textbox_A_2R, reagentsA)
            errors += self.error_management(self.path_textbox_B_2R, reagentsB)

            if step == 1:
                if modeA == 1:
                    errors += self.error_management(self.representative_A2R, inputA)
                if modeB == 1:
                    errors += self.error_management(self.representative_B2R, inputB)

                if modeA == 2:
                    errors += self.error_management(self.smiles_A2R, inputA)
                if modeB == 2:
                    errors += self.error_management(self.smiles_B2R, inputB)

            if step == 2:
                errors += self.error_management(self.selected_step2_2R, step_path)
                if not os.path.exists(step_path):
                    errors += self.error_management(self.selected_step2_2R, "")

        # Errors - 3 reagents
        if reagents_nb == 3:
            errors += self.error_management(self.path_textbox_A_3R, reagentsA)
            errors += self.error_management(self.path_textbox_B_3R, reagentsB)
            errors += self.error_management(self.path_textbox_C_3R, reagentsC)

            if step == 1:
                if modeA == 1:
                    errors += self.error_management(self.representative_step1_A3R, inputA)
                if modeB == 1:
                    errors += self.error_management(self.representative_step1_B3R, inputB)
                if modeC == 1:
                    errors += self.error_management(self.representative_step1_C3R, inputC)

                if modeA == 2:
                    errors += self.error_management(self.smiles_step1_A3R, inputA)
                if modeB == 2:
                    errors += self.error_management(self.smiles_step1_B3R, inputB)
                if modeC == 2:
                    errors += self.error_management(self.smiles_step1_C3R, inputC)

            if step == 2:
                errors += self.error_management(self.selected_step2_3R, step_path)
                if not os.path.exists(step_path):
                    errors += self.error_management(self.selected_step2_3R, "")

                if modeA == 1:
                    errors += self.error_management(self.representative_step2_A3R, inputA)
                if modeB == 1:
                    errors += self.error_management(self.representative_step2_B3R, inputB)
                if modeC == 1:
                    errors += self.error_management(self.representative_step2_C3R, inputC)

                if modeA == 2:
                    errors += self.error_management(self.smiles_step2_A3R, inputA)
                if modeB == 2:
                    errors += self.error_management(self.smiles_step2_B3R, inputB)
                if modeC == 2:
                    errors += self.error_management(self.smiles_step2_C3R, inputC)


            if step == 3:
                errors += self.error_management(self.selected_step3_3R, step_path)
                if not os.path.exists(step_path):
                    errors += self.error_management(self.selected_step3_3R, "")

        if errors > 0:
            return "ERROR"

        # Lock line edits
        self.smartsReaction.setEnabled(False)
        self.outputdir_textbox.setEnabled(False)
        self.outputdir_button.setEnabled(False)
        self.cpu_slider.setEnabled(False)
        self.reaction_tabs.setEnabled(False)

        # Launch main function
        launch_function(reagentsA, reagentsB, reagentsC, reaction, output_path, nb_cpu, modeA, modeB, modeC, inputA, inputB, inputC, step, step_path)
    
    def select_outputdir(self):
        directory = QFileDialog.getExistingDirectory(self, "Select output directory")
        self.outputdir_textbox.setText(directory)


##################################################################################################################################

def run_SATELLiTES():
    global app
    global installdir_SATELLiTES

    app = QApplication(sys.argv)
    qqq = QIcon(installdir_SATELLiTES + 'SATELLiTES_favicon.png')
    app.setWindowIcon(qqq)

    sys.path.insert(0, sys.argv[0])
    # print(sys.path)

    app.setStyle('Fusion')

    with open(installdir_SATELLiTES + "SATELLiTES_style.css", "r") as f:
        style = f.read()
    
    style_bgs = "QMainWindow {background-image: url(" + installdir_SATELLiTES + "SATELLiTES_bg01.png);} QTabWidget>QWidget>QWidget {background-image: url(" + installdir_SATELLiTES + "SATELLiTES_bg02.png);}\n"

    app.setStyleSheet(style_bgs + style)

    global window
    window = Window()
    window.show()

    sys.exit(app.exec())

##################################################################################################################################


##################################################################################################################################
