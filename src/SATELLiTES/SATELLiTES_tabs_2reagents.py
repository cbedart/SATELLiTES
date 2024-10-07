####################################################################################################
# SATELLiTES - version 1.0.10
# Copyright (C) 2024 Corentin BEDART
####################################################################################################

import sys
import os

from PyQt6.QtWidgets import *
from PyQt6.QtCore import Qt, QTimer, QSize
from PyQt6.QtGui import QFont, QPixmap, QImage, QTextCursor, QIcon

import time

from rdkit import Chem
from rdkit.Chem import Draw, Descriptors

import multiprocessing as mp

##################################################################################################################################

def CB_test(self):
    self.frame_step1_2R.hide()
    self.frame_step2_2R.show()
    self.button_step1_2R.setStyleSheet("background-color:none;")
    self.button_step2_2R.setStyleSheet("background-color:lightgray;")
    self.adjustSize()

def CB_adjust(self, element):
    self.resize(self.width(), element.sizeHint().height())

def generate_image_A2R(self):
    smiles = self.smiles_A2R.text()
    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        image = Draw.MolToImage(molecule, size=(300,100))
        image_data = image.convert("RGBA").tobytes("raw", "RGBA")
        q_image = QImage(image_data, image.size[0], image.size[1], QImage.Format.Format_RGBA8888_Premultiplied)

        pixmap = QPixmap.fromImage(q_image)
        self.image_A2R.setPixmap(pixmap)

        self.mw_A2R.setText(f"Molecular weight = {round(Descriptors.ExactMolWt(molecule), 1)} Da")
    else:
        self.image_A2R.setText("Wrong SMILES")
        self.mw_A2R.clear()

def refresh_image_A2R(self):
    generate_image_A2R(self)

def schedule_refresh_image_A2R(self):
    self.timer_A2R.start(1000)  # Start the timer with a 1-second delay

#####

def generate_image_B2R(self):
    smiles = self.smiles_B2R.text()
    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        image = Draw.MolToImage(molecule, size=(300,100))
        image_data = image.convert("RGBA").tobytes("raw", "RGBA")
        q_image = QImage(image_data, image.size[0], image.size[1], QImage.Format.Format_RGBA8888_Premultiplied)

        pixmap = QPixmap.fromImage(q_image)
        self.image_B2R.setPixmap(pixmap)

        self.mw_B2R.setText(f"Molecular weight = {round(Descriptors.ExactMolWt(molecule), 1)} Da")
    else:
        self.image_B2R.setText("Wrong SMILES")
        self.mw_B2R.clear()

def refresh_image_B2R(self):
    generate_image_B2R(self)

def schedule_refresh_image_B2R(self):
    self.timer_B2R.start(1000)  # Start the timer with a 1-second delay

#####

def open_file_reagent_step2_2R(self):
    file_dialog = QFileDialog(self, "Select your SMILES file containing your selected compounds")
    file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
    file_dialog.setNameFilter("SMI Files (*.smi)")

    if file_dialog.exec():
        file_path_reagent = file_dialog.selectedFiles()[0]
        self.selected_step2_2R.setText(file_path_reagent)

##################################################################################################################################

def tabs_2reagents(self):

    self.reaction_2reagents_inside = QVBoxLayout()
    self.reaction_2reagents_inside.setAlignment(Qt.AlignmentFlag.AlignTop)
    self.reaction_2reagents.setLayout(self.reaction_2reagents_inside)

    row_inputA_2R_layout = QHBoxLayout()
    self.reaction_2reagents_inside.addLayout(row_inputA_2R_layout)

    self.label_reagentsA_2R = QLabel("Reagents A - Path")
    self.label_reagentsA_2R.setFixedHeight(self.lineheight)
    row_inputA_2R_layout.addWidget(self.label_reagentsA_2R)

    self.path_textbox_A_2R = QLineEdit(self)
    self.path_textbox_A_2R.setPlaceholderText("File path for reagents A")
    row_inputA_2R_layout.addWidget(self.path_textbox_A_2R)

    self.browse_button_A_2R = QPushButton("Browse", self)
    self.browse_button_A_2R.clicked.connect(self.open_file_reagentA_2R)
    row_inputA_2R_layout.addWidget(self.browse_button_A_2R)

    #####

    row_inputB_2R_layout = QHBoxLayout()
    self.reaction_2reagents_inside.addLayout(row_inputB_2R_layout)

    self.label_reagentsB_2R = QLabel("Reagents B - Path")
    self.label_reagentsB_2R.setFixedHeight(self.lineheight)
    row_inputB_2R_layout.addWidget(self.label_reagentsB_2R)

    self.path_textbox_B_2R = QLineEdit(self)
    self.path_textbox_B_2R.setPlaceholderText("File path for reagents B")
    row_inputB_2R_layout.addWidget(self.path_textbox_B_2R)

    self.browse_button_B_2R = QPushButton("Browse", self)
    self.browse_button_B_2R.clicked.connect(self.open_file_reagentB_2R)
    row_inputB_2R_layout.addWidget(self.browse_button_B_2R)


    ##################################################################################################################################

    self.button_step1_2R = QPushButton("Step #1 - Basis Products Library")
    self.button_step2_2R = QPushButton("Step #2 - Focused Products Library")

    self.button_step1_2R.clicked.connect(lambda: (self.frame_step1_2R.show(), self.frame_step2_2R.hide(),
                                                  self.button_step1_2R.setStyleSheet("background-color:lightgray;"),
                                                  self.button_step2_2R.setStyleSheet("background-color:none;")
                                                  ))

    # self.button_step2_2R.clicked.connect(lambda: (self.frame_step1_2R.hide(), self.frame_step2_2R.show(),
    #                                               self.button_step1_2R.setStyleSheet("background-color:none;"),
    #                                               self.button_step2_2R.setStyleSheet("background-color:lightgray;"),
    #                                               self.resize(self.width(), self.frame_step2_2R.sizeHint().height())
    #                                               ))
    self.button_step2_2R.clicked.connect(lambda: (self.frame_step1_2R.hide(), self.frame_step2_2R.show(),
                                                  self.button_step1_2R.setStyleSheet("background-color:none;"),
                                                  self.button_step2_2R.setStyleSheet("background-color:lightgray;")
                                                  ))
    self.button_step2_2R.clicked.connect(lambda:(self.resize(self.width(), self.frame_step2_2R.sizeHint().height()), self.layout.update()))

    self.buttons_layout_2R = QHBoxLayout()
    self.buttons_layout_2R.addWidget(self.button_step1_2R)
    self.buttons_layout_2R.addWidget(self.button_step2_2R)
    self.reaction_2reagents_inside.addLayout(self.buttons_layout_2R)

    ###########

    self.frame_step1_2R = QFrame()
    self.layout_tabs_step1_2R = QHBoxLayout()
    self.frame_step1_2R.setLayout(self.layout_tabs_step1_2R)
    self.reaction_2reagents_inside.addWidget(self.frame_step1_2R)

    ###########

    self.frame_step2_2R = QFrame()
    self.layout_tabs_step2_2R = QHBoxLayout()
    self.frame_step2_2R.setLayout(self.layout_tabs_step2_2R)
    self.reaction_2reagents_inside.addWidget(self.frame_step2_2R)

    ###########    

    self.frame_step1_2R.hide()
    self.frame_step2_2R.hide()

    
    ##################################################################################################################################
    # STEP #1 - MEL
    ##################################################################################################################################

    # TABS - REAGENT A

    self.layout_tabsA_2R = QVBoxLayout()
    self.layout_tabs_step1_2R.addLayout(self.layout_tabsA_2R, 50)

    self.tabsA_2R_title = QLabel("Reagent A representative")
    self.tabsA_2R_title.setToolTip("The representative of the full reagents A dataset, named a, to generate Ba chemical libraries")
    # self.tabsA_2R_title.setToolTipDuration(100)
    self.tabsA_2R_title.setFixedHeight(self.lineheight)
    self.tabsA_2R_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
    self.tabsA_2R_title_font = QFont("Calibri", 10)
    self.tabsA_2R_title_font.setBold(True)
    self.tabsA_2R_title.setFont(self.tabsA_2R_title_font)
    self.layout_tabsA_2R.addWidget(self.tabsA_2R_title)

    ###########

    self.button1_A2R = QPushButton("Smallest")
    self.button2_A2R = QPushButton("Chosen ID")
    self.button3_A2R = QPushButton("Custom SMILES")

    self.button1_A2R.clicked.connect(lambda: (self.frame1_A2R.show(), self.frame2_A2R.hide(), self.frame3_A2R.hide(), 
                                              self.button1_A2R.setStyleSheet("background-color:lightgrey;"),
                                              self.button2_A2R.setStyleSheet("background-color:none;"),
                                              self.button3_A2R.setStyleSheet("background-color:none;"),
                                              ))
    
    self.button2_A2R.clicked.connect(lambda: (self.frame1_A2R.hide(), self.frame2_A2R.show(), self.frame3_A2R.hide(), 
                                              self.button1_A2R.setStyleSheet("background-color:none;"),
                                              self.button2_A2R.setStyleSheet("background-color:lightgrey;"),
                                              self.button3_A2R.setStyleSheet("background-color:none;"),
                                              ))

    self.button3_A2R.clicked.connect(lambda: (self.frame1_A2R.hide(), self.frame2_A2R.hide(), self.frame3_A2R.show(), 
                                              self.button1_A2R.setStyleSheet("background-color:none;"),
                                              self.button2_A2R.setStyleSheet("background-color:none;"),
                                              self.button3_A2R.setStyleSheet("background-color:lightgrey;"),
                                              ))
    
    self.buttons_A2R_layout = QHBoxLayout()
    self.buttons_A2R_layout.addWidget(self.button1_A2R)
    self.buttons_A2R_layout.addWidget(self.button2_A2R)
    self.buttons_A2R_layout.addWidget(self.button3_A2R)
    self.layout_tabsA_2R.addLayout(self.buttons_A2R_layout)

    ###########

    self.frame1_A2R = QFrame()
    # self.frame1_A2R.setFixedHeight(150)
    self.layout1_A2R = QVBoxLayout()

    self.frame1_A2R.setLayout(self.layout1_A2R)
    self.layout_tabsA_2R.addWidget(self.frame1_A2R)
    self.label1_A2R = QLabel("No parameters required")
    self.label1_A2R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter)
    self.layout1_A2R.addWidget(self.label1_A2R)
    
    ###########

    self.frame2_A2R = QFrame()
    # self.frame2_A2R.setFixedHeight(150)
    self.layout2_A2R = QVBoxLayout()
    self.layout2_A2R.setAlignment(Qt.AlignmentFlag.AlignVCenter)
    self.frame2_A2R.setLayout(self.layout2_A2R)
    self.layout_tabsA_2R.addWidget(self.frame2_A2R)

    self.label2_A2R = QLabel("Representative SMILES ID")
    self.label2_A2R.setAlignment(Qt.AlignmentFlag.AlignHCenter)
    self.layout2_A2R.addWidget(self.label2_A2R)

    self.representative_A2R = QLineEdit()
    self.representative_A2R.setPlaceholderText("Existing SMILES ID in reagents A file")
    self.layout2_A2R.addWidget(self.representative_A2R)

    # def lock(self):
    #     self.representative_ID.setEnabled(False)

    ###########

    self.frame3_A2R = QFrame()
    # self.frame3_A2R.setFixedHeight(150)
    self.layout3_A2R = QVBoxLayout()
    self.layout3_A2R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.frame3_A2R.setLayout(self.layout3_A2R)
    self.layout_tabsA_2R.addWidget(self.frame3_A2R)

    self.smiles_A2R = QLineEdit()
    self.smiles_A2R.setPlaceholderText("SMILES string for reagent A representative")
    
    self.timer_A2R = QTimer()
    self.timer_A2R.setSingleShot(True)
    self.layout3_A2R.addWidget(self.smiles_A2R)

    self.image_A2R = QLabel()
    self.image_A2R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_A2R.addWidget(self.image_A2R)

    self.mw_A2R = QLabel()
    self.mw_A2R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_A2R.addWidget(self.mw_A2R)
    self.smiles_A2R.textChanged.connect(lambda:schedule_refresh_image_A2R(self))
    self.timer_A2R.timeout.connect(lambda:refresh_image_A2R(self))

    ###########

    self.frame1_A2R.show()
    self.frame2_A2R.hide()
    self.frame3_A2R.hide()
    self.button1_A2R.setStyleSheet("background-color:lightgrey;")

    ##################################################################################################################################
    # SEPARATOR
    self.separator = QFrame()
    self.separator.setFrameShape(QFrame.Shape.VLine)
    # self.separator.setFrameShadow(QFrame.Sunken)
    self.separator.setLineWidth(2)
    self.separator.setStyleSheet("color: gray;")
    self.layout_tabs_step1_2R.addWidget(self.separator)

    ##################################################################################################################################
    # TABS - REAGENT B
    self.layout_tabsB_2R = QVBoxLayout()
    self.layout_tabs_step1_2R.addLayout(self.layout_tabsB_2R, 50)

    self.tabsB_2R_title = QLabel("Reagent B representative")
    self.tabsB_2R_title.setToolTip("The representative of the full reagents B dataset, named b, to generate Ab chemical libraries")
    self.tabsB_2R_title.setFixedHeight(self.lineheight)
    self.tabsB_2R_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
    self.tabsB_2R_title_font = QFont("Calibri", 10)
    self.tabsB_2R_title_font.setBold(True)
    self.tabsB_2R_title.setFont(self.tabsB_2R_title_font)
    self.layout_tabsB_2R.addWidget(self.tabsB_2R_title)

    ###########

    self.button1_B2R = QPushButton("Smallest")
    self.button2_B2R = QPushButton("Chosen ID")
    self.button3_B2R = QPushButton("Custom SMILES")
    self.button1_B2R.clicked.connect(lambda: (self.frame1_B2R.show(), self.frame2_B2R.hide(), self.frame3_B2R.hide(), 
                                              self.button1_B2R.setStyleSheet("background-color:lightgrey;"),
                                              self.button2_B2R.setStyleSheet("background-color:none;"),
                                              self.button3_B2R.setStyleSheet("background-color:none;"),
                                              ))


    self.button2_B2R.clicked.connect(lambda: (self.frame1_B2R.hide(), self.frame2_B2R.show(), self.frame3_B2R.hide(), 
                                              self.button1_B2R.setStyleSheet("background-color:none;"),
                                              self.button2_B2R.setStyleSheet("background-color:lightgrey;"),
                                              self.button3_B2R.setStyleSheet("background-color:none;"),
                                              ))


    self.button3_B2R.clicked.connect(lambda: (self.frame1_B2R.hide(), self.frame2_B2R.hide(), self.frame3_B2R.show(), 
                                              self.button1_B2R.setStyleSheet("background-color:none;"),
                                              self.button2_B2R.setStyleSheet("background-color:none;"),
                                              self.button3_B2R.setStyleSheet("background-color:lightgrey;"),
                                              ))


    self.buttons_B2R_layout = QHBoxLayout()
    self.buttons_B2R_layout.addWidget(self.button1_B2R)
    self.buttons_B2R_layout.addWidget(self.button2_B2R)
    self.buttons_B2R_layout.addWidget(self.button3_B2R)

    self.layout_tabsB_2R.addLayout(self.buttons_B2R_layout)
    
    ###########

    self.frame1_B2R = QFrame()
    self.layout1_B2R = QVBoxLayout()
    self.frame1_B2R.setLayout(self.layout1_B2R)
    self.layout_tabsB_2R.addWidget(self.frame1_B2R)
    self.label1_B2R = QLabel("No parameters required")
    self.label1_B2R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter)
    self.layout1_B2R.addWidget(self.label1_B2R)
    
    ###########

    self.frame2_B2R = QFrame()
    self.layout2_B2R = QVBoxLayout()
    self.layout2_B2R.setAlignment(Qt.AlignmentFlag.AlignVCenter)
    self.frame2_B2R.setLayout(self.layout2_B2R)
    self.layout_tabsB_2R.addWidget(self.frame2_B2R)

    self.label2_B2R = QLabel("Representative SMILES ID")
    self.label2_B2R.setAlignment(Qt.AlignmentFlag.AlignHCenter)
    self.layout2_B2R.addWidget(self.label2_B2R)

    self.representative_B2R = QLineEdit()
    self.representative_B2R.setPlaceholderText("Existing SMILES ID in reagents A file")
    self.layout2_B2R.addWidget(self.representative_B2R)

    # def lock(self):
    #     self.representative_ID.setEnabled(False)

    ###########

    self.frame3_B2R = QFrame()
    self.layout3_B2R = QVBoxLayout()
    self.layout3_B2R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.frame3_B2R.setLayout(self.layout3_B2R)
    self.layout_tabsB_2R.addWidget(self.frame3_B2R)

    self.smiles_B2R = QLineEdit()
    self.smiles_B2R.setPlaceholderText("SMILES string for reagent A representative")
    
    self.timer_B2R = QTimer()
    self.timer_B2R.setSingleShot(True)
    self.layout3_B2R.addWidget(self.smiles_B2R)

    self.image_B2R = QLabel()
    self.image_B2R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_B2R.addWidget(self.image_B2R)

    self.mw_B2R = QLabel()
    self.mw_B2R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_B2R.addWidget(self.mw_B2R)
    self.smiles_B2R.textChanged.connect(lambda:schedule_refresh_image_B2R(self))
    self.timer_B2R.timeout.connect(lambda:refresh_image_B2R(self))

    ###########

    self.frame1_B2R.show()
    self.frame2_B2R.hide()
    self.frame3_B2R.hide()
    self.button1_B2R.setStyleSheet("background-color:lightgrey;")

    ##################################################################################################################################
    # STEP #2 - FEL
    ##################################################################################################################################

    self.layout_step2_2R = QVBoxLayout()
    self.layout_tabs_step2_2R.addLayout(self.layout_step2_2R)

    self.label_step2 = QLabel("Selected reagents SMILES file to generate Focused Enumeration Libraries")
    self.label_step2.setAlignment(Qt.AlignmentFlag.AlignHCenter)
    self.layout_step2_2R.addWidget(self.label_step2)

    self.layout_step2_browse_2R = QHBoxLayout()
    self.layout_step2_2R.addLayout(self.layout_step2_browse_2R)

    self.selected_step2_2R = QLineEdit(self)
    self.selected_step2_2R.setPlaceholderText("File path for selected reagents")
    self.layout_step2_browse_2R.addWidget(self.selected_step2_2R)

    self.browse_step2_2R = QPushButton("Browse", self)
    self.browse_step2_2R.clicked.connect(lambda:open_file_reagent_step2_2R(self))
    self.layout_step2_browse_2R.addWidget(self.browse_step2_2R)
    
