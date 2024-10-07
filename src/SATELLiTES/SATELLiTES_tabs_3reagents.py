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

def generate_image_step1_A3R(self):
    smiles = self.smiles_step1_A3R.text()
    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        image = Draw.MolToImage(molecule, size=(300,100))
        image_data = image.convert("RGBA").tobytes("raw", "RGBA")
        q_image = QImage(image_data, image.size[0], image.size[1], QImage.Format.Format_RGBA8888_Premultiplied)

        pixmap = QPixmap.fromImage(q_image)
        self.image_step1_A3R.setPixmap(pixmap)

        self.mw_step1_A3R.setText(f"Molecular weight = {round(Descriptors.ExactMolWt(molecule), 1)} Da")
    else:
        self.image_step1_A3R.setText("Wrong SMILES")
        self.mw_step1_A3R.clear()

def refresh_image_step1_A3R(self):
    generate_image_step1_A3R(self)

def schedule_refresh_image_step1_A3R(self):
    self.timer_step1_A3R.start(1000)  # Start the timer with a 1-second delay

#####

def generate_image_step1_B3R(self):
    smiles = self.smiles_step1_B3R.text()
    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        image = Draw.MolToImage(molecule, size=(300,100))
        image_data = image.convert("RGBA").tobytes("raw", "RGBA")
        q_image = QImage(image_data, image.size[0], image.size[1], QImage.Format.Format_RGBA8888_Premultiplied)

        pixmap = QPixmap.fromImage(q_image)
        self.image_step1_B3R.setPixmap(pixmap)

        self.mw_step1_B3R.setText(f"Molecular weight = {round(Descriptors.ExactMolWt(molecule), 1)} Da")
    else:
        self.image_step1_B3R.setText("Wrong SMILES")
        self.mw_step1_B3R.clear()

def refresh_image_step1_B3R(self):
    generate_image_step1_B3R(self)

def schedule_refresh_image_step1_B3R(self):
    self.timer_step1_B3R.start(1000)  # Start the timer with a 1-second delay

#####

def generate_image_step1_C3R(self):
    smiles = self.smiles_step1_C3R.text()
    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        image = Draw.MolToImage(molecule, size=(300,100))
        image_data = image.convert("RGBA").tobytes("raw", "RGBA")
        q_image = QImage(image_data, image.size[0], image.size[1], QImage.Format.Format_RGBA8888_Premultiplied)

        pixmap = QPixmap.fromImage(q_image)
        self.image_step1_C3R.setPixmap(pixmap)

        self.mw_step1_C3R.setText(f"Molecular weight = {round(Descriptors.ExactMolWt(molecule), 1)} Da")
    else:
        self.image_step1_C3R.setText("Wrong SMILES")
        self.mw_step1_C3R.clear()

def refresh_image_step1_C3R(self):
    generate_image_step1_C3R(self)

def schedule_refresh_image_step1_C3R(self):
    self.timer_step1_C3R.start(1000)  # Start the timer with a 1-second delay

##################################################################################################################################

def generate_image_step2_A3R(self):
    smiles = self.smiles_step2_A3R.text()
    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        image = Draw.MolToImage(molecule, size=(300,100))
        image_data = image.convert("RGBA").tobytes("raw", "RGBA")
        q_image = QImage(image_data, image.size[0], image.size[1], QImage.Format.Format_RGBA8888_Premultiplied)

        pixmap = QPixmap.fromImage(q_image)
        self.image_step2_A3R.setPixmap(pixmap)

        self.mw_step2_A3R.setText(f"Molecular weight = {round(Descriptors.ExactMolWt(molecule), 1)} Da")
    else:
        self.image_step2_A3R.setText("Wrong SMILES")
        self.mw_step2_A3R.clear()

def refresh_image_step2_A3R(self):
    generate_image_step2_A3R(self)

def schedule_refresh_image_step2_A3R(self):
    self.timer_step2_A3R.start(1000)  # Start the timer with a 1-second delay

#####

def generate_image_step2_B3R(self):
    smiles = self.smiles_step2_B3R.text()
    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        image = Draw.MolToImage(molecule, size=(300,100))
        image_data = image.convert("RGBA").tobytes("raw", "RGBA")
        q_image = QImage(image_data, image.size[0], image.size[1], QImage.Format.Format_RGBA8888_Premultiplied)

        pixmap = QPixmap.fromImage(q_image)
        self.image_step2_B3R.setPixmap(pixmap)

        self.mw_step2_B3R.setText(f"Molecular weight = {round(Descriptors.ExactMolWt(molecule), 1)} Da")
    else:
        self.image_step2_B3R.setText("Wrong SMILES")
        self.mw_step2_B3R.clear()

def refresh_image_step2_B3R(self):
    generate_image_step2_B3R(self)

def schedule_refresh_image_step2_B3R(self):
    self.timer_step2_B3R.start(1000)  # Start the timer with a 1-second delay

#####

def generate_image_step2_C3R(self):
    smiles = self.smiles_step2_C3R.text()
    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        image = Draw.MolToImage(molecule, size=(300,100))
        image_data = image.convert("RGBA").tobytes("raw", "RGBA")
        q_image = QImage(image_data, image.size[0], image.size[1], QImage.Format.Format_RGBA8888_Premultiplied)

        pixmap = QPixmap.fromImage(q_image)
        self.image_step2_C3R.setPixmap(pixmap)

        self.mw_step2_C3R.setText(f"Molecular weight = {round(Descriptors.ExactMolWt(molecule), 1)} Da")
    else:
        self.image_step2_C3R.setText("Wrong SMILES")
        self.mw_step2_C3R.clear()

def refresh_image_step2_C3R(self):
    generate_image_step2_C3R(self)

def schedule_refresh_image_step2_C3R(self):
    self.timer_step2_C3R.start(1000)  # Start the timer with a 1-second delay

#####

def open_file_reagent_step2_3R(self):
    file_dialog = QFileDialog(self, "Select your SMILES file containing your selected compounds")
    file_dialog.setFileMode(QFileDialog.ExistingFile)
    file_dialog.setNameFilter("SMI Files (*.smi)")

    if file_dialog.exec_():
        file_path_reagent = file_dialog.selectedFiles()[0]
        self.selected_step2_3R.setText(file_path_reagent)
    
#####

def open_file_reagent_step3_3R(self):
    file_dialog = QFileDialog(self, "Select your SMILES file containing your selected compounds")
    file_dialog.setFileMode(QFileDialog.ExistingFile)
    file_dialog.setNameFilter("SMI Files (*.smi)")

    if file_dialog.exec_():
        file_path_reagent = file_dialog.selectedFiles()[0]
        self.selected_step3_3R.setText(file_path_reagent)


##################################################################################################################################

def tabs_3reagents(self):

    self.reaction_3reagents_inside = QVBoxLayout()
    self.reaction_3reagents_inside.setAlignment(Qt.AlignmentFlag.AlignTop)
    self.reaction_3reagents.setLayout(self.reaction_3reagents_inside)
    
    #####

    row_inputA_3R_layout = QHBoxLayout()
    self.reaction_3reagents_inside.addLayout(row_inputA_3R_layout)

    self.label_reagentsA_3R = QLabel("Reagents A - Path")
    self.label_reagentsA_3R.setFixedHeight(self.lineheight)
    row_inputA_3R_layout.addWidget(self.label_reagentsA_3R)

    self.path_textbox_A_3R = QLineEdit(self)
    self.path_textbox_A_3R.setPlaceholderText("File path for reagents A")
    row_inputA_3R_layout.addWidget(self.path_textbox_A_3R)

    self.browse_button_A_3R = QPushButton("Browse", self)
    self.browse_button_A_3R.clicked.connect(self.open_file_reagentA_3R)
    row_inputA_3R_layout.addWidget(self.browse_button_A_3R)

    #####

    row_inputB_3R_layout = QHBoxLayout()
    self.reaction_3reagents_inside.addLayout(row_inputB_3R_layout)

    self.label_reagentsB_3R = QLabel("Reagents B - Path")
    self.label_reagentsB_3R.setFixedHeight(self.lineheight)
    row_inputB_3R_layout.addWidget(self.label_reagentsB_3R)

    self.path_textbox_B_3R = QLineEdit(self)
    self.path_textbox_B_3R.setPlaceholderText("File path for reagents B")
    row_inputB_3R_layout.addWidget(self.path_textbox_B_3R)

    self.browse_button_B_3R = QPushButton("Browse", self)
    self.browse_button_B_3R.clicked.connect(self.open_file_reagentB_3R)
    row_inputB_3R_layout.addWidget(self.browse_button_B_3R)

    #####

    row_inputC_3R_layout = QHBoxLayout()
    self.reaction_3reagents_inside.addLayout(row_inputC_3R_layout)

    self.label_reagentsC_3R = QLabel("Reagents C - Path")
    self.label_reagentsC_3R.setFixedHeight(self.lineheight)
    row_inputC_3R_layout.addWidget(self.label_reagentsC_3R)

    self.path_textbox_C_3R = QLineEdit(self)
    self.path_textbox_C_3R.setPlaceholderText("File path for reagents C")
    row_inputC_3R_layout.addWidget(self.path_textbox_C_3R)

    self.browse_button_C_3R = QPushButton("Browse", self)
    self.browse_button_C_3R.clicked.connect(self.open_file_reagentC_3R)
    row_inputC_3R_layout.addWidget(self.browse_button_C_3R)

    ##################################################################################################################################
    
    self.button_step1_3R = QPushButton("Step #1 - Basis Products Library")
    self.button_step2_3R = QPushButton("Step #2 - Composite Products Library")
    self.button_step3_3R = QPushButton("Step #3 - Focused Products Library")

    self.button_step1_3R.clicked.connect(lambda: (self.frame_step1_3R.show(), self.frame_step2_3R.hide(), self.frame_step3_3R.hide(),
                                                  self.button_step1_3R.setStyleSheet("background-color:lightgray;"),
                                                  self.button_step2_3R.setStyleSheet("background-color:none;"),
                                                  self.button_step3_3R.setStyleSheet("background-color:none;")
                                                  ))
    self.button_step2_3R.clicked.connect(lambda: (self.frame_step1_3R.hide(), self.frame_step2_3R.show(), self.frame_step3_3R.hide(),
                                                  self.button_step1_3R.setStyleSheet("background-color:none;"),
                                                  self.button_step2_3R.setStyleSheet("background-color:lightgray;"),
                                                  self.button_step3_3R.setStyleSheet("background-color:none;")
                                                  ))
    self.button_step3_3R.clicked.connect(lambda: (self.frame_step1_3R.hide(), self.frame_step2_3R.hide(), self.frame_step3_3R.show(),
                                                  self.button_step1_3R.setStyleSheet("background-color:none;"),
                                                  self.button_step2_3R.setStyleSheet("background-color:none;"),
                                                  self.button_step3_3R.setStyleSheet("background-color:lightgray;")
                                                  ))

    self.buttons_layout_3R = QHBoxLayout()
    self.buttons_layout_3R.addWidget(self.button_step1_3R)
    self.buttons_layout_3R.addWidget(self.button_step2_3R)
    self.buttons_layout_3R.addWidget(self.button_step3_3R)
    self.reaction_3reagents_inside.addLayout(self.buttons_layout_3R)

    ###########

    self.frame_step1_3R = QFrame()
    self.layout_tabs_step1_3R = QHBoxLayout()
    self.frame_step1_3R.setLayout(self.layout_tabs_step1_3R)
    self.reaction_3reagents_inside.addWidget(self.frame_step1_3R)

    ###########

    self.frame_step2_3R = QFrame()
    self.layout_tabs_step2_3R_vertical = QVBoxLayout()
    self.frame_step2_3R.setLayout(self.layout_tabs_step2_3R_vertical)
    self.reaction_3reagents_inside.addWidget(self.frame_step2_3R)

    ###########

    self.frame_step3_3R = QFrame()
    self.layout_tabs_step3_3R = QVBoxLayout()
    self.frame_step3_3R.setLayout(self.layout_tabs_step3_3R)
    self.reaction_3reagents_inside.addWidget(self.frame_step3_3R)

    ###########

    self.frame_step1_3R.hide()
    self.frame_step2_3R.hide()
    self.frame_step3_3R.hide()

    ##################################################################################################################################
    # STEP #1 - MEL
    ##################################################################################################################################

    # TABS - REAGENT A

    self.layout_tabsA_step1_3R = QVBoxLayout()
    self.layout_tabs_step1_3R.addLayout(self.layout_tabsA_step1_3R, 50)

    self.tabsA_step1_3R_title = QLabel("Reagent A representative")
    self.tabsA_step1_3R_title.setToolTip("The representative of the full reagents A dataset, named a, to generate Ba chemical libraries")
    # self.tabsA_step1_3R_title.setToolTipDuration(100)
    self.tabsA_step1_3R_title.setFixedHeight(self.lineheight)
    self.tabsA_step1_3R_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
    self.tabsA_step1_3R_title_font = QFont("Calibri", 10)
    self.tabsA_step1_3R_title_font.setBold(True)
    self.tabsA_step1_3R_title.setFont(self.tabsA_step1_3R_title_font)
    self.layout_tabsA_step1_3R.addWidget(self.tabsA_step1_3R_title)

    ###########

    self.button1_step1_A3R = QPushButton("Smallest")
    self.button2_step1_A3R = QPushButton("Chosen ID")
    self.button3_step1_A3R = QPushButton("Custom SMILES")

    self.button1_step1_A3R.clicked.connect(lambda: (self.frame1_step1_A3R.show(), self.frame2_step1_A3R.hide(), self.frame3_step1_A3R.hide(), 
                                              self.button1_step1_A3R.setStyleSheet("background-color:lightgrey;"),
                                              self.button2_step1_A3R.setStyleSheet("background-color:none;"),
                                              self.button3_step1_A3R.setStyleSheet("background-color:none;"),
                                              ))

    self.button2_step1_A3R.clicked.connect(lambda: (self.frame1_step1_A3R.hide(), self.frame2_step1_A3R.show(), self.frame3_step1_A3R.hide(), 
                                              self.button1_step1_A3R.setStyleSheet("background-color:none;"),
                                              self.button2_step1_A3R.setStyleSheet("background-color:lightgrey;"),
                                              self.button3_step1_A3R.setStyleSheet("background-color:none;"),
                                              ))

    self.button3_step1_A3R.clicked.connect(lambda: (self.frame1_step1_A3R.hide(), self.frame2_step1_A3R.hide(), self.frame3_step1_A3R.show(), 
                                              self.button1_step1_A3R.setStyleSheet("background-color:none;"),
                                              self.button2_step1_A3R.setStyleSheet("background-color:none;"),
                                              self.button3_step1_A3R.setStyleSheet("background-color:lightgrey;"),
                                              ))
    
    self.buttons_step1_A3R_layout = QHBoxLayout()
    self.buttons_step1_A3R_layout.addWidget(self.button1_step1_A3R)
    self.buttons_step1_A3R_layout.addWidget(self.button2_step1_A3R)
    self.buttons_step1_A3R_layout.addWidget(self.button3_step1_A3R)
    self.layout_tabsA_step1_3R.addLayout(self.buttons_step1_A3R_layout)

    ###########

    self.frame1_step1_A3R = QFrame()
    # self.frame1_step1_A3R.setFixedHeight(150)
    self.layout1_step1_A3R = QVBoxLayout()

    self.frame1_step1_A3R.setLayout(self.layout1_step1_A3R)
    self.layout_tabsA_step1_3R.addWidget(self.frame1_step1_A3R)
    self.label1_step1_A3R = QLabel("No parameters required")
    self.label1_step1_A3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter)
    self.layout1_step1_A3R.addWidget(self.label1_step1_A3R)
    
    ###########

    self.frame2_step1_A3R = QFrame()
    # self.frame2_step1_A3R.setFixedHeight(150)
    self.layout2_step1_A3R = QVBoxLayout()
    self.layout2_step1_A3R.setAlignment(Qt.AlignmentFlag.AlignVCenter)
    self.frame2_step1_A3R.setLayout(self.layout2_step1_A3R)
    self.layout_tabsA_step1_3R.addWidget(self.frame2_step1_A3R)

    self.label2_step1_A3R = QLabel("Representative SMILES ID")
    self.label2_step1_A3R.setAlignment(Qt.AlignmentFlag.AlignHCenter)
    self.layout2_step1_A3R.addWidget(self.label2_step1_A3R)

    self.representative_step1_A3R = QLineEdit()
    self.representative_step1_A3R.setPlaceholderText("Existing SMILES ID in reagents A file")
    self.layout2_step1_A3R.addWidget(self.representative_step1_A3R)

    # def lock(self):
    #     self.representative_ID.setEnabled(False)

    ###########

    self.frame3_step1_A3R = QFrame()
    # self.frame3_step1_A3R.setFixedHeight(150)
    self.layout3_step1_A3R = QVBoxLayout()
    self.layout3_step1_A3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.frame3_step1_A3R.setLayout(self.layout3_step1_A3R)
    self.layout_tabsA_step1_3R.addWidget(self.frame3_step1_A3R)

    self.smiles_step1_A3R = QLineEdit()
    self.smiles_step1_A3R.setPlaceholderText("SMILES string for reagent A representative")
    
    self.timer_step1_A3R = QTimer()
    self.timer_step1_A3R.setSingleShot(True)
    self.layout3_step1_A3R.addWidget(self.smiles_step1_A3R)

    self.image_step1_A3R = QLabel()
    self.image_step1_A3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_step1_A3R.addWidget(self.image_step1_A3R)

    self.mw_step1_A3R = QLabel()
    self.mw_step1_A3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_step1_A3R.addWidget(self.mw_step1_A3R)
    self.smiles_step1_A3R.textChanged.connect(lambda:schedule_refresh_image_step1_A3R(self))
    self.timer_step1_A3R.timeout.connect(lambda:refresh_image_step1_A3R(self))

    ###########

    self.frame1_step1_A3R.show()
    self.frame2_step1_A3R.hide()
    self.frame3_step1_A3R.hide()
    self.button1_step1_A3R.setStyleSheet("background-color:lightgrey;")

    ##################################################################################################################################
    # SEPARATOR
    self.separator = QFrame()
    self.separator.setFrameShape(QFrame.Shape.VLine)
    # self.separator.setFrameShadow(QFrame.Sunken)
    self.separator.setLineWidth(2)
    self.separator.setStyleSheet("color: gray;")
    self.layout_tabs_step1_3R.addWidget(self.separator)

    ##################################################################################################################################
    # TABS - REAGENT B
    self.layout_tabsB_step1_3R = QVBoxLayout()
    self.layout_tabs_step1_3R.addLayout(self.layout_tabsB_step1_3R, 50)

    self.tabsB_step1_3R_title = QLabel("Reagent B representative")
    self.tabsB_step1_3R_title.setToolTip("The representative of the full reagents B dataset, named b, to generate Ab chemical libraries")
    self.tabsB_step1_3R_title.setFixedHeight(self.lineheight)
    self.tabsB_step1_3R_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
    self.tabsB_step1_3R_title_font = QFont("Calibri", 10)
    self.tabsB_step1_3R_title_font.setBold(True)
    self.tabsB_step1_3R_title.setFont(self.tabsB_step1_3R_title_font)
    self.layout_tabsB_step1_3R.addWidget(self.tabsB_step1_3R_title)

    ###########

    self.button1_step1_B3R = QPushButton("Smallest")
    self.button2_step1_B3R = QPushButton("Chosen ID")
    self.button3_step1_B3R = QPushButton("Custom SMILES")
    self.button1_step1_B3R.clicked.connect(lambda: (self.frame1_step1_B3R.show(), self.frame2_step1_B3R.hide(), self.frame3_step1_B3R.hide(), 
                                              self.button1_step1_B3R.setStyleSheet("background-color:lightgrey;"),
                                              self.button2_step1_B3R.setStyleSheet("background-color:none;"),
                                              self.button3_step1_B3R.setStyleSheet("background-color:none;"),
                                              ))


    self.button2_step1_B3R.clicked.connect(lambda: (self.frame1_step1_B3R.hide(), self.frame2_step1_B3R.show(), self.frame3_step1_B3R.hide(), 
                                              self.button1_step1_B3R.setStyleSheet("background-color:none;"),
                                              self.button2_step1_B3R.setStyleSheet("background-color:lightgrey;"),
                                              self.button3_step1_B3R.setStyleSheet("background-color:none;"),
                                              ))


    self.button3_step1_B3R.clicked.connect(lambda: (self.frame1_step1_B3R.hide(), self.frame2_step1_B3R.hide(), self.frame3_step1_B3R.show(), 
                                              self.button1_step1_B3R.setStyleSheet("background-color:none;"),
                                              self.button2_step1_B3R.setStyleSheet("background-color:none;"),
                                              self.button3_step1_B3R.setStyleSheet("background-color:lightgrey;"),
                                              ))


    self.buttons_step1_B3R_layout = QHBoxLayout()
    self.buttons_step1_B3R_layout.addWidget(self.button1_step1_B3R)
    self.buttons_step1_B3R_layout.addWidget(self.button2_step1_B3R)
    self.buttons_step1_B3R_layout.addWidget(self.button3_step1_B3R)

    self.layout_tabsB_step1_3R.addLayout(self.buttons_step1_B3R_layout)
    
    ###########

    self.frame1_step1_B3R = QFrame()
    self.layout1_step1_B3R = QVBoxLayout()
    self.frame1_step1_B3R.setLayout(self.layout1_step1_B3R)
    self.layout_tabsB_step1_3R.addWidget(self.frame1_step1_B3R)
    self.label1_step1_B3R = QLabel("No parameters required")
    self.label1_step1_B3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter)
    self.layout1_step1_B3R.addWidget(self.label1_step1_B3R)
    
    ###########

    self.frame2_step1_B3R = QFrame()
    self.layout2_step1_B3R = QVBoxLayout()
    self.layout2_step1_B3R.setAlignment(Qt.AlignmentFlag.AlignVCenter)
    self.frame2_step1_B3R.setLayout(self.layout2_step1_B3R)
    self.layout_tabsB_step1_3R.addWidget(self.frame2_step1_B3R)

    self.label2_step1_B3R = QLabel("Representative SMILES ID")
    self.label2_step1_B3R.setAlignment(Qt.AlignmentFlag.AlignHCenter)
    self.layout2_step1_B3R.addWidget(self.label2_step1_B3R)

    self.representative_step1_B3R = QLineEdit()
    self.representative_step1_B3R.setPlaceholderText("Existing SMILES ID in reagents A file")
    self.layout2_step1_B3R.addWidget(self.representative_step1_B3R)

    # def lock(self):
    #     self.representative_ID.setEnabled(False)

    ###########

    self.frame3_step1_B3R = QFrame()
    self.layout3_step1_B3R = QVBoxLayout()
    self.layout3_step1_B3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.frame3_step1_B3R.setLayout(self.layout3_step1_B3R)
    self.layout_tabsB_step1_3R.addWidget(self.frame3_step1_B3R)

    self.smiles_step1_B3R = QLineEdit()
    self.smiles_step1_B3R.setPlaceholderText("SMILES string for reagent A representative")
    
    self.timer_step1_B3R = QTimer()
    self.timer_step1_B3R.setSingleShot(True)
    self.layout3_step1_B3R.addWidget(self.smiles_step1_B3R)

    self.image_step1_B3R = QLabel()
    self.image_step1_B3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_step1_B3R.addWidget(self.image_step1_B3R)

    self.mw_step1_B3R = QLabel()
    self.mw_step1_B3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_step1_B3R.addWidget(self.mw_step1_B3R)
    self.smiles_step1_B3R.textChanged.connect(lambda:schedule_refresh_image_step1_B3R(self))
    self.timer_step1_B3R.timeout.connect(lambda:refresh_image_step1_B3R(self))

    ###########

    self.frame1_step1_B3R.show()
    self.frame2_step1_B3R.hide()
    self.frame3_step1_B3R.hide()
    self.button1_step1_B3R.setStyleSheet("background-color:lightgrey;")


    ##################################################################################################################################
    # SEPARATOR
    self.separator = QFrame()
    self.separator.setFrameShape(QFrame.Shape.VLine)
    # self.separator.setFrameShadow(QFrame.Sunken)
    self.separator.setLineWidth(2)
    self.separator.setStyleSheet("color: gray;")
    self.layout_tabs_step1_3R.addWidget(self.separator)

    ##################################################################################################################################
    # TABS - REAGENT C
    self.layout_tabsC_step1_3R = QVBoxLayout()
    self.layout_tabs_step1_3R.addLayout(self.layout_tabsC_step1_3R, 50)

    self.tabsC_step1_3R_title = QLabel("Reagent C representative")
    self.tabsC_step1_3R_title.setToolTip("The representative of the full reagents C dataset, named c, to generate Ab chemical libraries")
    self.tabsC_step1_3R_title.setFixedHeight(self.lineheight)
    self.tabsC_step1_3R_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
    self.tabsC_step1_3R_title_font = QFont("Calibri", 10)
    self.tabsC_step1_3R_title_font.setBold(True)
    self.tabsC_step1_3R_title.setFont(self.tabsC_step1_3R_title_font)
    self.layout_tabsC_step1_3R.addWidget(self.tabsC_step1_3R_title)

    ###########

    self.button1_step1_C3R = QPushButton("Smallest")
    self.button2_step1_C3R = QPushButton("Chosen ID")
    self.button3_step1_C3R = QPushButton("Custom SMILES")
    self.button1_step1_C3R.clicked.connect(lambda: (self.frame1_step1_C3R.show(), self.frame2_step1_C3R.hide(), self.frame3_step1_C3R.hide(), 
                                              self.button1_step1_C3R.setStyleSheet("background-color:lightgrey;"),
                                              self.button2_step1_C3R.setStyleSheet("background-color:none;"),
                                              self.button3_step1_C3R.setStyleSheet("background-color:none;"),
                                              ))


    self.button2_step1_C3R.clicked.connect(lambda: (self.frame1_step1_C3R.hide(), self.frame2_step1_C3R.show(), self.frame3_step1_C3R.hide(), 
                                              self.button1_step1_C3R.setStyleSheet("background-color:none;"),
                                              self.button2_step1_C3R.setStyleSheet("background-color:lightgrey;"),
                                              self.button3_step1_C3R.setStyleSheet("background-color:none;"),
                                              ))


    self.button3_step1_C3R.clicked.connect(lambda: (self.frame1_step1_C3R.hide(), self.frame2_step1_C3R.hide(), self.frame3_step1_C3R.show(), 
                                              self.button1_step1_C3R.setStyleSheet("background-color:none;"),
                                              self.button2_step1_C3R.setStyleSheet("background-color:none;"),
                                              self.button3_step1_C3R.setStyleSheet("background-color:lightgrey;"),
                                              ))


    self.buttons_step1_C3R_layout = QHBoxLayout()
    self.buttons_step1_C3R_layout.addWidget(self.button1_step1_C3R)
    self.buttons_step1_C3R_layout.addWidget(self.button2_step1_C3R)
    self.buttons_step1_C3R_layout.addWidget(self.button3_step1_C3R)

    self.layout_tabsC_step1_3R.addLayout(self.buttons_step1_C3R_layout)
    
    ###########

    self.frame1_step1_C3R = QFrame()
    self.layout1_step1_C3R = QVBoxLayout()
    self.frame1_step1_C3R.setLayout(self.layout1_step1_C3R)
    self.layout_tabsC_step1_3R.addWidget(self.frame1_step1_C3R)
    self.label1_step1_C3R = QLabel("No parameters required")
    self.label1_step1_C3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter)
    self.layout1_step1_C3R.addWidget(self.label1_step1_C3R)
    
    ###########

    self.frame2_step1_C3R = QFrame()
    self.layout2_step1_C3R = QVBoxLayout()
    self.layout2_step1_C3R.setAlignment(Qt.AlignmentFlag.AlignVCenter)
    self.frame2_step1_C3R.setLayout(self.layout2_step1_C3R)
    self.layout_tabsC_step1_3R.addWidget(self.frame2_step1_C3R)

    self.label2_step1_C3R = QLabel("Representative SMILES ID")
    self.label2_step1_C3R.setAlignment(Qt.AlignmentFlag.AlignHCenter)
    self.layout2_step1_C3R.addWidget(self.label2_step1_C3R)

    self.representative_step1_C3R = QLineEdit()
    self.representative_step1_C3R.setPlaceholderText("Existing SMILES ID in reagents A file")
    self.layout2_step1_C3R.addWidget(self.representative_step1_C3R)

    # def lock(self):
    #     self.representative_ID.setEnabled(False)

    ###########

    self.frame3_step1_C3R = QFrame()
    self.layout3_step1_C3R = QVBoxLayout()
    self.layout3_step1_C3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.frame3_step1_C3R.setLayout(self.layout3_step1_C3R)
    self.layout_tabsC_step1_3R.addWidget(self.frame3_step1_C3R)

    self.smiles_step1_C3R = QLineEdit()
    self.smiles_step1_C3R.setPlaceholderText("SMILES string for reagent A representative")
    
    self.timer_step1_C3R = QTimer()
    self.timer_step1_C3R.setSingleShot(True)
    self.layout3_step1_C3R.addWidget(self.smiles_step1_C3R)

    self.image_step1_C3R = QLabel()
    self.image_step1_C3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_step1_C3R.addWidget(self.image_step1_C3R)

    self.mw_step1_C3R = QLabel()
    self.mw_step1_C3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_step1_C3R.addWidget(self.mw_step1_C3R)
    self.smiles_step1_C3R.textChanged.connect(lambda:schedule_refresh_image_step1_C3R(self))
    self.timer_step1_C3R.timeout.connect(lambda:refresh_image_step1_C3R(self))

    ###########

    self.frame1_step1_C3R.show()
    self.frame2_step1_C3R.hide()
    self.frame3_step1_C3R.hide()
    self.button1_step1_C3R.setStyleSheet("background-color:lightgrey;")

    ##################################################################################################################################
    # STEP #2 - MEL 2
    ##################################################################################################################################
    # BEST CPDS FROM STEP #1
    self.layout_selected_step2_3R = QVBoxLayout()
    self.layout_tabs_step2_3R_vertical.addLayout(self.layout_selected_step2_3R)

    self.label_step2 = QLabel("Selected compounds from the step #1 to generate the second Minimal Enumeration Library")
    self.label_step2.setAlignment(Qt.AlignmentFlag.AlignHCenter)
    self.layout_selected_step2_3R.addWidget(self.label_step2)

    self.layout_selected_browse_step2_3R = QHBoxLayout()
    self.layout_selected_step2_3R.addLayout(self.layout_selected_browse_step2_3R)

    self.selected_step2_3R = QLineEdit(self)
    self.selected_step2_3R.setPlaceholderText("File path for selected compounds from step #1")
    self.layout_selected_browse_step2_3R.addWidget(self.selected_step2_3R)

    self.browse_step2_3R = QPushButton("Browse", self)
    self.browse_step2_3R.clicked.connect(lambda:open_file_reagent_step2_3R(self))
    self.layout_selected_browse_step2_3R.addWidget(self.browse_step2_3R)
    

    self.layout_tabs_step2_3R = QHBoxLayout()
    self.layout_tabs_step2_3R_vertical.addLayout(self.layout_tabs_step2_3R)

    ######################################################
    # TABS - REAGENT A
    self.layout_tabsA_step2_3R = QVBoxLayout()
    self.layout_tabs_step2_3R.addLayout(self.layout_tabsA_step2_3R, 50)

    self.tabsA_step2_3R_title = QLabel("Reagent A representative")
    self.tabsA_step2_3R_title.setToolTip("The representative of the full reagents A dataset, named a, to generate Ba chemical libraries")
    # self.tabsA_step2_3R_title.setToolTipDuration(100)
    self.tabsA_step2_3R_title.setFixedHeight(self.lineheight)
    self.tabsA_step2_3R_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
    self.tabsA_step2_3R_title_font = QFont("Calibri", 10)
    self.tabsA_step2_3R_title_font.setBold(True)
    self.tabsA_step2_3R_title.setFont(self.tabsA_step2_3R_title_font)
    self.layout_tabsA_step2_3R.addWidget(self.tabsA_step2_3R_title)

    ###########

    self.button1_step2_A3R = QPushButton("Smallest")
    self.button2_step2_A3R = QPushButton("Chosen ID")
    self.button3_step2_A3R = QPushButton("Custom SMILES")

    self.button1_step2_A3R.clicked.connect(lambda: (self.frame1_step2_A3R.show(), self.frame2_step2_A3R.hide(), self.frame3_step2_A3R.hide(), 
                                              self.button1_step2_A3R.setStyleSheet("background-color:lightgrey;"),
                                              self.button2_step2_A3R.setStyleSheet("background-color:none;"),
                                              self.button3_step2_A3R.setStyleSheet("background-color:none;"),
                                              ))

    self.button2_step2_A3R.clicked.connect(lambda: (self.frame1_step2_A3R.hide(), self.frame2_step2_A3R.show(), self.frame3_step2_A3R.hide(), 
                                              self.button1_step2_A3R.setStyleSheet("background-color:none;"),
                                              self.button2_step2_A3R.setStyleSheet("background-color:lightgrey;"),
                                              self.button3_step2_A3R.setStyleSheet("background-color:none;"),
                                              ))

    self.button3_step2_A3R.clicked.connect(lambda: (self.frame1_step2_A3R.hide(), self.frame2_step2_A3R.hide(), self.frame3_step2_A3R.show(), 
                                              self.button1_step2_A3R.setStyleSheet("background-color:none;"),
                                              self.button2_step2_A3R.setStyleSheet("background-color:none;"),
                                              self.button3_step2_A3R.setStyleSheet("background-color:lightgrey;"),
                                              ))
    
    self.buttons_step2_A3R_layout = QHBoxLayout()
    self.buttons_step2_A3R_layout.addWidget(self.button1_step2_A3R)
    self.buttons_step2_A3R_layout.addWidget(self.button2_step2_A3R)
    self.buttons_step2_A3R_layout.addWidget(self.button3_step2_A3R)
    self.layout_tabsA_step2_3R.addLayout(self.buttons_step2_A3R_layout)

    ###########

    self.frame1_step2_A3R = QFrame()
    # self.frame1_step2_A3R.setFixedHeight(150)
    self.layout1_step2_A3R = QVBoxLayout()

    self.frame1_step2_A3R.setLayout(self.layout1_step2_A3R)
    self.layout_tabsA_step2_3R.addWidget(self.frame1_step2_A3R)
    self.label1_step2_A3R = QLabel("No parameters required")
    self.label1_step2_A3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter)
    self.layout1_step2_A3R.addWidget(self.label1_step2_A3R)
    
    ###########

    self.frame2_step2_A3R = QFrame()
    # self.frame2_step2_A3R.setFixedHeight(150)
    self.layout2_step2_A3R = QVBoxLayout()
    self.layout2_step2_A3R.setAlignment(Qt.AlignmentFlag.AlignVCenter)
    self.frame2_step2_A3R.setLayout(self.layout2_step2_A3R)
    self.layout_tabsA_step2_3R.addWidget(self.frame2_step2_A3R)

    self.label2_step2_A3R = QLabel("Representative SMILES ID")
    self.label2_step2_A3R.setAlignment(Qt.AlignmentFlag.AlignHCenter)
    self.layout2_step2_A3R.addWidget(self.label2_step2_A3R)

    self.representative_step2_A3R = QLineEdit()
    self.representative_step2_A3R.setPlaceholderText("Existing SMILES ID in reagents A file")
    self.layout2_step2_A3R.addWidget(self.representative_step2_A3R)

    # def lock(self):
    #     self.representative_ID.setEnabled(False)

    ###########

    self.frame3_step2_A3R = QFrame()
    # self.frame3_step2_A3R.setFixedHeight(150)
    self.layout3_step2_A3R = QVBoxLayout()
    self.layout3_step2_A3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.frame3_step2_A3R.setLayout(self.layout3_step2_A3R)
    self.layout_tabsA_step2_3R.addWidget(self.frame3_step2_A3R)

    self.smiles_step2_A3R = QLineEdit()
    self.smiles_step2_A3R.setPlaceholderText("SMILES string for reagent A representative")
    
    self.timer_step2_A3R = QTimer()
    self.timer_step2_A3R.setSingleShot(True)
    self.layout3_step2_A3R.addWidget(self.smiles_step2_A3R)

    self.image_step2_A3R = QLabel()
    self.image_step2_A3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_step2_A3R.addWidget(self.image_step2_A3R)

    self.mw_step2_A3R = QLabel()
    self.mw_step2_A3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_step2_A3R.addWidget(self.mw_step2_A3R)
    self.smiles_step2_A3R.textChanged.connect(lambda:schedule_refresh_image_step2_A3R(self))
    self.timer_step2_A3R.timeout.connect(lambda:refresh_image_step2_A3R(self))

    ###########

    self.frame1_step2_A3R.show()
    self.frame2_step2_A3R.hide()
    self.frame3_step2_A3R.hide()
    self.button1_step2_A3R.setStyleSheet("background-color:lightgrey;")

    ##################################################################################################################################
    # SEPARATOR
    self.separator = QFrame()
    self.separator.setFrameShape(QFrame.Shape.VLine)
    # self.separator.setFrameShadow(QFrame.Sunken)
    self.separator.setLineWidth(2)
    self.separator.setStyleSheet("color: gray;")
    self.layout_tabs_step2_3R.addWidget(self.separator)

    ##################################################################################################################################
    # TABS - REAGENT B
    self.layout_tabsB_step2_3R = QVBoxLayout()
    self.layout_tabs_step2_3R.addLayout(self.layout_tabsB_step2_3R, 50)

    self.tabsB_step2_3R_title = QLabel("Reagent B representative")
    self.tabsB_step2_3R_title.setToolTip("The representative of the full reagents B dataset, named b, to generate Ab chemical libraries")
    self.tabsB_step2_3R_title.setFixedHeight(self.lineheight)
    self.tabsB_step2_3R_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
    self.tabsB_step2_3R_title_font = QFont("Calibri", 10)
    self.tabsB_step2_3R_title_font.setBold(True)
    self.tabsB_step2_3R_title.setFont(self.tabsB_step2_3R_title_font)
    self.layout_tabsB_step2_3R.addWidget(self.tabsB_step2_3R_title)

    ###########

    self.button1_step2_B3R = QPushButton("Smallest")
    self.button2_step2_B3R = QPushButton("Chosen ID")
    self.button3_step2_B3R = QPushButton("Custom SMILES")
    self.button1_step2_B3R.clicked.connect(lambda: (self.frame1_step2_B3R.show(), self.frame2_step2_B3R.hide(), self.frame3_step2_B3R.hide(), 
                                              self.button1_step2_B3R.setStyleSheet("background-color:lightgrey;"),
                                              self.button2_step2_B3R.setStyleSheet("background-color:none;"),
                                              self.button3_step2_B3R.setStyleSheet("background-color:none;"),
                                              ))


    self.button2_step2_B3R.clicked.connect(lambda: (self.frame1_step2_B3R.hide(), self.frame2_step2_B3R.show(), self.frame3_step2_B3R.hide(), 
                                              self.button1_step2_B3R.setStyleSheet("background-color:none;"),
                                              self.button2_step2_B3R.setStyleSheet("background-color:lightgrey;"),
                                              self.button3_step2_B3R.setStyleSheet("background-color:none;"),
                                              ))


    self.button3_step2_B3R.clicked.connect(lambda: (self.frame1_step2_B3R.hide(), self.frame2_step2_B3R.hide(), self.frame3_step2_B3R.show(), 
                                              self.button1_step2_B3R.setStyleSheet("background-color:none;"),
                                              self.button2_step2_B3R.setStyleSheet("background-color:none;"),
                                              self.button3_step2_B3R.setStyleSheet("background-color:lightgrey;"),
                                              ))


    self.buttons_step2_B3R_layout = QHBoxLayout()
    self.buttons_step2_B3R_layout.addWidget(self.button1_step2_B3R)
    self.buttons_step2_B3R_layout.addWidget(self.button2_step2_B3R)
    self.buttons_step2_B3R_layout.addWidget(self.button3_step2_B3R)

    self.layout_tabsB_step2_3R.addLayout(self.buttons_step2_B3R_layout)
    
    ###########

    self.frame1_step2_B3R = QFrame()
    self.layout1_step2_B3R = QVBoxLayout()
    self.frame1_step2_B3R.setLayout(self.layout1_step2_B3R)
    self.layout_tabsB_step2_3R.addWidget(self.frame1_step2_B3R)
    self.label1_step2_B3R = QLabel("No parameters required")
    self.label1_step2_B3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter)
    self.layout1_step2_B3R.addWidget(self.label1_step2_B3R)
    
    ###########

    self.frame2_step2_B3R = QFrame()
    self.layout2_step2_B3R = QVBoxLayout()
    self.layout2_step2_B3R.setAlignment(Qt.AlignmentFlag.AlignVCenter)
    self.frame2_step2_B3R.setLayout(self.layout2_step2_B3R)
    self.layout_tabsB_step2_3R.addWidget(self.frame2_step2_B3R)

    self.label2_step2_B3R = QLabel("Representative SMILES ID")
    self.label2_step2_B3R.setAlignment(Qt.AlignmentFlag.AlignHCenter)
    self.layout2_step2_B3R.addWidget(self.label2_step2_B3R)

    self.representative_step2_B3R = QLineEdit()
    self.representative_step2_B3R.setPlaceholderText("Existing SMILES ID in reagents A file")
    self.layout2_step2_B3R.addWidget(self.representative_step2_B3R)

    # def lock(self):
    #     self.representative_ID.setEnabled(False)

    ###########

    self.frame3_step2_B3R = QFrame()
    self.layout3_step2_B3R = QVBoxLayout()
    self.layout3_step2_B3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.frame3_step2_B3R.setLayout(self.layout3_step2_B3R)
    self.layout_tabsB_step2_3R.addWidget(self.frame3_step2_B3R)

    self.smiles_step2_B3R = QLineEdit()
    self.smiles_step2_B3R.setPlaceholderText("SMILES string for reagent A representative")
    
    self.timer_step2_B3R = QTimer()
    self.timer_step2_B3R.setSingleShot(True)
    self.layout3_step2_B3R.addWidget(self.smiles_step2_B3R)

    self.image_step2_B3R = QLabel()
    self.image_step2_B3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_step2_B3R.addWidget(self.image_step2_B3R)

    self.mw_step2_B3R = QLabel()
    self.mw_step2_B3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_step2_B3R.addWidget(self.mw_step2_B3R)
    self.smiles_step2_B3R.textChanged.connect(lambda:schedule_refresh_image_step2_B3R(self))
    self.timer_step2_B3R.timeout.connect(lambda:refresh_image_step2_B3R(self))

    ###########

    self.frame1_step2_B3R.show()
    self.frame2_step2_B3R.hide()
    self.frame3_step2_B3R.hide()
    self.button1_step2_B3R.setStyleSheet("background-color:lightgrey;")


    ##################################################################################################################################
    # SEPARATOR
    self.separator = QFrame()
    self.separator.setFrameShape(QFrame.Shape.VLine)
    # self.separator.setFrameShadow(QFrame.Sunken)
    self.separator.setLineWidth(2)
    self.separator.setStyleSheet("color: gray;")
    self.layout_tabs_step2_3R.addWidget(self.separator)

    ##################################################################################################################################
    # TABS - REAGENT C
    self.layout_tabsC_step2_3R = QVBoxLayout()
    self.layout_tabs_step2_3R.addLayout(self.layout_tabsC_step2_3R, 50)

    self.tabsC_step2_3R_title = QLabel("Reagent C representative")
    self.tabsC_step2_3R_title.setToolTip("The representative of the full reagents C dataset, named c, to generate Ab chemical libraries")
    self.tabsC_step2_3R_title.setFixedHeight(self.lineheight)
    self.tabsC_step2_3R_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
    self.tabsC_step2_3R_title_font = QFont("Calibri", 10)
    self.tabsC_step2_3R_title_font.setBold(True)
    self.tabsC_step2_3R_title.setFont(self.tabsC_step2_3R_title_font)
    self.layout_tabsC_step2_3R.addWidget(self.tabsC_step2_3R_title)

    ###########

    self.button1_step2_C3R = QPushButton("Smallest")
    self.button2_step2_C3R = QPushButton("Chosen ID")
    self.button3_step2_C3R = QPushButton("Custom SMILES")
    self.button1_step2_C3R.clicked.connect(lambda: (self.frame1_step2_C3R.show(), self.frame2_step2_C3R.hide(), self.frame3_step2_C3R.hide(), 
                                              self.button1_step2_C3R.setStyleSheet("background-color:lightgrey;"),
                                              self.button2_step2_C3R.setStyleSheet("background-color:none;"),
                                              self.button3_step2_C3R.setStyleSheet("background-color:none;"),
                                              ))


    self.button2_step2_C3R.clicked.connect(lambda: (self.frame1_step2_C3R.hide(), self.frame2_step2_C3R.show(), self.frame3_step2_C3R.hide(), 
                                              self.button1_step2_C3R.setStyleSheet("background-color:none;"),
                                              self.button2_step2_C3R.setStyleSheet("background-color:lightgrey;"),
                                              self.button3_step2_C3R.setStyleSheet("background-color:none;"),
                                              ))


    self.button3_step2_C3R.clicked.connect(lambda: (self.frame1_step2_C3R.hide(), self.frame2_step2_C3R.hide(), self.frame3_step2_C3R.show(), 
                                              self.button1_step2_C3R.setStyleSheet("background-color:none;"),
                                              self.button2_step2_C3R.setStyleSheet("background-color:none;"),
                                              self.button3_step2_C3R.setStyleSheet("background-color:lightgrey;"),
                                              ))


    self.buttons_step2_C3R_layout = QHBoxLayout()
    self.buttons_step2_C3R_layout.addWidget(self.button1_step2_C3R)
    self.buttons_step2_C3R_layout.addWidget(self.button2_step2_C3R)
    self.buttons_step2_C3R_layout.addWidget(self.button3_step2_C3R)

    self.layout_tabsC_step2_3R.addLayout(self.buttons_step2_C3R_layout)
    
    ###########

    self.frame1_step2_C3R = QFrame()
    self.layout1_step2_C3R = QVBoxLayout()
    self.frame1_step2_C3R.setLayout(self.layout1_step2_C3R)
    self.layout_tabsC_step2_3R.addWidget(self.frame1_step2_C3R)
    self.label1_step2_C3R = QLabel("No parameters required")
    self.label1_step2_C3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter)
    self.layout1_step2_C3R.addWidget(self.label1_step2_C3R)
    
    ###########

    self.frame2_step2_C3R = QFrame()
    self.layout2_step2_C3R = QVBoxLayout()
    self.layout2_step2_C3R.setAlignment(Qt.AlignmentFlag.AlignVCenter)
    self.frame2_step2_C3R.setLayout(self.layout2_step2_C3R)
    self.layout_tabsC_step2_3R.addWidget(self.frame2_step2_C3R)

    self.label2_step2_C3R = QLabel("Representative SMILES ID")
    self.label2_step2_C3R.setAlignment(Qt.AlignmentFlag.AlignHCenter)
    self.layout2_step2_C3R.addWidget(self.label2_step2_C3R)

    self.representative_step2_C3R = QLineEdit()
    self.representative_step2_C3R.setPlaceholderText("Existing SMILES ID in reagents A file")
    self.layout2_step2_C3R.addWidget(self.representative_step2_C3R)

    # def lock(self):
    #     self.representative_ID.setEnabled(False)

    ###########

    self.frame3_step2_C3R = QFrame()
    self.layout3_step2_C3R = QVBoxLayout()
    self.layout3_step2_C3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.frame3_step2_C3R.setLayout(self.layout3_step2_C3R)
    self.layout_tabsC_step2_3R.addWidget(self.frame3_step2_C3R)

    self.smiles_step2_C3R = QLineEdit()
    self.smiles_step2_C3R.setPlaceholderText("SMILES string for reagent A representative")
    
    self.timer_step2_C3R = QTimer()
    self.timer_step2_C3R.setSingleShot(True)
    self.layout3_step2_C3R.addWidget(self.smiles_step2_C3R)

    self.image_step2_C3R = QLabel()
    self.image_step2_C3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_step2_C3R.addWidget(self.image_step2_C3R)

    self.mw_step2_C3R = QLabel()
    self.mw_step2_C3R.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
    self.layout3_step2_C3R.addWidget(self.mw_step2_C3R)
    self.smiles_step2_C3R.textChanged.connect(lambda:schedule_refresh_image_step2_C3R(self))
    self.timer_step2_C3R.timeout.connect(lambda:refresh_image_step2_C3R(self))

    ###########

    self.frame1_step2_C3R.show()
    self.frame2_step2_C3R.hide()
    self.frame3_step2_C3R.hide()
    self.button1_step2_C3R.setStyleSheet("background-color:lightgrey;")

    ##################################################################################################################################
    # STEP #3 - FEL
    ##################################################################################################################################

    self.layout_step3_3R = QVBoxLayout()
    self.layout_tabs_step3_3R.addLayout(self.layout_step3_3R)

    self.label_step3 = QLabel("Selected compounds from the step #2 to generate the final Focused Enumeration Library")
    self.label_step3.setAlignment(Qt.AlignmentFlag.AlignHCenter)
    self.layout_step3_3R.addWidget(self.label_step3)

    self.layout_step3_browse_3R = QHBoxLayout()
    self.layout_step3_3R.addLayout(self.layout_step3_browse_3R)

    self.selected_step3_3R = QLineEdit(self)
    self.selected_step3_3R.setPlaceholderText("File path for selected compounds from step #2")
    self.layout_step3_browse_3R.addWidget(self.selected_step3_3R)

    self.browse_step3_3R = QPushButton("Browse", self)
    self.browse_step3_3R.clicked.connect(lambda:open_file_reagent_step3_3R(self))
    self.layout_step3_browse_3R.addWidget(self.browse_step3_3R)



    ##################################################################################################################################