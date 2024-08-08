##############################################################################################################################################
# SATELLiTES - KNIME extension
# Copyright (C) Corentin BEDART, 2024
##############################################################################################################################################

import knime.extension as knext
from satellites_extension import *

import os
import shutil
import pandas as pd

import SATELLiTES.SATELLiTES_enumeration as stenum
import SATELLiTES.SATELLiTES_gui as stgui

##############################################################################################################################################

main_category = knext.category(path="/", level_id="satellites", name="SATELLiTES", description="SATELLiTES Extension", icon="icon.png")

#############################~

category_utilities = knext.category(path=main_category, level_id="satellites_utilities", name="Utilities", description="SATELLiTES - Utilities", icon="icon.png")
category_2reagents = knext.category(path=main_category, level_id="satellites_2reagents", name="2-reagents", description="SATELLiTES - 2-reagents", icon="icon.png")
category_3reagents = knext.category(path=main_category, level_id="satellites_3reagents", name="3-reagents", description="SATELLiTES - 3-reagents", icon="icon.png")
category_multiprocessing = knext.category(path=main_category, level_id="satellites_multiprocessing", name="Multiprocessing nodes", description="SATELLiTES - Multiprocessing", icon="icon.png")

##############################################################################################################################################

class SatellitesObjectSpec(knext.PortObjectSpec):
    def __init__(self, df: pd.DataFrame) -> None:
        self._df = df

    def serialize(self) -> dict:
        return {"spec_data": self._df.to_json()}

    @classmethod
    def deserialize(cls, data: dict) -> "SatellitesObjectSpec":
        cls(pd.read_json(data["spec_data"]))

    @property
    # def df(self) -> knext.Table:
    #     return knext.Table.from_pandas(self._df)
    def df(self) -> pd.DataFrame:
        return self._df

#############################~

class SatellitesObject(knext.PortObject):
    def __init__(self, spec: SatellitesObjectSpec, model) -> None:
        super().__init__(spec)
        self._model = model

    def serialize(self) -> bytes:
        return pickle.dumps(self._model)
    
    @property
    def spec(self) -> SatellitesObjectSpec:
        return super().spec
    
    @classmethod
    def deserialize(cls, spec: SatellitesObjectSpec, data: bytes) -> "SatellitesObject":
        model = pickle.loads(data)
        return cls(spec, model)
    
    def predict(self, data):
        return self._model.predict(data)

#############################~

port_type_SATELLiTES = knext.port_type(name="SATELLiTES port type", object_class=SatellitesObject, spec_class=SatellitesObjectSpec)

##############################################################################################################################################

# Utilities
import nodes.satellites_chosen
import nodes.satellites_representative_smallest
import nodes.satellites_representative_byid
import nodes.satellites_representative_custom

# 2-reagents
import nodes.satellites_2reagents_step1
import nodes.satellites_2reagents_step2

# 3-reagents
import nodes.satellites_3reagents_step1
import nodes.satellites_3reagents_step2
import nodes.satellites_3reagents_step3

# Multiprocessing
import nodes.satellites_2reagents_step1_MP
import nodes.satellites_2reagents_step2_MP
import nodes.satellites_3reagents_step1_MP
import nodes.satellites_3reagents_step2_MP
import nodes.satellites_3reagents_step3_MP

##############################################################################################################################################


