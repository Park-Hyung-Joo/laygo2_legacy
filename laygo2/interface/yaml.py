#!/usr/bin/python
########################################################################################################################
#
# Copyright (c) 2020, Nifty Chips Laboratory, Hanyang University
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
# following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following
#   disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
#    following disclaimer in the documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
########################################################################################################################

"""
This module implements interfaces with yaml files.

"""

__author__ = "Jaeduk Han"
__maintainer__ = "Jaeduk Han"
__status__ = "Prototype"

import yaml
import os.path
import laygo2
import numpy as np

def export_template(template, filename, mode='append'):
    """Export a template to a yaml file.

    Parameters
    ----------
    template: laygo2.object.template.Template
        The template object to be exported.        
    filename: str
        The name of the yaml file.
    mode: str
        If 'append', it adds the template entry without erasing the 
        preexisting file.

    Example
    -------
    >>> import laygo2
    >>> from laygo2.object.physical import Pin
    >>> from laygo2.object.template import NativeInstanceTemplate
    >>> p = dict()
    >>> p['i'] = Pin(xy=[[0, 0], [10, 10]], layer=['M1', 'drawing'],
    >>>                  netname='i')
    >>> p['o'] = Pin(xy=[[90, 90], [100, 100]], layer=['M1', 'drawing'],
    >>>                  netname='o')
    >>> nt = NativeInstanceTemplate(libname='mylib', cellname='mytemp',
    >>>                                 bbox=[[0, 0], [100, 100]], pins=p)
    >>> laygo2.interface.yaml.export_template(nt, filename="mytemplates.yaml")    
    Your design was translated into YAML format.
    {'mylib': {
        'mytemp': {
            'libname': 'mylib', 
            'cellname': 'mytemp', 
            'bbox': [[0, 0], [100, 100]], 
            'pins': {
                'i': {
                    'xy': [[0, 0], [10, 10]], 
                    'layer': ['M1', 'drawing'], 
                    'name': None, 
                    'netname': 'i'
                    }, 
                'o': {
                    'xy': [[90, 90], [100, 100]], 
                    'layer': ['M1', 'drawing'], 
                    'name': None, 
                    'netname': 'o'
    }}}}}
    """
    libname = template.libname
    cellname = template.cellname
    pins = template.pins()
    # sub_blocks = template.sub_blocks # sub_blocks[instName] == inst(not template)
    # obstacles = template.obstacles # obstacles = [Rect0, Rect1, Rect2 ...]
    db = dict()
    if mode == 'append':  # in append mode, the template is appended to 'filename' if the file exists.
        if os.path.exists(filename):
            with open(filename, 'r') as stream:
                db = yaml.load(stream, Loader=yaml.FullLoader)
        else:
            f_new = open(filename, "w")
            f_new.write(f"{libname}:\n")
            f_new.write(f"    dummy:\n")
            f_new.write(f"        bbox:\n")
            f_new.write(f"        - - 0\n")
            f_new.write(f"          - 0\n")
            f_new.write(f"        - - 0\n")
            f_new.write(f"          - 0\n")
            f_new.write(f"        cellname: dummy\n")
            f_new.write(f"        libname: {libname}\n")

            f_new.close()
            with open(filename, 'r') as stream:
                db = yaml.load(stream, Loader=yaml.FullLoader)
    if db == None:
        db[libname] = dict()
    elif libname not in db:
        db[libname] = dict()
    db[libname][cellname] = template.export_to_dict()
    with open(filename, 'w') as stream:
        yaml.dump(db, stream)
    #print("Your design was translated into YAML format.")
    return db

#filename=libname+'_templates.yaml'
def export_design(design, filename, obs_layers:list=None, mode='append'):

    libname = design.libname
    cellname = design.cellname
    pins:dict = design.pins

    db = dict()
    if mode == 'append':  # in append mode, the template is appended to 'filename' if the file exists.
        if os.path.exists(filename):
            with open(filename, 'r') as stream:
                db = yaml.load(stream, Loader=yaml.FullLoader)
        else:
            f_new = open(filename, "w")
            f_new.write(f"{libname}:\n")
            f_new.write(f"    dummy:\n")
            f_new.write(f"        bbox:\n")
            f_new.write(f"        - - 0\n")
            f_new.write(f"          - 0\n")
            f_new.write(f"        - - 0\n")
            f_new.write(f"          - 0\n")
            f_new.write(f"        cellname: dummy\n")
            f_new.write(f"        libname: {libname}\n")

            f_new.close()
            with open(filename, 'r') as stream:
                db = yaml.load(stream, Loader=yaml.FullLoader)
    if libname not in db:
        db[libname] = dict()
    # export basic information
    nat_temp = design.export_to_template(libname=libname, cellname=cellname)
    db[libname][cellname] = nat_temp.export_to_dict()
    # export subblocks
    db[libname][cellname]['sub_blocks'] = dict()
    for _instName, _inst in design.instances.items():
        if "NoName" in _instName: # via or other instances which is not sub-block
            continue
        inst = dict()
        inst['name'] = _instName
        inst['cellname'] = _inst.cellname
        inst['libname'] = _inst.libname
        inst['xy'] = _inst.xy.tolist()
        inst['transform'] = _inst.transform
        db[libname][cellname]['sub_blocks'][_instName] = inst
    # export obstacles
    db[libname][cellname]['obstacles'] = list()
    for _layer in obs_layers:
        for _metal in design.get_matched_rects_by_layer([_layer,'drawing']):
            metal = dict()
        #    print(_metal.xy)
            metal['xy'] = _metal.xy.tolist()
            metal['layer'] = _layer
            if _metal.netname is not None:
                metal['netname'] = _metal.netname
            db[libname][cellname]['obstacles'].append(metal)

    with open(filename, 'w') as stream:
        yaml.dump(db, stream)
    #print("Your design was translated into YAML format.")
    return db

def import_template(filename):
    """Import templates from a yaml file.

    Parameters
    ----------
    filename: str
        The name of the yaml file.
    
    Example
    -------
    >>> import laygo2
    >>> from laygo2.object.physical import Pin
    >>> from laygo2.object.template import NativeInstanceTemplate
    >>> p = dict()
    >>> p['i'] = Pin(xy=[[0, 0], [10, 10]], layer=['M1', 'drawing'],
    >>>                  netname='i')
    >>> p['o'] = Pin(xy=[[90, 90], [100, 100]], layer=['M1', 'drawing'],
    >>>                  netname='o')
    >>> nt = NativeInstanceTemplate(libname='mylib', cellname='mytemp',
    >>>                                 bbox=[[0, 0], [100, 100]], pins=p)
    >>> laygo2.interface.yaml.export_template(nt, filename="mytemplates.yaml")    
    >>> # Import the template back to python.
    >>> my_tlib = laygo2.interface.yaml.import_template("mytemplates.yaml")
    >>> print(my_tlib)
    <laygo2.object.database.TemplateLibrary object at 0x000001FE3440A410> 
        name: mylib, 
        params: None       
        elements: {
            'mytemp': <laygo2.object.template.NativeInstanceTemplate object at 0x000001FE3440A2C0>
        }
    """
    # load yaml file
    if os.path.exists(filename):
        with open(filename, 'r') as stream:
            db = yaml.load(stream, Loader=yaml.FullLoader)
    else:
        print("no such file: "+filename)
    libname = list(db.keys())[0]  # assuming there's only one library defined in each file.
    # create template library
    tlib = laygo2.object.database.TemplateLibrary(name=libname)
    # read out the yaml file entries and build template objects
    for tn, tdict in db[libname].items():
        pins = dict()
        if 'pins' in tdict:
            for pinname, pdict in tdict['pins'].items():
                pins[pinname] = laygo2.object.Pin(xy=pdict['xy'], layer=pdict['layer'], netname=pdict['netname'])
        t = laygo2.object.NativeInstanceTemplate(libname=libname, cellname=tn, bbox=tdict['bbox'], pins=pins)
        tlib.append(t)
    return tlib

