from collections import deque
import yaml
import os
from pickle import FALSE
import laygo2
import numpy as np
import laygo2.util.transform as tf
from laygo2.object.physical import Rect, Pin

class rcNode:
    def __init__(self, rc_num):
        self.rc_num = rc_num
        self.dot_list = list()
        self.metal_list = list()

class metalNode:
    def __init__(self, layer_name, mn, metal_ref=None):
        self.layer=layer_name
        self.mn = mn
        self.net_name = set([None])
        self.metal_seg = list()
        self.via_list = list()
        self.visited = False
        self.metal_ref = metal_ref
class viaNode:
    def __init__(self, layer_pair, mn):
        self.mn = mn
        self.layer_pair=layer_pair

class netMap_basic:
    def __init__(self, rc, layer_name, grid):
        self.rc = rc
        self.grid = grid
        self.layer = layer_name
        self.deviants = list()
    def search_rc(self, rc_num):
        start=0
        end=len(self.rc)-1
        while start <= end:
            mid = (start+end)//2
            if self.rc[mid].rc_num == rc_num:
                return mid
            elif self.rc[mid].rc_num > rc_num:
                end = mid-1
            else:
                start = mid+1
        return None
    
    def get_rc_index(self, new_rc_num):
        start = 0
        end = len(self.rc)-1
        if len(self.rc) == 0:
            return 0
        
        while start+1 < end:
            mid = (start+end)//2
            if self.rc[mid].rc_num > new_rc_num:
                end = mid
            elif self.rc[mid].rc_num < new_rc_num:
                start = mid
            else:
                return mid
            
        if self.rc[start].rc_num > new_rc_num:
            return start
        elif self.rc[end].rc_num < new_rc_num:
            return end+1
        else:
            return end


class netMap_vert(netMap_basic):
    def __init__(self, layer_name, grid):
        self.cols = list()
        self.cols_unmerged = list()
        super().__init__(self.cols, layer_name, grid)
        # self.rc = self.cols
        # self.layer = layer_name
        # self.grid = grid
    def search_metal_index(self, col_index, obj_y1): #col -> self.col 
        col = self.cols[col_index]
        start = 0
        end = len(col.metal_list)-1

        if len(col.metal_list) == 0:
            return 0
        
        while start+1 < end:
            mid = (start+end)//2
            if col.metal_list[mid].mn[0][1] > obj_y1:       
                end = mid
            elif col.metal_list[mid].mn[0][1] < obj_y1:      
                start = mid
            else:
                return mid+1

        if col.metal_list[start].mn[0][1] > obj_y1:       
            return start
        elif col.metal_list[end].mn[0][1] <= obj_y1:      
            return end+1
        else:
            return end

    def insert_metal(self, mn, net_name=None, metal_ref=None):
        if mn[0][0] < mn[1][0]:
            x1=mn[0][0]
            x2=mn[1][0]
        else:
            x1=mn[1][0]
            x2=mn[0][0]
        
        if mn[0][1] < mn[1][1]:
            y1=mn[0][1]
            y2=mn[1][1]
        else:
            y1=mn[1][1]
            y2=mn[0][1]


        for i in range (x1, x2+1):
            col_index = self.search_rc(i)
            if col_index is None:
                #self.cols에 rcNode(i)를 순서에 맞춰 삽입
                new_col_index = self.get_rc_index(i)
                self.cols.insert(new_col_index, rcNode(i))
                self.cols_unmerged.insert(new_col_index, rcNode(i))

                #red, blue list 추가
                self.cols[new_col_index].dot_list.append((y1, 'b'))
                self.cols[new_col_index].dot_list.append((y2, 'r'))

                new_node = metalNode(self.layer, [[i, y1], [i, y2]], metal_ref)
                new_node.net_name.add(net_name)
                #새로운 metalnode를 unmerged list에 추가
                self.cols_unmerged[new_col_index].metal_list.append(new_node)

            else:
                #red, blue list 추가
                self.cols[col_index].dot_list.append((y1, 'b'))
                self.cols[col_index].dot_list.append((y2, 'r'))

                new_node = metalNode(self.layer, [[i, y1], [i, y2]], metal_ref)
                new_node.net_name.add(net_name)
                #새로운 metalnode를 unmerged list에 추가
                self.cols_unmerged[col_index].metal_list.append(new_node)

                
    #insert가 끝난 뒤 마지막에 호출해주면 된다.
    def merge(self):
        #blue랑 red 비교/병합
        for col in self.cols:
            #dot들 정렬
            col.dot_list = sorted(col.dot_list, key = lambda x : x[0] + ord(x[1])/128)
            br_diff = 0
            start = 0
            end = 0
            for dot in col.dot_list:
                if br_diff == 0:
                    start = dot[0]
                if dot[1] == 'b':
                    br_diff = br_diff + 1
                else:
                    br_diff = br_diff - 1

                if br_diff < 0:
                    print("error")
                    break
                if br_diff == 0:
                    end = dot[0]
                    col.metal_list.append(metalNode(self.layer, [[col.rc_num, start], [col.rc_num, end]]))
        col_index = 0
        #unmerged list와 병합해 기타 데이터 추가
        for col_unmerged in self.cols_unmerged:
            for unmerged_metal in col_unmerged.metal_list:
                metal_index = self.search_metal_index(col_index, unmerged_metal.mn[0][1])-1
                self.cols[col_index].metal_list[metal_index].metal_seg.append(unmerged_metal)
                for _netname in unmerged_metal.net_name:
                    self.cols[col_index].metal_list[metal_index].net_name.add(_netname)
            col_index = col_index + 1

    def is_occupied(self, mn):
        col_index = self.search_rc(mn[0])
        if col_index is None:
            # check deviants
            for deviant in self.deviants:
                _mn = self.grid.mn(deviant)
                if _mn[0][0] <= mn[0] and mn[0] <= _mn[1][0] and _mn[0][1] <= mn[1] and mn[1] <= _mn[1][1]:
                    # one of deviants is occupying mn
                    return [_mn, [deviant.netname]]
            # Not found
            return None
        else:
            col = self.cols[col_index]
            if len(col.metal_list) == 0:
                return None
            start = 0
            end = len(col.metal_list)-1 # end >= 0
            mid = end
            mid_pre = mid
            while start < end:
                mid = (start+end)//2
                if mid_pre == mid:
                    break
                elif col.metal_list[mid].mn[0][1] > mn[1]:       
                    end = mid
                elif col.metal_list[mid].mn[1][1] < mn[1]:      
                    start = mid
                else:
                    break
                mid_pre = mid
            if col.metal_list[mid].mn[0][1] <= mn[1] and mn[1] <= col.metal_list[mid].mn[1][1]:
                _metal = col.metal_list[mid]
                return [_metal.mn, list(_metal.net_name)]
                # found occupying metal
            else:
                # not found -> check deviant
                for deviant in self.deviants:
                    _mn = self.grid.mn(deviant)
                    if _mn[0][0] <= mn[0] and mn[0] <= _mn[1][0] and _mn[0][1] <= mn[1] and mn[1] <= _mn[1][1]:
                        # one of deviants is occupying mn
                        return [_mn, [deviant.netname]]
                # Not found
                return None

class netMap_hor(netMap_basic):
    def __init__(self, layer_name, grid):
        self.rows = list()
        self.rows_unmerged = list()
        super().__init__(self.rows, layer_name, grid)
        # self.rc = self.rows
        # self.layer=layer_name
        # self.grid = grid
    def search_metal_index(self, row_index, obj_x1): #row -> self.row
        row = self.rows[row_index]
        start = 0
        end = len(row.metal_list)-1

        if len(row.metal_list) == 0:
            return 0
        
        while start+1 < end:
            mid = (start+end)//2
            if row.metal_list[mid].mn[0][0] > obj_x1:       
                end = mid
            elif row.metal_list[mid].mn[0][0] < obj_x1:      
                start = mid
            else:
                return mid+1

        if row.metal_list[start].mn[0][0] > obj_x1:       
            return start
        elif row.metal_list[end].mn[0][0] <= obj_x1:      
            return end+1
        else:
            return end

    def insert_metal(self, mn, net_name=None, metal_ref=None):
        if mn[0][0] < mn[1][0]:
            x1=mn[0][0]
            x2=mn[1][0]
        else:
            x1=mn[1][0]
            x2=mn[0][0]
        
        if mn[0][1] < mn[1][1]:
            y1=mn[0][1]
            y2=mn[1][1]
        else:
            y1=mn[1][1]
            y2=mn[0][1]

        for i in range (y1, y2+1):
            row_index = self.search_rc(i)
            if row_index is None:
                #self.rows에 rowNode(i)를 순서에 맞춰 삽입
                new_row_index = self.get_rc_index(i)
                self.rows.insert(new_row_index, rcNode(i))
                self.rows_unmerged.insert(new_row_index, rcNode(i))

                #red, blue list 추가
                self.rows[new_row_index].dot_list.append((x1, 'b'))
                self.rows[new_row_index].dot_list.append((x2, 'r'))

                new_node = metalNode(self.layer, [[x1, i], [x2, i]], metal_ref)
                new_node.net_name.add(net_name)
                #새로운 metalnode를 unmerged list에 추가
                self.rows_unmerged[new_row_index].metal_list.append(new_node)

            else:
                #red, blue list 추가
                self.rows[row_index].dot_list.append((x1, 'b'))
                self.rows[row_index].dot_list.append((x2, 'r'))

                new_node = metalNode(self.layer, [[x1, i], [x2, i]], metal_ref)
                new_node.net_name.add(net_name)
                #새로운 metalnode를 unmerged list에 추가
                self.rows_unmerged[row_index].metal_list.append(new_node)

    def merge(self):
        #blue랑 red 비교/병합
        for row in self.rows:
            #dot들 정렬
            row.dot_list = sorted(row.dot_list, key = lambda x : x[0] + ord(x[1])/128)
            br_diff = 0
            start = 0
            end = 0
            for dot in row.dot_list:
                if br_diff == 0:
                    start = dot[0]
                if dot[1] == 'b':
                    br_diff = br_diff + 1
                else:
                    br_diff = br_diff - 1

                if br_diff < 0:
                    print("error")
                    break
                if br_diff == 0:
                    end = dot[0]
                    row.metal_list.append(metalNode(self.layer, [[start, row.rc_num], [end, row.rc_num]]))

        row_index = 0
        #unmerged list와 병합해 기타 데이터 추가
        for row_unmerged in self.rows_unmerged:
            for unmerged_metal in row_unmerged.metal_list:
                metal_index = self.search_metal_index(row_index, unmerged_metal.mn[0][0])-1
                self.rows[row_index].metal_list[metal_index].metal_seg.append(unmerged_metal)
                for _netname in unmerged_metal.net_name:
                    self.rows[row_index].metal_list[metal_index].net_name.add(_netname)
            row_index = row_index+1

    def is_occupied(self, mn):
        row_index = self.search_rc(mn[1])
        if row_index is None:
            # check deviants
            for deviant in self.deviants:
                _mn = self.grid.mn(deviant)
                if _mn[0][0] <= mn[0] and mn[0] <= _mn[1][0] and _mn[0][1] <= mn[1] and mn[1] <= _mn[1][1]:
                    # one of deviants is occupying mn
                    return [_mn, [deviant.netname]]
            # Not found
            return None
        else:
            row = self.rows[row_index]
            if len(row.metal_list) == 0:
                return None
            start = 0
            end = len(row.metal_list)-1 # end >= 0
            mid = end
            mid_pre = mid
            while start < end:
                mid = (start+end)//2
                if mid_pre == mid:
                    break
                elif row.metal_list[mid].mn[0][0] > mn[0]:       
                    end = mid
                elif row.metal_list[mid].mn[1][0] < mn[0]:      
                    start = mid
                else:
                    break
                mid_pre = mid
                
            if row.metal_list[mid].mn[0][0] <= mn[0] and mn[0] <= row.metal_list[mid].mn[1][0]:
                _metal = row.metal_list[mid]
                return [_metal.mn, list(_metal.net_name)]
                # found occupying metal
            else:
                # not found -> check deviant
                for deviant in self.deviants:
                    _mn = self.grid.mn(deviant)
                    if _mn[0][0] <= mn[0] and mn[0] <= _mn[1][0] and _mn[0][1] <= mn[1] and mn[1] <= _mn[1][1]:
                        # one of deviants is occupying mn
                        return [_mn, [deviant.netname]]
                # Not found
                return None

class NetMap:
#member : type
#layers : dict
#layers_orient : dict
#pins : list
#_via_table : dict
#grid : laygo2.object.grid
#orient : bool
    def __init__(self, via_table:dict, grid_table:dict, orient_first="vertical", layer_names:list=['M1','M2','M3','M4','M5'], net_ignore = [], libs:dict = None):
        self.layers=dict()
        self.layers_orient=dict()
        self.pins=list()
        self.terminals=list()
        self.vias=list()
        self._via_table=via_table
        self.grids=grid_table
        self.deviants = list() # those metals do not follow the orient of their layer TODO: move this list to vertical & horizontal map
        self.net_ignore = set(net_ignore)
        self.libs = libs
        if orient_first=="vertical":
            self.orient_order = True
            for idx in range(len(layer_names)//2):
                self.layers[layer_names[idx*2]] = netMap_vert(layer_names[idx*2], grid_table[layer_names[idx*2]])
                self.layers[layer_names[idx*2+1]] = netMap_hor(layer_names[idx*2+1], grid_table[layer_names[idx*2+1]])
                self.layers_orient[layer_names[idx*2]] = "vertical"
                self.layers_orient[layer_names[idx*2+1]] = "horizontal"
            if len(layer_names)%2 == 1:
                self.layers[layer_names[-1]] = netMap_vert(layer_names[-1], grid_table[layer_names[-1]])
                self.layers_orient[layer_names[-1]] = "vertical"

        elif orient_first=="horizontal":
            self.orient_order = False
            for idx in range(len(layer_names)//2):
                self.layers[layer_names[idx*2]] = netMap_hor(layer_names[idx*2], grid_table[layer_names[idx*2]])
                self.layers[layer_names[idx*2+1]] = netMap_vert(layer_names[idx*2+1], grid_table[layer_names[idx*2+1]])
                self.layers_orient[layer_names[idx*2]] = "horizontal"
                self.layers_orient[layer_names[idx*2+1]] = "vertical"
            if len(layer_names)%2 == 1:
                self.layers[layer_names[-1]] = netMap_hor(layer_names[-1], grid_table[layer_names[-1]])
                self.layers_orient[layer_names[-1]] = "horizontal"
        else:
            print("wrong orient")
    
    def translate_rect_vinst(self, elem, master): # translate mn of metals of (virtual)instances
    # leaf node of export tree
        mxy = master.xy
        mtf = master.transform
        _xy = np.sort(elem.xy, axis=0)  # make sure obj.xy is sorted
        _xy = mxy + np.dot(_xy, tf.Mt(mtf).T)
        return self.grids[elem.layer[0]].mn(_xy)
        #self.layers[elem.layer[0]].insert_metal(self.grid.mn(_xy), netname=elem.netname) # TODO: insert metal로 바꿀 것 
    
    def insert_metal(self, metal):
        _mn = self.grids[metal.layer[0]].mn(metal)
        if _mn[0][0] == None or _mn[0][1] == None:
            print("mn error: layer: ", metal.layer[0],"xy: ",metal.xy.tolist(),"mn: ",_mn.tolist())
            self.layers[metal.layer[0]].deviants.append(metal)
            return None
        if self.layers_orient[metal.layer[0]] == "vertical" and abs(_mn[1][0]-_mn[0][0]) <= abs(_mn[1][1]-_mn[0][1]): # mn_vetical >= mn_horizontal
            return self.layers[metal.layer[0]].insert_metal(self.grids[metal.layer[0]].mn(metal), net_name=metal.netname, metal_ref=metal)
        elif self.layers_orient[metal.layer[0]] == "horizontal" and abs(_mn[1][0]-_mn[0][0]) >= abs(_mn[1][1]-_mn[0][1]): # mn_vetical <= mn_horizontal
            return self.layers[metal.layer[0]].insert_metal(self.grids[metal.layer[0]].mn(metal), net_name=metal.netname, metal_ref=metal)
        else:
            self.layers[metal.layer[0]].deviants.append(metal)
            return None
        
    def merge(self):
        for layer in self.layers.values():
            layer.merge()
    
    def insert_instance_blackbox(self, inst):
        pin_list = list()
        for pin in inst.pins.values():
            # _mn = [[-9999,-9999],[9999,9999]]
            # _mn[0][0] = round(pin.xy[0][0]/72)
            # _mn[0][1] = round(pin.xy[0][1]/72)
            # _mn[1][0] = round(pin.xy[1][0]/72)
            # _mn[1][1] = round(pin.xy[1][1]/72)
            if pin.layer[0] in self.layers and pin.netname not in self.net_ignore:
                self.layers[pin.layer[0]].insert_metal(self.grids[pin.layer[0]].mn(pin), net_name=pin.netname)
                pin_list.append(pin)
        return pin_list
    
    def insert_instance_rect(self, inst:dict, Mxy, Mtransform, Mname, net_extern:dict=None):
        # loop invariant: the xy value of 'inst' is not inserted in the stack yet
        # start with insert the xy value of 'inst'
        _xy = np.array(inst['xy'])
        _transform = inst['transform']
        transform = tf.Mt(_transform) @ Mtransform
        libname = inst['libname']
        cellname = inst['cellname']
        if isinstance(self.libs, dict) and libname in self.libs:
            print(Mname, inst['cellname'], inst['name'], Mxy, Mtransform.tolist(), _transform, transform.tolist())
            lib = self.libs[libname]
            if cellname in lib.keys():
                cell = lib[cellname]
            else:
                print("Error! no such cell: "+cellname)
                return
            # inserting sub-blocks
            if 'sub_blocks' in cell:
                sub_blocks = cell['sub_blocks']
                if sub_blocks.__class__ is dict:     
                    if net_extern.__class__ is dict:    
                        for block in sub_blocks.values():
                            if 'pins' in block.keys():
                                _pins_extern = block['pins']
                                net_block = dict()
                                for pin in _pins_extern.values():
                                    _netname = pin['netname']
                                    if _netname in net_extern.keys():
                                        net_block[pin['termName']] = net_extern[_netname]
                                    else:
                                        net_block[pin['termName']] = pin['netname']
                            else:
                                net_block = None
                            xy = Mxy + _xy @ Mtransform
                            blockName = Mname+'/'+ block['name']
                            self.insert_instance_rect(block, xy, transform, blockName, net_block)
                    else:
                        for block in sub_blocks.values():
                            if 'pins' in block.keys():
                                _pins_extern = block['pins']
                                net_block = dict()
                                for pin in _pins_extern.values():
                                    net_block[pin['termName']] = pin['netname']
                            else:
                                net_block = None
                            xy = Mxy + _xy @ Mtransform
                            blockName = Mname+'/'+ block['name']
                            self.insert_instance_rect(block, xy, transform, blockName, net_block)
            # transform_child = R0 @ transform = transform
            # case Pin: just insert corresponding metal
            # -> can be deleted! B.C the corresponding metals are in port, the Pin itself is just for netname propagation
            # pins = cell['pins']
            # for _pinName, _pin in pins.items():
            #     if _pin['netname'] in self.net_ignore:
            #         continue
            #     _layerName = _pin['layer'][0]
            #     if _layerName in self.layers:
            #         offset = Mxy + _xy @ Mtransform
            #         rxy = offset + (np.array(_pin['xy']) @ Mtransform) # not ordered, but will be at insert_metal function
            #         _pin0 = laygo2.object.physical.Pin(xy=rxy, layer=_pin['layer'], netname = Mname+'__'+_pinName)
            #         _rmn = self.grids[_layerName].mn(_pin0)
            #         self.insert_metal(Rect(xy=self.grids[_layerName].xy(_rmn), layer=[_layerName,'drawing'], netname = Mname+'__'+_pinName))
            # deal with additional metals
            if 'ports' in cell:
    #            if net_extern.__class__ is dict: # sub-block with external pins
                for netName in cell['ports'].keys():
                    ports = cell['ports'][netName]
                    if net_extern[netName] in self.net_ignore:
                        continue
                    _netname_extern = net_extern[netName]
                    for metal in ports:
                        _layer = metal['layer']
                        _grid = self.grids[_layer]
                        offset = Mxy + _xy @ Mtransform
                        rxy = offset + (np.array(metal['xy']) @ transform) # not ordered, but will be at insert_metal function
                        if self.layers_orient[_layer] == 'horizontal':
                            _Rect = Rect(xy=rxy, vextension=int(_grid.hwidth[_grid.mn(rxy)[0][0]]), layer=[_layer,'drawing'], netname=_netname_extern)
                        else:
                            _Rect = Rect(xy=rxy, hextension=int(_grid.vwidth[_grid.mn(rxy)[0][0]]), layer=[_layer,'drawing'], netname=_netname_extern)
                        # print("layer:",_layer,"xy:",self.grids[_layer].mn(_Rect).tolist(),"internal:",netName, "extern:",_netname_extern)
                        self.insert_metal(_Rect)                    
                # else: # not sub-block or sub-block without external pins
                #     for pinName in cell['ports'].keys():
                #         ports = cell['ports'][pinName]
                #         if len(ports) < 1:
                #             continue
                #         elif ports[0]['netname'] in self.net_ignore:
                #             continue
                #         _netname = ports[0]['netname']
                #         for metal in ports:
                #             _layer = metal['layer']
                #             _grid = self.grids[_layer]
                #             offset = Mxy + _xy @ Mtransform
                #             rxy = offset + (np.array(metal['xy']) @ transform) # not ordered, but will be at insert_metal function
                #             if self.layers_orient[_layer] == 'horizontal':
                #                 _Rect = Rect(xy=rxy, vextension=int(_grid.hwidth[_grid.mn(rxy)[0][0]]), layer=[_layer,'drawing'], netname = Mname+'__'+_netname)
                #             else:
                #                 _Rect = Rect(xy=rxy, hextension=int(_grid.vwidth[_grid.mn(rxy)[0][0]]), layer=[_layer,'drawing'], netname = Mname+'__'+_netname)
                #             self.insert_metal(_Rect)
            # TODO: change assigning netname as __obtacle__ into assigning netname as masterName/instName/netname
            if 'obstacles' in cell:
                obstacles = cell['obstacles']
                for metal in obstacles:
                    _layer = metal['layer']
                    _grid = self.grids[_layer]
                    offset = Mxy + _xy @ Mtransform
                    rxy = offset + (np.array(metal['xy']) @ transform) # not ordered, but will be at insert_metal function
                    if self.layers_orient[_layer] == 'horizontal':
                        _Rect = Rect(xy=rxy, vextension=int(_grid.hwidth[_grid.mn(rxy)[0][0]]), layer=[_layer,'drawing'], netname = Mname+'__obstacle__')
                    else:
                        _Rect = Rect(xy=rxy, hextension=int(_grid.vwidth[_grid.mn(rxy)[0][0]]), layer=[_layer,'drawing'], netname = Mname+'__obstacle__')
                    self.insert_metal(_Rect)
        elif "_microtemplates_" in libname:
            return
        else:
            print("no such library: "+libname)
            return
        
    def insert_virtual_instance(self, vinst): 
        via_list=list()
        pin_list=list()
        mxy = vinst.xy
        mtf = vinst.transform
        if 'rtrunk_' in vinst.cellname: # special case1: routingmesh
            print("rtrunk found:",vinst.name)
            for velem in vinst.native_elements.values():
                if isinstance(velem, Rect): #insert metals of a virtual instance
                    _mn = self.grids[velem.layer[0]].mn(velem)
                    self.layers[velem.layer[0]].insert_metal(_mn, net_name=velem.netname, metal_ref=velem)
                elif isinstance(velem, Pin):
                    continue
                elif self.is_via(velem):
                    via_list.append(viaNode(self._via_table[velem.cellname], self.grids[ self._via_table[velem.cellname][0] ].mn(velem.xy)))
        else: # ordinary virtual instance
            for velem in vinst.native_elements.values():
                if isinstance(velem, Rect): #insert metals of a virtual instance
                    if velem.layer[0] in self.layers:
                        _xy = mxy + np.dot(np.sort(velem.xy, axis=0), tf.Mt(mtf).T)
                        _mn = self.grids[velem.layer[0]].mn(Rect(xy=_xy, layer=velem.layer)) # temporal Rect implemented
                        if _mn[0][0] is None or _mn[0][1] is None:
                            print("grid mismatch: ", velem.layer[0], self.grids[velem.layer[0]].name, mxy, mtf, (tf.Mt(mtf).T).tolist(), _xy.tolist())
                            continue
                        self.layers[velem.layer[0]].insert_metal(_mn, net_name=velem.netname, metal_ref=velem)
                elif isinstance(velem, Pin):
                    continue
                elif self.is_via(velem):
                    via_list.append(viaNode(self._via_table[velem.cellname], self.grids[ self._via_table[velem.cellname][0] ].mn(velem.xy)))
            #TODO: implement insert instance code
            for vpin in vinst.pins.values():
                if vpin.netname in self.net_ignore:
                    continue
                pin_list.append(vpin)
        return pin_list, via_list

    def insert_pin(self, pin):
        _layer_name = pin.layer[0]
        _layer = self.layers[_layer_name]
        # somehow grid.mn(pin.xy) dosen't work (it returns 'None' for numbers those are not multiple of 72)
        # _mn = [[-9999,-9999],[9999,9999]]
        # _mn[0][0] = round(pin.xy[0][0]/72)
        # _mn[0][1] = round(pin.xy[0][1]/72)
        # _mn[1][0] = round(pin.xy[1][0]/72)
        # _mn[1][1] = round(pin.xy[1][1]/72)
        _mn = self.grids[pin.layer[0]].mn(pin)
        # self.insert_metal(Rect(xy=self.grids[_layer_name].xy(_mn), layer=[_layer_name,'drawing'], netname = pin.netname))
        # print(_layer_name, _mn, pin.netname, pin.name)
        if _mn[0][0] > _mn[1][0]:
            _mn[0][0], _mn[1][0]=_mn[1][0], _mn[0][0]       
        if _mn[0][1] > _mn[1][1]:
            _mn[0][1], _mn[1][1]=_mn[1][1], _mn[0][1]

        if self.layers_orient[_layer_name] == "horizontal":
            for i in range(_mn[0][1], _mn[1][1]+1):
                _row_idx = self.layers[_layer_name].search_rc(i)
                if _row_idx is None: # error case
                    print("pin error: No metal on "+_layer_name+', y='+str(_mn[0][1]))
                    return None
                if _layer.rows[_row_idx].metal_list[0].mn[0][0] > _mn[0][0]:
                    print("pin error: No metal on "+_layer_name+', x='+str(_mn[0][0])+', y='+str(_mn[0][1]))
                    return None
                elif _layer.rows[_row_idx].metal_list[len(_layer.rows[_row_idx].metal_list)-1].mn[1][0] < _mn[1][0]:
                    print("pin error: No metal on "+_layer_name+', x='+str(_mn[1][0])+', y='+str(_mn[0][1]))
                    return None
                # invariant: node[0].x1 <= pinx1 <= pinx2 <= node[last].x2
                _pin_idx = _layer.search_metal_index(_row_idx, _mn[0][0])-1
                if _layer.rows[_row_idx].metal_list[_pin_idx].mn[1][0] < _mn[1][0]:
                    print("pin error: No metal on "+_layer+', x='+str(_mn[1][0])+', y='+str(_mn[0][1]))
                    return None
                self.pins.append((pin,_layer.rows[_row_idx].metal_list[_pin_idx]))
        else: # layer orient == vertical
            for i in range(_mn[0][0], _mn[1][0]+1):
                _col_idx = self.layers[_layer_name].search_rc(i)
                if _col_idx is None: # error case
                    print("pin error: No metal on "+_layer_name+', x= %d' % (_mn[0][1]))
                    return None
                if _layer.cols[_col_idx].metal_list[0].mn[0][1] > _mn[0][1]:
                    print("pin error: No metal on "+_layer_name+', x= %d, y= %d' % (_mn[0][0],_mn[0][1]))
                    return None
                elif _layer.cols[_col_idx].metal_list[len(_layer.cols[_col_idx].metal_list)-1].mn[1][1] < _mn[1][1]:
                    print("pin error: No metal on "+_layer_name+', x= %d, y= %d' % (_mn[1][0],+_mn[0][1]))
                    return None
                # invariant: node[0].y1 <= piny1 <= piny2 <= node[last].y2
                _pin_idx = _layer.search_metal_index(_col_idx, _mn[0][1])-1
                if _layer.cols[_col_idx].metal_list[_pin_idx].mn[1][1] < _mn[1][1]:
                    print("pin error: No metal on "+_layer_name+', x= %d, y= %d' % (_mn[1][0],_mn[1][1]))
                    return None
                self.pins.append((pin,_layer.cols[_col_idx].metal_list[_pin_idx]))
        return self.pins
    def is_via(self, _inst):
        if _inst.cellname in self._via_table:
            return True
        else: 
            return False

    def insert_via(self, via):
        # if not self.is_via(inst):
        #     print("error:"+inst+"is not_via")
        #via_mn, layer1(vertical), layer2(horizontal) mapping
        via_mn = via.mn
        if self.layers_orient[via.layer_pair[0]] == "vertical":
            layer1_name, layer2_name = via.layer_pair
        else:
            layer2_name, layer1_name = via.layer_pair
        layer1, layer2 = self.layers[layer1_name], self.layers[layer2_name]
        col_index = layer1.search_rc(via_mn[0])
        if col_index is None:
            print("Error: no metal on layer= "+layer1_name+" x="+str(via_mn[0]))
            return
        metal_index = layer1.search_metal_index(col_index, via_mn[1])-1
        if metal_index == -1:
            print("Error: no metal on via: layer="+layer1_name+" x="+str(via_mn[0]))
            return
        layer1_metal = layer1.cols[col_index].metal_list[metal_index]
        if layer1_metal.mn[1][1] < via_mn[1]:
            print("Error: no metal on via: layer="+layer1_name+" x="+str(via_mn[0])+" y="+str(via_mn[1]))
            return

        row_index = layer2.search_rc(via_mn[1])
        if row_index is None:
            print("Error: no metal on layer= "+layer2_name+" y="+str(via_mn[1]))
            return
        metal_index = layer2.search_metal_index(row_index, via_mn[0])-1
        if metal_index == -1:
            print("Error: no metal on via: layer= "+layer2_name+" y="+str(via_mn[1]))
            return
        layer2_metal = layer2.rows[row_index].metal_list[metal_index]
        if layer2_metal.mn[1][0] < via_mn[0]:
            print("Error: no metal on via: layer="+layer2_name+" x="+str(via_mn[0])+" y="+str(via_mn[1]))
            return
        
        layer1_metal.via_list.append([layer2_metal, layer2, via_mn])
        layer2_metal.via_list.append([layer1_metal, layer1, via_mn])
    
    def is_occupied(self, layer, mn):
        return self.layers[layer].is_occupied(mn)

    def check_node(self, ref_net_name, metal):
        _netname = metal.net_name
        if len(_netname) == 1 and ref_net_name not in _netname:
            if _netname == {None}:
                metal.net_name.remove(None)
                metal.net_name.add(ref_net_name)
            else:    
                print(_netname,end='')
                print(' is connected to',end=' ')
                print(set([ref_net_name]))
                print(metal.layer, metal.mn)
        elif len(_netname) == 2:
            if None in _netname and _netname - {None} == {ref_net_name}:
                metal.net_name.remove(None)
            else:
                print(_netname,end='')
                print(' are connected to',end=' ')
                print(set([ref_net_name]))
                print(metal.layer, metal.mn)
                metal.net_name.add(ref_net_name)
        elif len(_netname) > 2:
            print(_netname,end='')
            print(' are connected to',end=' ')
            print(set([ref_net_name]))
            print(metal.layer, metal.mn)
            metal.net_name.add(ref_net_name)
        else:
            pass

    def preview(self):
        libtemp = laygo2.object.database.Library(name="tempLib")
        dsntemp = laygo2.object.database.Design(name="preview")
        libtemp.append(dsntemp)
        for _layerName, _layer in self.layers.items():
            if self.layers_orient[_layerName] == 'vertical':
                for col in _layer.cols:
                    for metal in col.metal_list:
                        _xy = self.grids[metal.layer].abs2phy(metal.mn)
                        dsntemp.append(Rect(xy = _xy, layer = [_layerName,'drawing'], ))
            elif self.layers_orient[_layerName] == 'horizontal':
                for row in _layer.rows:
                    for metal in row.metal_list:
                        _xy = self.grids[metal.layer].abs2phy(metal.mn)
                        dsntemp.append(Rect(xy = _xy, layer = [_layerName,'drawing']))
            else:
                print("grid error")
        laygo2.interface.gds.export(libtemp, filename=None, cellname="preview", scale = 1e9, layermapfile="laygo2_example/prj_db/gds_sky130.layermap",
           physical_unit=1e-9, logical_unit=0.5, pin_label_height=0.1,
           svg_filename="netMap_preview.svg",png_filename="netMap_preview.png")
    
    def net_traverse(self, pin, pin_net_name,pin_set): # travel net graph in BFS order. Source node is pin
        # if pin.visited == True: # pin nodes also could be visited by previous traverse
        #     return
        if pin_net_name in pin_set:
            if pin.visited == False:
                print("open error: %s [(%d %d),(%d %d)] is not connected to same named net" %(pin_net_name,pin.mn[0][0],pin.mn[0][1],pin.mn[1][0],pin.mn[1][1]))
            #     pin.visited = True
            else:
            #     pass
            # #    print("Warning: %s is repeated but connected to same named net"%(pin_net_name))
                return
        else:
            pin_set.add(pin_net_name)
        queue = deque([pin])
        # ref_net = set()
        # ref_net.add(pin_net_name)
        # set the pin node visited
        pin.visited = True
        # repeat until the queue is empty
        print("start queue: ",pin_net_name)
        while queue:
            # pop v and check netname of v
            v = queue.popleft()
            # print(v.net_name, v.mn)
            self.check_node(pin_net_name, v)
            # insert nodes that connected to v and not visited 
            for node, layer, via_mn in v.via_list:
                if not node.visited:
                    queue.append(node)
                    node.visited = True

    def trimNone(self):
        for layer in self.layers.values():
            for mlist in layer.rc:
                for metal in mlist.metal_list:
                    if len(metal.net_name) > 1 and None in metal.net_name:
                        metal.net_name.remove(None)

    @classmethod
    def import_from_design(cls, dsn, grid_table, via_table, orient_first="vertical", layer_names=['M3','M4','M5'], net_ignore = [],lib_ref:str = None):
        if os.path.exists(lib_ref):
            libraries = dict()
            with open(lib_ref, 'r') as stream:
                db = yaml.load(stream, Loader=yaml.FullLoader)                
            for _libname, _lib in db['libraries'].items():
                filename = _lib['path']
                if os.path.exists(filename):
                    with open(filename, 'r') as stream:
                        buffer = yaml.load(stream, Loader=yaml.FullLoader)
                        libraries[_libname] = buffer[_libname]
                else:
                    print("library file not found:", _libname, filename)
        else:
            print("project file not found")
            libraries = None
        nMap = cls(via_table=via_table, grid_table = grid_table, orient_first=orient_first, layer_names=layer_names, net_ignore = net_ignore, libs = libraries)
        # import top cell layer elements(including rects & vias of virtual-instances)
        pin_list = list()
        pin_set = set()   
        for pin in dsn.pins.values():
            pin_list.append(pin)
            nMap.terminals.append(pin) # for exporting ports
        for vinst in dsn.virtual_instances.values():
            # special case1: routingMesh
            # special case2: via
            _pin_list, _via_list = nMap.insert_virtual_instance(vinst)
            pin_list.extend(_pin_list)
            nMap.vias.extend(_via_list)
        for inst in dsn.instances.values():
            if nMap.is_via(inst):
                nMap.vias.append(viaNode(nMap._via_table[inst.cellname], nMap.grids[via_table[inst.cellname][0]].mn(inst.xy)))
            else:
                pin_list.extend(nMap.insert_instance_blackbox(inst))
        for rect in dsn.rects.values():
            if rect.layer[0] in layer_names:
                _metal = nMap.insert_metal(rect)
            else:
                continue
        # DFS forest
        for _instName, _inst in dsn.instances.items():
            if _inst.shape is None:
                inst = dict()
                inst['xy'] = _inst.xy
                inst['transform'] = _inst.transform
                inst['libname'] = _inst.libname
                inst['cellname'] = _inst.cellname
                inst['name'] = _instName
                if inst['libname'] in libraries.keys() and libraries[inst['libname']][inst['cellname']]['pins'] is not None:
                    _pins = libraries[inst['libname']][inst['cellname']]['pins']
                    net_extern = dict()
                    for _pinName, _pin in _pins.items():
                        if _pin['netname'] is not None or _pin['netname'] is not 'None':
                            net_extern[_pin['netname']] = _inst.pins[_pinName].netname
                else:
                    net_extern = None
                print(net_extern)
                print("start implementing "+inst['name']+" xy: ",inst['xy'])
                nMap.insert_instance_rect(inst, np.array([0,0]), np.array([[1,0],[0,1]]), inst['name'], net_extern)
                print("finish implementing "+inst['name'])
            else:
                _transform = _inst.transform
                transform = tf.Mt(_transform)
                for col in range(_inst.shape[0]):
                    for row in range(_inst.shape[1]):
                        inst = dict()
                        inst['transform'] = _transform
                        inst['libname'] = _inst.libname
                        inst['cellname'] = _inst.cellname
                        inst['name']:str = _instName + '_'+str(col*_inst.shape[1]+row)
                        inst['xy'] = _inst.xy + np.array([_inst.pitch[0]*col,_inst.pitch[1]*row]) @ transform
                        if inst['libname'] in libraries.keys() and libraries[inst['libname']][inst['cellname']]['pins'] is not None:
                            _pins = libraries[inst['libname']][inst['cellname']]['pins']
                            net_extern = dict()
                            for _pinName, _pin in _pins.items():
                                if _pin['netname'] is not None or _pin['netname'] is not 'None':
                                    net_extern[_pin['netname']] = _inst.pins[_pinName].netname
                        else:
                            net_extern = None
                        print("start implementing "+inst['name']+" xy: ",inst['xy'])
                        nMap.insert_instance_rect(inst, np.array([0,0]), np.array([[1,0],[0,1]]), inst['name'], net_extern)
                        print("finish implementing "+inst['name'])
        #unmerged->merge
        nMap.merge()
        for via in nMap.vias:
            nMap.insert_via(via)
        _layers = set(layer_names)  
        for pin in pin_list:
            if pin.layer[0] in _layers and pin.netname not in nMap.net_ignore:
                if pin.netname is None:
                    print("error, pinname:", pin.name, "netname: None","layer:",pin.layer[0])
                    continue
                _pins = nMap.insert_pin(pin)

        # update netname (propagation)
        for pin, pin_node in nMap.pins:
        #    print("\npin name: %s, netname: %s, layer: %s, xy:[[%d, %d],[%d, %d]], pin_node.xy:"\
        #         % (pin.name, pin.netname,pin.layer[0],nMap.grids[pin.layer[0]].mn(pin)[0][0],nMap.grids[pin.layer[0]].mn(pin)[0][1],nMap.grids[pin.layer[0]].mn(pin)[1][0],nMap.grids[pin.layer[0]].mn(pin)[1][1]))#,end='')
            nMap.net_traverse(pin_node,pin.netname,pin_set)
        nMap.trimNone()
        for layer in nMap.layers.values():
            for mlist in layer.rc:
                for metal in mlist.metal_list:
                    print(metal.layer, metal.mn, metal.net_name)
                    if metal.metal_ref is None: # pin or merged metal
                        for segment in metal.metal_seg:
                            if segment.metal_ref is not None:
                                for _name in metal.net_name: # for accessing set type data
                                    segment.metal_ref.netname = _name
                    else: # unmerged metal
                        for _name in metal.net_name: # for accessing set type data
                            segment.metal_ref.netname = _name
        # # get rect set in dsn for adding port metals
        # _rects_dsn = set()
        # for rect in dsn.rects.values():
        #     _rects_dsn.add(rect)
        # _termNames = set()
        # _newName = "NoName_port_"
        # _idx = 0
        # for terminal in nMap.terminals:
        #     _termNames.add(terminal.netname)
        # for layer in nMap.layers.values():
        #     for mlist in layer.rc:
        #         for metal in mlist.metal_list:
        #             if metal.metal_ref is None: # pin or merged metal
        #                 for segment in metal.metal_seg:
        #                     if segment.metal_ref is not None:
        #                         for _name in metal.net_name: # for accessing set type data
        #                             if _name in _termNames and segment.metal_ref not in _rects_dsn:
        #                                 dsn.rects[_newName+str(_idx)] = segment.metal_ref
        #                                 _idx += 1
        #             else: # unmerged metal
        #                 for _name in metal.net_name: # for accessing set type data
        #                     if _name in _termNames and segment.metal_ref not in _rects_dsn:
        #                         dsn.rects[_newName+str(_idx)] = segment.metal_ref
        #                         _idx += 1   
        return nMap

    @classmethod
    # deprecated
    def lvs_check(cls, dsn, grid, via_table, orient_first="vertical", layer_names=['M1','M2','M3','M4','M5']):
        nMap = cls(grid=grid, via_table=via_table, orient_first=orient_first, layer_names=layer_names)
        pin_list = list()
        pin_set = set()
        for pin in dsn.pins.values():
            pin_list.append(pin)
        for vinst in dsn.virtual_instances.values():
            _pin_list, _via_list = nMap.insert_virtual_instance(vinst)
            pin_list.extend(_pin_list)
            nMap.vias.extend(_via_list)
        for inst in dsn.instances.values():
            if nMap.is_via(inst):
                nMap.vias.append(viaNode(nMap._via_table[inst.cellname], nMap.grid.mn(inst.xy)))
            else:
                pin_list.extend(nMap.insert_instance_blackbox(inst))
        for rect in dsn.rects.values():
            _metal = nMap.insert_metal(rect)
        #unmerged->merge
        nMap.merge()
        for via in nMap.vias:
            nMap.insert_via(via)
        for pin in pin_list:
            nMap.insert_pin(pin)
        #lvs test by bfs
        for pin, pin_node in nMap.pins:
            print("pin name: %s, netname: %s, layer: %s, xy:[[%d %d][%d %d]], pin_node.xy:"\
                % (pin.name, pin.netname,pin.layer[0],nMap.grid.mn(pin)[0][0],nMap.grid.mn(pin)[0][1],nMap.grid.mn(pin)[1][0],nMap.grid.mn(pin)[1][1]),end='')
            print(pin_node.mn)
            nMap.net_traverse(pin_node,pin.netname,pin_set)
        for layer in nMap.layers.values():
            for mlist in layer.rc:
                for metal in mlist.metal_list:
                    if metal.visited is not True:
                        print(metal.layer, metal.mn, metal.net_name)