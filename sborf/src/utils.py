#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  parser.py
#
#  Copyright 2019 u0127326 <u0127326@gbw-l-l0039>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
import numpy as np
import math
def buildTrainingData(diz):
    stop_codons = ["TAG","TAA","TGA"]
    x = []
    y = []
    prot_names = []
    skip = 0
    for i in sorted(list(diz.keys())):
        if diz[i]['abn_new'] == 'NA' or float(diz[i]['abn_new']) < 0.001:
            continue
        yt = [math.log10(float(diz[i]['abn_new']))]
        try:
            assert  len(diz[i]['dna'])%3==0
            for sc in stop_codons:
                assert not sc in diz[i]['dna']
            xt = [diz[i]['dna']]
        except:
            skip += 1
            continue
        x += xt
        y += yt
        prot_names += [i]

    # scaling  ###
    y = np.array(y)
    y_std = (y - y.min()) / (y.max() - y.min())
    y = y_std * (1 - 0) + 0
    return x,y,prot_names