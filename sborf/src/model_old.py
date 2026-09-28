#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  RRN_attention.py
#
#  Copyright 2019 Gabriel Orlando <orlando.gabriele89@gmail.com>
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
import torch
from torch import nn
import torch as t
from torch.nn.utils.rnn import pack_padded_sequence, pad_packed_sequence


class attention(nn.Module):
    def __init__(self, nfea, hinner=8, heads=4, dev=t.device('cuda'), dropout=0.1):

        super(attention, self).__init__()

        self.rnn = nn.LSTM(input_size=nfea, hidden_size=int(heads / 2), bidirectional=True, batch_first=True,
                           num_layers=2, dropout=dropout).to(device=dev)

        self.batchnorm2 = nn.BatchNorm1d(hinner).to(device=dev)
        self.dev = dev
        self.dropout = nn.Dropout(p=dropout)
        self.softmax = nn.Softmax(dim=1)

        self.nonlinear = nn.ReLU().to(device=dev)

        self.final_linear = nn.Sequential(
            nn.Linear(heads, heads, bias=True).to(device=dev),
            nn.LayerNorm(heads),
            nn.Tanh().to(device=dev),
            nn.Linear(heads, 1, bias=True).to(device=dev)
        ).to(device=dev)
        self.sig = nn.Sigmoid().to(device=dev)

    # self.apply(self.init_weights)

    def init_weights(self, m):
        if isinstance(m, t.nn.Conv2d) or isinstance(m, t.nn.Conv1d) or isinstance(m, t.nn.Linear) or isinstance(m,
                                                                                                                t.nn.Bilinear):
            print("Initializing weights...", m.__class__.__name__)
            t.nn.init.xavier_uniform_(m.weight)
            if type(m.bias) != type(None):
                m.bias.data.fill_(0.01)
        elif isinstance(m, t.nn.Embedding):
            print("Initializing weights...", m.__class__.__name__)
            t.nn.init.xavier_uniform_(m.weight)

    def forward(self, x, xlens):

        # assert len(x.shape)==3
        if type(xlens) == type(None):
            xlens = torch.Tensor([x.shape[1]] * x.shape[0]).type(torch.LongTensor)

        # x=pack_padded_sequence(x,xlens,batch_first=True,enforce_sorted=False)
        x, h = self.rnn(x)
        x = pad_packed_sequence(x, batch_first=True)[0]
        x = self.nonlinear(x)
        mask = masking(x, xlens).to(device=self.dev)
        out = x.masked_fill(mask, 0)

        out1 = out.contiguous().view((out.shape[0] * out.shape[1], out.shape[2]))
        out1 = self.final_linear(out1)
        out1 = out1.view(out.shape[0], out.shape[1])
        mask_soft = (torch.arange(out1.shape[1])[None, :] < xlens[:, None]).to(device=self.dev)
        out1 = out1.masked_fill(~mask_soft, -float("inf"))
        out1 = self.softmax(out1)
        return out1


class prediction(nn.Module):
    def __init__(self, nfea, hinner=10, out=8, kernel_size=5, dev=t.device('cpu'), dropout=0.1):
        super(prediction, self).__init__()
        self.rnn = nn.LSTM(input_size=nfea, hidden_size=int(out / 2), bidirectional=True, batch_first=True,
                           num_layers=2, dropout=dropout).to(device=dev)
        self.batchnorm2 = nn.BatchNorm1d(hinner).to(device=dev)
        self.dev = dev
        self.dropout = nn.Dropout(p=dropout)

        self.nonlinear = nn.ReLU().to(device=dev)

        self.final_linear = nn.Sequential(
            nn.Linear(out, out, bias=True).to(device=dev),
            nn.LayerNorm(out),
            nn.Tanh().to(device=dev),
            nn.Linear(out, 1, bias=True).to(device=dev),

        ).to(device=dev)
        self.sig = nn.Sigmoid().to(device=dev)

    def init_weights(self, m):
        if isinstance(m, t.nn.Conv2d) or isinstance(m, t.nn.Conv1d) or isinstance(m, t.nn.Linear) or isinstance(m,
                                                                                                                t.nn.Bilinear):
            print("Initializing weights...", m.__class__.__name__)
            t.nn.init.normal_(m.weight)
            if type(m.bias) != type(None):
                m.bias.data.fill_(0.01)
        elif isinstance(m, t.nn.Embedding):
            print("Initializing weights...", m.__class__.__name__)
            t.nn.init.xavier_uniform_(m.weight)

    def forward(self, x, xlens):

        if type(xlens) == type(None):
            xlens = torch.Tensor([x.shape[1]] * x.shape[0]).type(torch.LongTensor)
        x, h = self.rnn(x)
        x = pad_packed_sequence(x, batch_first=True)[0]
        x = self.nonlinear(x)
        mask = masking(x, xlens).to(device=self.dev)
        out = x.masked_fill(mask, 0)
        out1 = out.contiguous().view((out.shape[0] * out.shape[1], out.shape[2]))
        out1 = self.final_linear(out1)
        out1 = out1.view(out.shape[0], out.shape[1])
        return out1


class selfatt_RRN(nn.Module):

    def __init__(self, nfea, hinner=100, out=100, heads=20, initial_channels=10, embedding_size=25, dev=t.device('cpu'),
                 compress_initial_features=True, name="noname"):
        super(selfatt_RRN, self).__init__()
        hinner = 10
        out = 20
        heads = 20
        dropout = 0.0
        self.compress_initial_features = compress_initial_features
        self.nonlinear = nn.ReLU
        self.name = name
        self.pred = prediction(initial_channels, hinner=hinner, out=out, dev=dev, dropout=dropout)
        self.att = attention(initial_channels, hinner=hinner, heads=heads, dev=dev, dropout=dropout)
        self.dev = dev
        self.device = dev
        self.name = 'self att'

        if compress_initial_features:
            self.mapping_fea = nn.Sequential(
                nn.Linear(nfea, initial_channels),
                nn.Dropout(p=dropout),
                self.nonlinear(),
                nn.Linear(initial_channels, initial_channels),
                self.nonlinear()
            ).to(device=dev)
            self.pred = prediction(initial_channels, hinner=hinner, out=out, dev=dev)
            self.att = attention(initial_channels, hinner=hinner, heads=heads, dev=dev)
        else:

            self.pred = prediction(nfea, hinner=hinner, out=out, dev=dev, dropout=dropout)
            self.att = attention(nfea, hinner=hinner, heads=heads, dev=dev, dropout=dropout)

        self.apply(self.init_weights)

    def init_weights(self, m):

        if isinstance(m, t.nn.Linear) or isinstance(m, t.nn.Conv1d) or isinstance(m, t.nn.Linear) or isinstance(m,
                                                                                                                t.nn.Bilinear):
            #print("Initializing weights...", m.__class__.__name__)
            t.nn.init.normal_(m.weight)
            if type(m.bias) != type(None):
                m.bias.data.fill_(0.01)
        elif isinstance(m, t.nn.Embedding):
            #print("Initializing weights...", m.__class__.__name__)
            t.nn.init.xavier_uniform_(m.weight)

    def forward(self, x, xlens=None, return_attentionPrediction=False):

        x = self.mapping_fea(x)

        if type(xlens) == type(None):
            xlens = torch.Tensor([x.shape[1]] * x.shape[0]).type(torch.LongTensor)

        x = pack_padded_sequence(x, xlens, batch_first=True, enforce_sorted=False)
        p = self.pred(x, xlens)
        a = self.att(x, xlens)
        if return_attentionPrediction:
            return a, p
        else:
            out = torch.bmm(torch.unsqueeze(p, 1), torch.unsqueeze(a, -1)).squeeze()
        return out


def masking(X, X_len):
    maxlen = X.size(1)
    mask = torch.arange(maxlen)[None, :] < X_len[:, None]

    mask = torch.unsqueeze(mask, 2)

    mask = mask.expand(-1, -1, X.shape[2])

    return ~mask




