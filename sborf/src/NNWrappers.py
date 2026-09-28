#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
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

import os, copy, random, time
from sys import stdout
import numpy as np
import torch as t
import torch

from torch.nn.utils.rnn import pad_sequence
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import QuantileTransformer
from sborf.src.model_old import selfatt_RRN
from sborf.src.dna_features import onehot_codon

import warnings
warnings.filterwarnings("ignore")

class NNwrapper():

	def __init__(self, device='cpu',load=False):
		self.model = selfatt_RRN(nfea=61, dev = device)
		self.device = device
		self.explain_scaler = None
		if load:
			self.load(load)

	class myDataset(Dataset):
		def __init__(self, X, Y = None):
			self.X = X
			self.Y = Y

		def __len__(self):
			return len(self.X)

		def __getitem__(self, idx):
			if self.Y is not None:
				return self.X[idx], self.Y[idx]
			else:
				return self.X[idx]

	def load(self,organism="cerevisiae"):
		path = os.path.dirname(os.path.abspath(__file__))+"/../models/"+organism+"/"
		try:
			assert os.path.isdir(path)
		except:
			raise(ValueError,"path need to be a folder containing the stateDict of the NN (named stateDiz.m)  and the scaler (named probabilityScaler.m)")

		self.model.load_state_dict(torch.load(path+"stateDiz.m",map_location="cuda" if "cuda" in self.device else "cpu"))
		self.probability_scaler = torch.load(path+"probabilityScaler.m",weights_only=False)

		explain_scaler_path = path+"explainScaler.m"
		if os.path.isfile(explain_scaler_path):
			self.explain_scaler = torch.load(explain_scaler_path,weights_only=False)
		else:
			self.explain_scaler = None

	def collateFunction(self,batch):

		x1=[]
		x2=[]

		for i in range(len(batch)):
			x1 += [batch[i][0]]
			x2 += [batch[i][1]]

		return x1,x2

	def collateFunctionPredict(self,batch):
		x1=[]
		for i in range(len(batch)):
			x1 += [batch[i]]
		return x1

	def get_params(self,deep):
		return {}
			
	def predict(self,X,scale_output=True):

		tokens = []
		for k in X:
			tokens+= [torch.tensor([onehot_codon[i] for i in k],device=self.device,dtype=torch.float)]

		dataset = self.myDataset(tokens)
		loader = DataLoader(dataset, batch_size=150, shuffle=False, sampler=None, num_workers=0,collate_fn=self.collateFunctionPredict)

		self.model.training = False
		self.model.eval()

		pred = []
		for xOrig in loader:

			padding_lens = []
			for i in xOrig:
				padding_lens += [len(i)]

			x = pad_sequence(xOrig, padding_value=0, batch_first=True).to(self.device)
			padding_lens = torch.tensor(padding_lens, device="cpu", dtype=torch.long)

			yp = self.model(x, padding_lens)

			pred += yp.data.cpu().numpy().flatten().tolist()

		if scale_output:
			pred = self.probability_scaler.transform(np.array(pred).reshape((-1,1))).reshape(-1)
		return pred

	def explain(self,X,scale_output=True):
		# scale_output=True returns, for every codon, a value in [-1, 1]: the organism-specific
		# quantile transformation of the raw score rescaled from [0, 1] (0 = median contribution)

		tokens = []
		for k in X:
			tokens+= [torch.tensor([onehot_codon[i] for i in k],device=self.device,dtype=torch.float)]

		dataset = self.myDataset(tokens)
		loader = DataLoader(dataset, batch_size=150, shuffle=False, sampler=None, num_workers=0,collate_fn=self.collateFunctionPredict)

		self.model.training = False
		self.model.eval()

		pred = []
		for xOrig in loader:

			padding_lens = []
			for i in xOrig:
				padding_lens += [len(i)]

			x = pad_sequence(xOrig, padding_value=0, batch_first=True).to(self.device)
			padding_lens = torch.tensor(padding_lens, device="cpu", dtype=torch.long)

			a,p = self.model(x, padding_lens,return_attentionPrediction=True)
			residue_scores = (a*p).data.cpu().numpy()
			for i in range(len(xOrig)):
				pred += [residue_scores[i,:len(xOrig[i])].tolist()]

		if scale_output:
			if self.explain_scaler is None:
				raise ValueError("this organism has no explainScaler.m in its models folder: cannot scale the explanation")
			pred = [(self.explain_scaler.transform(np.array(s).reshape((-1,1))).reshape(-1)*2-1).tolist() for s in pred]

		return pred



