#!/usr/bin/env python3
"""Reproduce the constructed mathematical controls reported for MRTOE.

This standalone subset tests only mathematically defined cooling curves.
It does not access original experimental records or archived laboratory software.
"""
from __future__ import annotations
import json
import argparse
from pathlib import Path
import numpy as np

TAU, EPS, K = 60., .05, .31
TH = K * EPS * np.sqrt(TAU)


def trap(y, t):
 return float(np.trapezoid(y,t))

def calc(x,t,mode='zero',gradient='mean',interior_only=False):
 x=np.asarray(x,dtype=float);t=np.asarray(t,dtype=float)
 w=min(int(TAU*10),len(x)//3)
 if w<2: w=2
 if w%2==0: w+=1
 kern=np.ones(w)/w; m=np.convolve(x,kern,'same')
 if mode=='normalized': m=m/np.convolve(np.ones(len(x)),kern,'same')
 if gradient=='mean': gate=np.abs(np.gradient(x,np.mean(np.diff(t))))>EPS
 elif gradient=='actual': gate=np.abs(np.gradient(x,t))>EPS
 else: raise ValueError(gradient)
 y=(x-m)*gate/TAU
 p=w//2; bd=np.zeros(len(x),dtype=bool)
 bd[:p]=True; bd[-p:]=True
 # For interior_only, integrate only the contiguous full-window subset
 # using the same trapezoid rule (does not cross the excluded boundary).
 masked_total=trap(y[~bd],t[~bd]) if (~bd).sum()>=2 else None
 return dict(score=trap(y,t),boundary=trap(y*bd,t),interior=trap(y*(~bd),t),
             interior_integral=masked_total,window=w,gate=int(sum(gate)),duration=float(t[-1]-t[0]),
             decision='ACTIVE' if trap(y,t)>TH else 'INACTIVE')

def main():
 parser = argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--output', type=Path, default=None,
                     help='Save the constructed-control results to a JSON file.')
 args = parser.parse_args()
 n=600; t=np.arange(n)*.1; end=t[-1]; quadratic=50.-25.*(t/end)**2
 curves={'linear_cooling':np.linspace(50,25,n), 'exponential_cooling':25+25*np.exp(-np.arange(n)/90), 'concave_quadratic_cooling':quadratic, 'constant_temperature':np.full(n,26.)}
 synthetic={}
 for name,x in curves.items():
  z=calc(x,t); normalized=calc(x,t,'normalized'); synthetic[name]=dict(original=z,normalized=normalized,strictly_decreasing=bool(np.all(np.diff(x)<0)),
                                                        y_start=float(x[0]),y_end=float(x[-1]))
  assert abs(z['score']-z['boundary']-z['interior'])<1e-10
  assert abs(normalized['score']-normalized['boundary']-normalized['interior'])<1e-10
 assert abs(synthetic['linear_cooling']['original']['score']-3.1094527363184)<1e-10
 assert abs(synthetic['exponential_cooling']['original']['score']-1.7866289632213)<1e-10
 assert abs(synthetic['exponential_cooling']['normalized']['score']-0.204860504465693)<1e-10
 # Exact analytic interior residual for T0 - b*t^2
 w=201;p=100; b=25/end**2; mu2=(p*(p+1)/3)*.1**2
 analytic=(b*mu2)*(t[-p-1]-t[p])/TAU
 numerical=synthetic['concave_quadratic_cooling']['normalized']['interior_integral']
 assert abs(analytic-numerical)<1e-11, (analytic,numerical)
 assert numerical>TH and synthetic['concave_quadratic_cooling']['strictly_decreasing']
 # this is a numerical property of the score on constructed input, not physical proof
 out=dict(threshold=TH,synthetic=synthetic,analytic_interior_quadratic={'b':b,'mu2_s2':mu2,'score':analytic,'threshold_exceeded':bool(analytic>TH)},original_records={})
 out['source_provenance'] = 'Synthetic only'
 text=json.dumps(out,indent=2)
 if args.output:
  args.output.parent.mkdir(parents=True, exist_ok=True)
  args.output.write_text(text+'\n',encoding='utf-8')
 print(text)

if __name__=='__main__':
 main()
