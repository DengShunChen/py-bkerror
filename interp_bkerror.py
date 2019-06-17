#!/usr/bin/env python
#---documentation block 
#
# Purpose : the interpolation of GSI background error covariance  
#
# Modified by Deng-Shun Chen
#
# log : 
#   2018-12-04  D-S.Chen  Created
#   2019-01-14  D-S.Chen  add vertical interpolting capability  
#   2019-01-17  D-S.Chen  change 2d interpolation to "linear" instead of "cubic".
#                         the cubic scheme cause negitive values. 
#---documentation block 

import sys
import numpy as np
from scipy.interpolate import interp1d, interp2d
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from bkerror import bkerror
from splat import splat
from scipy import array
from matplotlib import pyplot as plt

class GSIbkgerr(object):
    '''
    Object containing GSI static background error information
    '''
    def __init__(self,filename):
        '''
        Read and store GSI background error file.
        '''
        nsig,nlat,nlon = bkerror.get_header(filename)
        ivar,agvin,bgvin,wgvin,corzin,hscalesin,vscalesin,corq2in,corsstin,hsstin,corpin,hscalespin = bkerror.get_bkerror(filename,nsig,nlat,nlon)
        var = (ivar.tostring()).replace('\x00','')[:-1].split('|')

        self.filename = filename

        self.nsig = nsig
        self.nlat = nlat
        self.nlon = nlon

        self.ivar = ivar
        self.var = var

        self.agvin = agvin
        self.bgvin = bgvin
        self.wgvin = wgvin
        self.corzin = corzin
        self.hscalesin = hscalesin
        self.vscalesin = vscalesin
        self.corq2in = corq2in
        self.corsstin = corsstin
        self.hsstin = hsstin
        self.corpin = corpin
        self.hscalespin = hscalespin

        return

    def print_summary(self):
        '''
        Print a summary of the GSI background error file
        '''

        print
        print 'file = %s' % self.filename
        print 'nsig = %d, nlat = %d, nlon = %d, nvar = %d' % (self.nsig,self.nlat,self.nlon,len(self.var))
        print 'variables = %s' % ', '.join(self.var)
        print 'agv.shape: ', self.agvin.shape
        print 'bgv.shape: ', self.bgvin.shape
        print 'wgv.shape: ', self.wgvin.shape
        print 'corz.shape: ', self.corzin.shape
        print 'hscales.shape: ', self.hscalesin.shape
        print 'vscales.shape: ', self.vscalesin.shape
        print 'corq2.shape: ', self.corq2in.shape
        print 'corsst.shape: ', self.corsstin.shape
        print 'hsst.shape: ', self.hsstin.shape
        print 'corp.shape: ', self.corpin.shape
        print 'hscalesp.shape: ', self.hscalespin.shape
        print

        return

def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(map(pointwise, array(xs)))

    return ufunclike

def read_siglevel(filename):

  f = open(filename,'r')
  sis=[]
  nsig = float(f.readline())
  for line in f:
    line = line.strip().split()
    sis.append(float(line[0]))
  sis.insert(0,1.)
  sis.append(0.)
  sigma = np.arange(nsig)
  for k,sig in enumerate(sigma):
    sigma[k] = 0.5*(sis[k]+sis[k+1])
  f.close()

  sigma = np.array(sigma,dtype=np.float32)

  return nsig,sigma

def berror_interp2d(x,y,berr,x_n,y_n,sig_flip=1,positive=0):
 
  # cause sigma values are not ascending, should be fliped. 
  if sig_flip == 1:
    x = np.flip(x,0)
    x_n = np.flip(x_n,0)
    berr = np.flip(berr,0)

  # 2d interpolation
  interp_kind='linear'
  f = interp2d(x,y,berr,kind=interp_kind)
  tmp_n = f(x_n,y_n)
  berr_n = np.array(tmp_n,dtype=np.float32)

  diff_max = np.abs(np.amax(berr_n) - np.amax(berr))
  diff_min = np.abs(np.amin(berr_n) - np.amin(berr))
  
  if diff_max > 1 or diff_min > 1:
    print('Warnning!! errors are too large ...')
    print('berr   max/min',np.amax(berr), np.amin(berr))
    print('berr_n max/min',np.amax(berr_n), np.amin(berr_n))
    
  # flip back to original order
  if sig_flip == 1:
    berr_n = np.flip(berr_n,0)

  # positive defined 
  if positive == 1:
    berr_n[berr_n<0] = 0.

  return berr_n

def berror_interp1d(x,berr,x_n):

  f_i = interp1d(x,berr,kind=interp_kind)
  f = extrap1d(f_i)
  tmp_n = f(x_n)
  berr_n = np.array(tmp_n,dtype=np.float32)

  return berr_n

def plot_comp(x,y,berr,x_n,y_n,berr_n,var_str='berr'):
  cmapdiv = 'Spectral_r'
  plt.figure()
  plt.subplot(2,1,1)
  x1,y1 = np.meshgrid(x,y)
  plt.contourf(x1,y1,berr,21,vmin=-berr.max(),cmap=cmapdiv,extend='both')
  plt.colorbar()
  plt.title('%s before' % (var_str),fontsize=12,fontweight='normal')
  plt.subplots_adjust(bottom=0.2)

  plt.subplot(2,1,2)
  x2,y2 = np.meshgrid(x_n,y_n)
  plt.contourf(x2,y2,berr_n,21,vmin=-berr_n.max(),cmap=cmapdiv,extend='both')
  plt.colorbar()
  plt.title('%s after' % (var_str),fontsize=12,fontweight='normal')

  plt.show()


if __name__ == '__main__': 
  # bkerror file to read; e.g. global_berror.l64y258.f77
  parser = ArgumentParser(description='read and interpolate background error file',formatter_class=ArgumentDefaultsHelpFormatter)
  parser.add_argument('--filename',help='background error file to read',type=str,required=True)
  parser.add_argument('--imax',help='interpolate to imax',type=int,required=True)
  parser.add_argument('--jmax',help='interpolate to jmax',type=int,required=True)
  parser.add_argument('--nsig',help='interpolate to nsig',type=int,required=False,default=72)
  parser.add_argument('--interpkind',help='interpolation kind ',type=str,required=False,default='cubic',choices=['linear','cubic'])
  args = parser.parse_args()

  jmax = args.jmax
  imax = args.imax
  nsig = args.nsig
  interp_kind = args.interpkind

  # Read and print some info about the file we are reading from
  gsi = GSIbkgerr(args.filename)
  gsi.print_summary()

  if imax is None or jmax is None:
    sys.exit(0)

  # Read in the same file for interpolation
  gsi_n = GSIbkgerr(args.filename)
  gsi_n.filename = 'berror_stats'

  Ps=1000.
  Pt=0.1
  Pt_ncep = 0.2  # have to check this 

  # Read in siglevel infor 
  sigfile_cwb='global_siglevel.l72.txt'
  nsig_cwb, sigs_cwb = read_siglevel(sigfile_cwb)
  pres_cwb = sigs_cwb * (Ps-Pt) + Pt
 
  sigfile_ncep='global_siglevel.l64.txt'
  nsig_ncep, sigs_ncep = read_siglevel(sigfile_ncep)
  pres_ncep = sigs_ncep * (Ps-Pt_ncep) + Pt_ncep

  # set up old grid dimension
  idrt = 4
  glat,wlat = splat(idrt,gsi.nlat)
  slon = np.linspace(0.,360.,gsi.nlon,endpoint=False)
  slat = 180. / np.arccos(-1.) * np.arcsin(glat[::-1])

  # set up new grid dimension
  gsi_n.nlon = imax
  gsi_n.nlat = jmax
  gsi_n.nsig = nsig

  glat_n,wlat_n = splat(idrt,gsi_n.nlat)
  slon_n = np.linspace(0.,360.,gsi_n.nlon,endpoint=False)
  slat_n = 180. / np.arccos(-1.) * np.arcsin(glat_n[::-1])

  if gsi.nsig != nsig_ncep:
    raise ValueError('levels inconsistant! gsi.nsig=%d, ngsi_ncep=%d' % (gsi.nsig,nsig_ncep))
  if gsi_n.nsig != nsig_cwb:
    raise ValueError('levels inconsistant! gsi_n.nsig=%d, ngsi_cwb=%d' % (gsi_n.nsig,nsig_cwb))
 
  case = 1
  if case == 1 :
    ssig = sigs_ncep
    ssig_n = sigs_cwb
  elif case == 2:
    ssig = pres_ncep
    ssig_n = pres_cwb
  else:
    ssig = np.arange(gsi.nsig)
    ssig_n = np.arange(gsi_n.nsig)
   
 
  print 'Interpolate from %d to %d' % (gsi.nlat, gsi_n.nlat)

  # agvin
  sigmax = max(gsi_n.nsig,gsi.nsig)
  tmp_n = np.zeros((gsi_n.nlat,sigmax,sigmax))

  for isig in np.arange(gsi.nsig):
    tmp_n[:,:,isig] = berror_interp2d(ssig,slat,gsi.agvin[:,:,isig],ssig_n,slat_n)

  for isig_n in np.arange(gsi_n.nsig):
    tmp_n[:,isig_n,:] = berror_interp2d(ssig,slat_n,tmp_n[:,isig_n,:gsi.nsig],ssig_n,slat_n)

  gsi_n.agvin = np.array(tmp_n.reshape(gsi_n.nlat,gsi_n.nsig,gsi_n.nsig),dtype=np.float32)

  # bgvin
  gsi_n.bgvin = berror_interp2d(ssig,slat,gsi.bgvin,ssig_n,slat_n) 
 #plot_comp(ssig,slat,gsi.bgvin,ssig_n,slat_n,gsi_n.bgvin,'bgvin')

  # wgvin
  gsi_n.wgvin = berror_interp2d(ssig,slat,gsi.wgvin,ssig_n,slat_n) 
 #plot_comp(ssig,slat,gsi.wgvin,ssig_n,slat_n,gsi_n.wgvin,'wgvin')

  # corzin
  tmp_n = np.zeros((gsi_n.nlat,gsi_n.nsig,6))
  for ivar in np.arange(6):
    tmp_n[:,:,ivar] = berror_interp2d(ssig,slat,gsi.corzin[:,:,ivar],ssig_n,slat_n,positive=0)
  gsi_n.corzin = np.array(tmp_n.reshape(gsi_n.nlat,gsi_n.nsig,6),dtype=np.float32)

 #for ivar in np.arange(6):
 #  plot_comp(ssig,slat,gsi.corzin[:,:,ivar],ssig_n,slat_n,gsi_n.corzin[:,:,ivar],'corzin')

  # hscalesin
  tmp_n = np.zeros((gsi_n.nlat,gsi_n.nsig,6))
  for ivar in np.arange(6):
    tmp_n[:,:,ivar] = berror_interp2d(ssig,slat,gsi.hscalesin[:,:,ivar],ssig_n,slat_n,positive=0)
  gsi_n.hscalesin = np.array(tmp_n.reshape(gsi_n.nlat,gsi_n.nsig,6),dtype=np.float32)

 #for ivar in np.arange(6):
 #  plot_comp(ssig,slat,gsi.hscalesin[:,:,ivar],ssig_n,slat_n,gsi_n.hscalesin[:,:,ivar],'hscalesin')

  # vscalesin
  tmp_n = np.zeros((gsi_n.nlat,gsi_n.nsig,6))
  for ivar in np.arange(6):
    tmp_n[:,:,ivar] = berror_interp2d(ssig,slat,gsi.vscalesin[:,:,ivar],ssig_n,slat_n,positive=0)
  gsi_n.vscalesin = np.array(tmp_n.reshape(gsi_n.nlat,gsi_n.nsig,6),dtype=np.float32)
  
 #for ivar in np.arange(6):
 #  plot_comp(ssig,slat,gsi.vscalesin[:,:,ivar],ssig_n,slat_n,gsi_n.vscalesin[:,:,ivar],'vscalesin')

  # corq2in
  nlat = 25
  gsi_n.corq2in = np.zeros((gsi_n.nlat,gsi_n.nsig))
  gsi_n.corq2in[0:nlat,:] = berror_interp2d(ssig,slat,gsi.corq2in,ssig_n,slat)[0:nlat,:]
  #plot_comp(ssig,slat[0:110],gsi.corq2in[0:110],ssig_n,slat_n[0:110],gsi_n.corq2in[0:100],'corq2in')

  # corsstin
  gsi_n.corsstin = berror_interp2d(slon,slat,gsi.corsstin,slon_n,slat_n,sig_flip=0)
 #plot_comp(slon,slat,gsi.corsstin,slon_n,slat_n,gsi_n.corsstin,'corsstin')

  # hsstin
  gsi_n.hsstin = berror_interp2d(slon,slat,gsi.hsstin,slon_n,slat_n,sig_flip=0)
 #plot_comp(slon,slat,gsi.hsstin,slon_n,slat_n,gsi_n.hsstin,'hsstin')

  # corpin
  gsi_n.corpin = berror_interp1d(slat,gsi.corpin,slat_n)

  # hscalespin
  gsi_n.hscalespin = berror_interp1d(slat,gsi.hscalespin,slat_n) 

  # Print some info about the interpolated data
  gsi_n.print_summary()

  bkerror.put_bkerror(gsi_n.filename,gsi_n.ivar,\
        gsi_n.agvin,gsi_n.bgvin,gsi_n.wgvin,\
        gsi_n.corzin,gsi_n.hscalesin,gsi_n.vscalesin,\
        gsi_n.corq2in,gsi_n.corsstin,gsi_n.hsstin,gsi_n.corpin,gsi_n.hscalespin)

  # check data in berror_stats
  gsi_rn = GSIbkgerr(gsi_n.filename)
  gsi_rn.print_summary()

  print 'differences'
  print np.abs(gsi_n.agvin-gsi_rn.agvin).max()
  print np.abs(gsi_n.bgvin-gsi_rn.bgvin).max()
  print np.abs(gsi_n.wgvin-gsi_rn.wgvin).max()
  print np.abs(gsi_n.hscalesin-gsi_rn.hscalesin).max()
