#!/usr/bin/env python

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from bkerror import bkerror
import numpy as np

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
        print(ivar.tostring())
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
        print 'agv.shape/max/min: ', self.agvin.shape, np.amax(self.agvin), np.amin(self.agvin)
        print 'bgv.shape/max/min: ', self.bgvin.shape, np.amax(self.bgvin), np.amin(self.bgvin)
        print 'wgv.shape/max/min: ', self.wgvin.shape, np.amax(self.wgvin), np.amin(self.wgvin)
        print 'corz.shape/max/min: ', self.corzin.shape, np.amax(self.corzin), np.amin(self.corzin)
        print 'hscales.shape/max/min: ', self.hscalesin.shape, np.amax(self.hscalesin), np.amin(self.hscalesin)
        print 'vscales.shape/max/min: ', self.vscalesin.shape, np.amax(self.vscalesin), np.amin(self.vscalesin)
        print 'corq2.shape/max/min: ', self.corq2in.shape, np.amax(self.corq2in), np.amin(self.corq2in)
        print 'corsst.shape/max/min: ', self.corsstin.shape, np.amax(self.corsstin), np.amin(self.corsstin)
        print 'hsst.shape/max/min: ', self.hsstin.shape, np.amax(self.hsstin), np.amin(self.hsstin)
        print 'corp.shape/max/min: ', self.corpin.shape, np.amax(self.corpin), np.amin(self.corpin)
        print 'hscalesp.shape/max/min: ', self.hscalespin.shape, np.amax(self.hscalespin), np.amin(self.hscalespin)
        print

        return


# bkerror file to read; e.g. global_berror.l64y258.f77
parser = ArgumentParser(description='read background error file',formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-f','--filename',help='background error file to read',type=str,required=True)
args = parser.parse_args()

gsi = GSIbkgerr(args.filename)
gsi.print_summary()

sys.exit(0)
