#!/bin/bash

 ./interp_bkerror.py --filename global_berror.l64y1154.f77 --imax 2576 --jmax 1282 
 cp berror_stats global_berror.l72y1282.f77.tco 
 ./interp_bkerror.py --filename global_berror.l64y1154.f77 --imax 1552 --jmax 770 
 cp berror_stats global_berror.l72y770.f77.tco

