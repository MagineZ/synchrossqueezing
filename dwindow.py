import numpy as np

def dwindow(h):
    """
    %DWINDOW Derive a window.
    %	DH=DWINDOW(H) derives a window H.
    %
    %	Example : 
    %	 plot(dwindow(tftb_window(210,'hanning')))
    %
    %	See also WINDOW.
    
    %	F. Auger, August 1994, July 1995.
    %	Copyright (c) 1996 by CNRS (France).
    %
    %  This program is free software; you can redistribute it and/or modify
    %  it under the terms of the GNU General Public License as published by
    %  the Free Software Foundation; either version 2 of the License, or
    %  (at your option) any later version.
    %
    %  This program is distributed in the hope that it will be useful,
    %  but WITHOUT ANY WARRANTY; without even the implied warranty of
    %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %  GNU General Public License for more details.
    %
    %  You should have received a copy of the GNU General Public License
    %  along with this program; if not, write to the Free Software
    %  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
    """

    [hrow,hcol]=h.shape 
    if not (hrow==1):
        raise ValueError('h must have only one row')
    
    Lh = int((hcol-1)/2)
    step_height = (h[0,0]+h[0,-1])/2
    ramp = (h[0,-1]-h[0,0])/(hcol-1)
    h2 = np.append([0],h-step_height-ramp*np.arange(-Lh,Lh+1))
    h2 = np.append(h2,[0])
    Dh = (h2[2:hcol+2]-h2[0:hcol])/2+ramp
    Dh[0]+=step_height 
    Dh[-1]-=step_height
    return np.array([Dh])

