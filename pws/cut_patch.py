import healpy as hp
import numpy as np
from matplotlib import pyplot
import os, sys
import pyfits
import datetime
from drizzlib import healpix2wcs

def rotate_coords(nside, theta, phi):
  r=hp.Rotator(coord=['C', 'G'])
  trot, prot=r(theta, phi)
  pix=hp.ang2pix(nside, trot, prot)
  return np.arange(len(theta))[pix]

freqs=[100, 143, 217, 353, 545, 857]
nfreq=len(freqs)
# Units in all the maps seem to be K_CMB, converted to MJy/sr for patches
units=['MJy/sr']*len(freqs)
conv_fact=[244.1, 371.74, 483.690, 287.450, 58.04, 2.27] # from http://adsabs.harvard.edu/abs/2014A%26A...571A...9P Table 6

# Output directory for files
outdir=sys.argv[1]
# Patch coordinates in decimal degrees
patch_RA=sys.argv[2]
patch_dec=sys.argv[3]

patch_num=0
nside=2048
npix=512
npix_mask=128
npatch=1

patch_cen_j=[patch_RA, patch_dec]
cen=[np.pi/2.-patch_dec*np.pi/180., patch_RA*np.pi/180.] # (theta,phi) in radians

r=hp.Rotator(coord=['C', 'G'])
cen_g=r(cen)
patch_cen_g=[cen_g[1]*180./np.pi, 90-cen_g[0]*180./np.pi]

bad_pix=-1.6375e+30

if not os.path.exists(outdir): os.mkdir(outdir)

if not os.path.exists(outdir+'/hp_proj'): os.mkdir(outdir+'/hp_proj')

fits_temp=pyfits.open('template.fits')
try:
  fits_temp[1].header.update('naxis1', npix)
  fits_temp[1].header.update('naxis2', npix)
  fits_temp[1].header.update('crval1', patch_cen_j[0])
  fits_temp[1].header.update('crpix1', npix/2+1)
  fits_temp[1].header.update('crval2', patch_cen_j[1])
  fits_temp[1].header.update('crpix2', npix/2+1)
except ValueError:
  fits_temp[1].header.update([('naxis1', npix), ('naxis2', npix), ('crval1', patch_cen_j[0]), ('crpix1', npix/2+1), ('crval2', patch_cen_j[1]), ('crpix2', npix/2+1)])

for i, freq in enumerate(freqs):
  hp_mapname='Maps/PlnkDX11d_%d.fits' % (freq)
  total_map=hp.read_map(hp_mapname)*conv_fact[i]
  if i==0:
    if hp.get_nside(total_map)!=nside:
      print 'Warning: specified nside (%d) not consistent with nside from healpix map (%d)' % (nside, hp.get_nside(total_map))
      nside=hp.get_nside(total_map)
    pix=hp.nside2resol(nside)*180./np.pi
    try:
      fits_temp[1].header.update('cdelt1', pix)
      fits_temp[1].header.update('cdelt2', pix)
    except ValueError:
      fits_temp[1].header.update([('cdelt1', pix), ('cdelt2', pix)])
    theta, phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    pix_J=rotate_coords(nside, theta, phi)
    theta_array=hp.gnomview(theta, rot=patch_cen_j, coord='C', xsize=npix, reso=pix*60., return_projected_map=True)
    phi_array=hp.gnomview(phi, rot=patch_cen_j, coord='C', xsize=npix, reso=pix*60., return_projected_map=True)
    pyplot.close('all')
    patchprops=open(outdir+'/PatchGeomProps.csv', 'w')
    patchprops.write(str(nside)+', '+str(npix)+', '+str(npix_mask)+', 0, '+str(npatch)+', -1, 0.0001\n')
    patchprops.write('    '+str(patch_num)+',-1,   -1,   -1,      0,'+str(np.pi/2.-patch_cen_j[1]*np.pi/180)+','+str(patch_cen_j[0]*np.pi/180)+','+str(theta_array[0,-1])+','+str(phi_array[0,-1])+','+str(theta_array[0,0])+','+str(phi_array[0,0])+','+str(theta_array[-1,-1])+','+str(phi_array[-1,-1])+','+str(theta_array[-1,0])+','+str(phi_array[-1,0])+',      0,      0, 0,  0.003,      7,      7,     -1,0.000429,     -1, 0,     -1,     -1,     -1,     -1,  0.003,      7,'+str(np.pi/2.-patch_cen_j[1]*np.pi/180)+','+str(patch_cen_j[0]*np.pi/180)+', -2e+10, -2e+10,'+str(bad_pix)+','+str(bad_pix)+','+str(bad_pix)+'\n')
    patchprops.close()
    # Construct pixel mask
    cenpix=npix/2
    theta_mask=theta_array[cenpix-npix_mask/2:cenpix+npix_mask/2,cenpix-npix_mask/2:cenpix+npix_mask/2]
    theta_mask=theta_mask.filled(fill_value=theta_mask.fill_value)
    phi_mask=phi_array[cenpix-npix_mask/2:cenpix+npix_mask/2,cenpix-npix_mask/2:cenpix+npix_mask/2]
    phi_mask=phi_mask.filled(fill_value=phi_mask.fill_value)
    corners=np.array([[theta_mask[0,0],phi_mask[0,0]], [theta_mask[-1,0],phi_mask[-1,0]], [theta_mask[-1,-1],phi_mask[-1,-1]], [theta_mask[0,-1], phi_mask[0,-1]]])
    vec_corners=hp.ang2vec(corners[:,0], corners[:,1])
    pix_mask=hp.query_polygon(nside, vec_corners)
    mask=np.zeros_like(total_map)
    mask[pix_mask]=1.
    hp.write_map(outdir+'/PwS_RejectionMask.fits', mask, coord='C')
  #hp.gnomview(total_map, rot=patch_cen_g, xsize=npix, reso=pix*60.)
  freqmap=hp.gnomview(total_map, rot=patch_cen_j, coord=['G','C'], xsize=npix, reso=pix*60., return_projected_map=True)
  pyplot.close('all')
  freqmap=freqmap.filled(fill_value=freqmap.fill_value)
  freqmap=freqmap[:,::-1]
  now=datetime.datetime.now()
  try:
    fits_temp[0].header.update('date', now.strftime('%Y-%m-%dT%H:%M:%S'))
    fits_temp[1].header.update('freq', freq*1e9)
    fits_temp[1].header.update('bunit', units[i])
  except ValueError:
    fits_temp[0].header.update([('date', now.strftime('%Y-%m-%dT%H:%M:%S'))])
    fits_temp[1].header.update([('freq', freq*1e9), ('bunit', units[i])])
  fits_temp[1].data=freqmap
  fits_temp.verify('fix')
  mapout='%s/hp_proj/map%d_%04d_%04d.fits' % (outdir,nside,freq,patch_num)
  fits_temp.writeto(mapout, clobber=True)
  mapout2='%s/map%d_%04d_%04d.fits' % (outdir,nside,freq,patch_num)
  freqmap2=healpix2wcs.healpix2wcs(hp_mapname, header=mapout, header_hdu=1, output=mapout2, clobber=True)
  # PwS needs two hdus so copy from the template
  f=pyfits.open(mapout)
  f[1].data=freqmap2*conv_fact[i]
  f.writeto(mapout2, clobber=True)

