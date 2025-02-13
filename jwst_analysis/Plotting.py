from astropy.coordinates import SkyCoord
import astropy.units as u

def get_FOV(cube, pos=(69.9744522, 26.0526275), width=0.00119):
    """Get the pixel bounds to match the provided FOV

    Args:
        cube (Cube): cube to get the FOV for
        pos (tuple): ra and dec of the center of the FOV in deg
        width (float): size of the FOV in degrees

    Returns:
        left, right, down up: pixel bounds of given FOV
    """
    # get center of FOV as a sky coordinate
    pos_sky = SkyCoord(pos[0]*u.um, pos[1]*u.um, frame='ircs')
    # convert to a pixel position
    pos_pix = pos_sky.to_pixel(cube.wcs.celestial)
    # convert the size of FOV to pixels
    deg_per_pixel = cube.header['CDELT1']
    width_pix = width / deg_per_pixel
    # get pixel bounds of FOV
    left = pos_pix[0] - (0.5* width_pix)
    right = pos_pix[0] + (0.5* width_pix)
    down = pos_pix[1] - (0.5* width_pix)
    up = pos_pix[1] + (0.5* width_pix)

    return left, right, down, up

def get_scale_bar_size(cube, size_au, dpc):
    """Get the pixel size given a size in au to make a scale bar

    Args:
        cube (Cube): cube to get size of scale bar for
        size_au (float): size of the scale bar in au
        dpc (float): distance to the source in pc

    Returns:
        size_pix (float): size of the scale bar in pixels
    """
    # get pixel size in degrees
    deg_per_pixel = cube.header['CDELT1']
    # convert to arcsec per pixel
    arcsec_per_pixel = (deg_per_pixel*u.deg).to(u.arcsec)
    # use distance to get au per pixel
    au_per_pix = arcsec_per_pixel.value * dpc
    # get size of the scale bar in pixels
    size_pix = size_au / au_per_pix
    return size_pix
