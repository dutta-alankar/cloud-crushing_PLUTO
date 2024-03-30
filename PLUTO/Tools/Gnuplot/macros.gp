#
# File: macros.gp
#
# Gnuplot script to define some useful macros for:
#
# - binary array format specifications
# - 1D slices
# - colormaps
# 
# Note: it requires the data type (dtype = "dbl" or dtype = "flt") to be set.
#
# Last Modified: 19 Aug, 2020 by A. Mignone (mignone@to.infn.it)
# 

# ----------------------------------------------------------------------
# Define the BINARR macro for plotting 2D single or double 
# precision data files.
# These macros can be used together with splot to 
# provide the necessary size and grid information, e.g.,
#
# gnuplot> splot "data.0002.dbl"  @BINARR
# ----------------------------------------------------------------------
if (dtype eq 'dbl') {
  dsize   = 8
  dformat = "%lf"
  print "> Data size set to double precision (dtype = 'dbl')"
}
if (dtype eq 'flt') {
  dsize = 4
  dformat = "%f"
  print "> Data size set to single precision (dtype = 'flt')"
}
print "> Grid size:"
print "  - nx1 = ",nx1
print "  - nx2 = ",nx2
print "  - nx3 = ",nx3

if (nx2 == 1){
  print "> Setting macro @BINARR (1D)"
  str1 = sprintf("bin array=%d format='%s' ",nx1, dformat)
  str2 = sprintf("dx=dx1 ");
  str3 = sprintf("skip=(nx1*dsize*nvar)");
  BINARR = str1.str2.str3    # concatenate strings
}
if (nx2 > 1){
  print "> Setting macro @BINARR (2D)"
  str1 = sprintf("bin array=%dx%d format='%s' ",nx1,nx2, dformat)
  str2 = sprintf("dx=dx1 dy=dx2 origin=(x1beg, x2beg, 0) ");
  str3 = sprintf("skip=(nx1*nx2*nx3*dsize*nvar)");
  BINARR = str1.str2.str3    # concatenate strings
}

# ----------------------------------------------------------------------
# Define the XSLICE and YSLICE macros for plotting 1D profiles along,
# respectively, rows or columns of 2D binary data files.
#
# gnuplot> jcut = 2; plot "data.0002.dbl"  @XSLICE
# gnuplot> icut = 5; plot "data.0002.dbl"  @YSLICE
# ----------------------------------------------------------------------
Lx1  = x1end - x1beg
Lx2  = x2end - x2beg
Lx3  = x3end - x3beg

icut = 0;
jcut = 0;

print "> Setting macro @XSLICE, @YSLICE"

str1 = sprintf("bin array=%d format='%s' ",nx1*nx2,dformat)
str2 = sprintf("dx=dx1 origin=(-jcut*Lx1 + 0.5*dx1, 0.0) ");
str3 = sprintf("skip=(nx1*nx2*dsize*nvar) every 1:1:(nx1*jcut):0:(nx1*jcut+nx1-1):0");
XSLICE = str1.str2.str3

str1 = sprintf("bin array=%d format='%s' ",nx1*nx2,dformat)
str2 = sprintf("dx=dx2/nx1 "); # two consecutive points in y are spaced by nx zones in x
str3 = sprintf("skip=(nx1*nx2*dsize*nvar) every nx1::icut:0:(icut+(nx2-1)*nx1):0");
YSLICE = str1.str2.str3      # concatenate strings

# ----------------------------------------------------------------------
# Define a few colormap macros
# ----------------------------------------------------------------------
HOT = "rgbformulae 22,13,-31"
RED = "rgbformulae 21,22,23"
RYG = 'model RGB defined ( 0 "red", 0.5 "yellow", 1 "green" )'
RAINBOW = "rgbformulae 33,13,10 " # rainbow (blue-green-yellow-red)
JET = "defined (0  0.0 0.0 0.5, \
                1  0.0 0.0 1.0, \
                2  0.0 0.5 1.0, \
                3  0.0 1.0 1.0, \
                4  0.5 1.0 0.5, \
                5  1.0 1.0 0.0, \
                6  1.0 0.5 0.0, \
                7  1.0 0.0 0.0, \
                8  0.5 0.0 0.0 )"

GREEN_CT = 'model RGB defined (0 "black", 1 "slateblue1", 2 "white")'
print "  - available colormap macros: HOT, RED, RYG, RAINBOW, JET, GREEN_CT"

