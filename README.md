## Fortran codes for vertical deflection computation on or outside geoid using Vening-Meinesz numerical integral
https://www.zcyphygeodesy.com/en/h-nd-138.html
## [Algorithm purpose]
    Using the generalized Vening-Meinesz numerical integral, compute the vertical deflection (″, SW, to south, to west) on or outside the geoid from the ellipsoidal height grid of the equipotential surface and its gravity anomaly or disturbance (mGal) grid.
    The generalized Vening-Meinesz formula is derived from the generalized Stokes/Hotine formula and belongs to the solution of the Stokes boundary value problem. Which requires the integrand gravity anomaly/disturbance to be on the equipotential surface.
    It is usually necessary to employ the remove-restore scheme with a reference geopotential model to use the finite radius for gravity field integral. Firstly, remove model gravity anomaly/disturbance on the boundary surface, then integrate to obtain the residual height anomaly at the calculation point, and finally restore the model height anomaly at the calculation point.
    The equipotential surface can be constructed from a global geopotential model (not greater than 360 degrees), which can also be represent by a normal (orthometric) equiheight surface with the altitude of not more than ten kilometers.
![](https://24192633.s21i.faiusr.com/2/ABUIABACGAAg8OzltwYo-OvcowIwpQ047gg.jpg)
## [Main program for test entrance]
    Vening_Meinesznumintegral.f90 The record format of the input calculation point file: ID (point no / point name), longitude (decimal degrees), latitude (decimal degrees), ellipsoidal height (m)......
    Input the ellipsoidal height grid file of the equipotential boundary surface, which employed to calculate the integral distance.
    Input the residual gravity anomaly/disturbance (mGal) grid file on the equipotential surface.
    Input parameter mode - when mode =0 from gravity anomaly, and when mode = 1 from gravity disturbance.
    The record format of the output file reslt.txt: Behind the record of the calculation point file, appends two columns of residual vertical deflection southward and westward calculated.
## (1) Algorithm module for Vening-Meinesz numerical integral from gravity anomaly
    VManomalyBLH(BLH,gra,sfh,nlat,nlon,hd,dr,vm,GRS)
    Input parameters: BLH(3) - longitude (decimal degrees), latitude (decimal degrees), ellipsoidal height (m) of the calculation point.
    Input parameters: sfh(nlat,nlon) - the ellipsoidal height grid of the equipotential boundary surface, which employed to calculate the integral distance.
    Input parameters: gra(nlat,nlon) - the residual gravity anomaly (mGal) grid on the equipotential surface.
    Input parameters: dr, hd(6) - the integral radius (m) and grid specification parameters (minimum and maximum longitude, minimum and maximum latitude, longitude and latitude intervals of a cell grid).
    Input parameters: GRS(6) - gm, ae, j2, omega, 1/f, default value
    Return parameters: vm(2) - the calculated residual vertical deflection (″, SW, to south, to west).
## (2) Algorithm module for Vening-Meinesz numerical integral from gravity disturbance
    VMdisturbBLH(BLH,gra,sfh,nlat,nlon,hd,dr,vm,GRS)
    Input parameters: rga(nlat,nlon) - the residual gravity disturbance (mGal) grid on the equipotential surface.
    Return parameters: vm(2) - the calculated residual vertical deflection (″, SW, to south, to west).
## (3) Calculation module for the normal gravity field
    normdjn(GRS,djn); GNormalfd(BLH,NFD,GRS)
    Return parameters: NFD(5) - the normal geopotential (m2/s2), normal gravity (mGal), normal gravity gradient (E), normal gravity line direction (', expressed by its north declination relative to the center of the Earth center of mass) or normal gravity gradient direction (', expressed by its north declination relative to the Earth center of mass).
## (4) Calculation module for Legendre functions and their derivatives to ψ
    LegPn_dt2(pn,dp1,dp2,n,t) ! t=cos ψ
## (5) Algorithm library for transforming of geodetic coordinates
    BLH_RLAT(GRS, BLH, RLAT); BLH_XYZ(GRS, BLH, XYZ)
    RLAT_BLH(GRS, RLAT, BLH)
## (6) Algorithm library for interpolation point value from numerical grid
    CGrdPntD(lon,lat,dt,row,col,hd); CGrdPntD2(lon,lat,dt,row,col,hd)
    CShepard(lon,lat,dt,row,col,hd); Gauss2D(lon,lat,dt,row,col,hd)
## (7) Other auxiliary modules
    PickRecord(str0, kln, rec, nn)
## [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler for any operating system. No external link library required.
## [Algorithmic formula] PAGravf4.5 User Reference https://www.zcyphygeodesy.com/en/
    1.4.1 Format convention for geodetic data file
    7.9.2 Vening-Meinesz integral formulas outside geoid
    7.1(4) Low-dgree Legendre function and its first and second derivative algorithms
The zip compression package includes the test project in visual studio 2017 - intel fortran integrated environment, DOS executable test file and all input and output data.
![](https://24192633.s21i.faiusr.com/2/ABUIABACGAAg8OzltwYo1-HBwgYwpQ047gg.jpg)
