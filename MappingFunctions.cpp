//
// Created by pier on 9/21/16.
//


#include <cmath>
#include "MappingFunctions.hpp"
#include <Eigen>

using namespace Eigen;

/*!
 * Converts geodetic coordinate to geocentric ECEF coordinates
 % ECEF  Earth-Centre-Earth-Fixed
 *
 * Input(s) :
 * float lat            :- the latitude to be converted (rad)
 * float lon            :- the longitutde to be converted (rad)
 * float alt            :- the altitude; the height of platform (m)
 *
 * Return:
 * double3               :- the converted ECEF coordinates in x, y, z.
 *
 */
double3 geodetic2geocentric(double lat, double lon, double alt)
{
    double3 xyz;

    double a = 6378137.0;
    double e = 8.1819190842622e-2;
    double N = a / sqrt(1.0 - e * e * sin(lat) * sin(lat));

    xyz.x = (N + alt) * cos(lat) * cos(lon);
    xyz.y = (N + alt) * cos(lat) * sin(lon);
    xyz.z = ((1.0 - e * e) * N + alt) * sin(lat);

    return xyz;
}

/*!
 * Calculate the rotation matrix to rotate local vertical axes to
 * East-North-Up orientation
 *
 * @param phi0 Corresponds to lat, in radians
 * @param lambda0 Corresponds to lon, in radians
 * @return matrix multiplied rotation matrix
 */
Matrix3d geocentric2lvRotationMatrix(double phi0, double lambda0)
{
    Matrix3d m1, m2;
    m1 << 1.0, 0, 0,
            0, sin(phi0), cos(phi0),
            0, -cos(phi0), sin(phi0);

    m2 << -sin(lambda0), cos(lambda0), 0,
            -cos(lambda0), -sin(lambda0), 0,
            0, 0, 1.0;

    Matrix3d m3 = m1 * m2;
    return m1 * m2;
}


/*!
 * Converts (X,Y,Z) in geocentric (ECEF) to local vertical coordinates
 * East-North-Up (ENU)wrt to the origin defined by
 * (lat,long,height) = (phi0, lambda0, h0)
 *
 * @param xp
 * @param yp
 * @param zp
 * @param ref_lat
 * @param ref_lon
 * @param height
 * @return
 */
double3 geocentric2lv(double xp, double yp, double zp, double ref_lat, double ref_lon, double height)
{
    double3 origin = geodetic2geocentric(ref_lat, ref_lon, height);
    Matrix3d rotationMatrix = geocentric2lvRotationMatrix(ref_lat, ref_lon);

    Vector3d P{xp - origin.x, yp - origin.y, zp - origin.z};
    Vector3d ans = rotationMatrix * P;
    double3 returnAns{ans[0], ans[1], ans[2]};
    return returnAns;
}
