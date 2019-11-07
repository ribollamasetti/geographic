const math = require( 'mathjs' );
const polyval = require( 'compute-polynomial' );
function geographic () {
{ // Initialize and construct setters and getters
    var me = this;

    /**
     * Initialize properties
     */
    var init = ()=>{
            ellipsoid = { radius:0, flattening:0, eccentricity:0 };
            me.setDefaultEllipsoid();
            origin = { lat:0, lon:0, h:0 };
    }

    var ellipsoid;
    var origin;

    /**
     * Set the ellipsoid to WGS84
     */
    this.setDefaultEllipsoid = () => {
            ellipsoid.radius = 6378137;
            ellipsoid.flattening = 1 / 298.257223563;
            ellipsoid.eccentricity = flat2ecc(1 / 298.257223563);
    }

    /**
     * Set the ellipsoid to the one with given earth radius and flattening
     */
    this.setEllipsoid = (radius, flattening) => {
            ellipsoid.radius = radius;
            ellipsoid.flattening = flattening;
            ellipsoid.eccentricity = flat2ecc(flattening);
    }

    /**
     * Get the ellipsoid object
     */
    this.getEllipsoid = () => { return ellipsoid; }

    /**
     * Set the projection origin object
     */
    this.setOrigin = (lat0, lon0, h0) => {
            origin.lat = lat0;
            origin.lon = lon0;
            origin.h = h0;
    }

    /**
     * Get the projection origin object
     */
    this.getOrigin = () => { return origin; }
}
{ // Conversion functions
        /**
     * Inverse transverse Mercator projection.
     * Convert (x[m],y[m]) coordinates from transverse Mercator projection to geografic position (lat[°],lon[°])
     */
    this.tranmerc_inv = ( x, y )=>{
        var S = [1,1];
        var Z = 0;
        var a = ellipsoid.radius;
        var e2 = math.re(ellipsoid.eccentricity**2); // verificar se e2 não pode ser um parametro global
        var f = ellipsoid.flattening; // verificar de trocar os usos de f para ellipsoid.flattening
        var e2m = 1 - e2;
        var cc = Math.sqrt(e2m) * Math.exp(eatanhe(1,e2));
        var n = f / ( 2 - f );
        var bet = betf( n );
        var bl = (1-f)*(A1m1f(n)+1);
        var a1 = bl * a;
        if ( origin.lat == 0 ) {
            var y0 = 0;
        } else {
            var sbet0 = NaN;
            var cbet0 = NaN;
            [sbet0, cbet0] = sincosdx( LatFix( origin.lat ) );
            [sbet0, cbet0] = norm2((1-f) * sbet0, cbet0);
            var y0 = a1 * (Math.atan2(sbet0, cbet0) + SinCosSeries(true, sbet0, cbet0, C1f(n)));
        }
        y = y + y0 + Z;
        var xi = y / a1;
        var eta = x / a1 + Z;
        var xisign = 1 - 2 * (xi<0);
        var etasign = 1 - 2 * (eta<0);
        xi = xi * xisign;
        eta = eta * etasign;
        var backside = xi>(Math.PI/2);
        if (backside) {
            xi = Math.PI - xi;
        }
        var c0 = Math.cos(2*xi);
        var ch0 = Math.cosh(2*eta);
        var s0 = Math.sin(2*xi);
        var sh0 = Math.sinh(2*eta);
        a = math.multiply(2,math.complex(c0*ch0,-s0*sh0));
        var j = 6;
        y0 = math.complex(Z);
        var y1 = y0;
        var z0 = y0;
        var z1 = y0;
        if (j%2) {
            y0 = y0 + bet[j];
        }
        for (let index = j; index >= 1; index = index - 2) {
            y1 = math.subtract(math.subtract(math.multiply(a,y0),y1),bet[index-1]);
            z1 = math.subtract(math.subtract(math.multiply(a,z0),z1),2 * index * bet[index-1]);
            y0 = math.subtract(math.subtract(math.multiply(a,y1),y0),bet[index-2]);
            z0 = math.subtract(math.subtract(math.multiply(a,z1),z0),2 * (index-1) * bet[index-2]);
        }
        a = math.divide(a,2);
        z1 = math.sum(math.subtract(1,z1),math.multiply(z0,a));
        a = math.complex(s0*ch0,c0*sh0);
        y1 = math.sum(math.complex(xi,eta),math.multiply(y0,a));
        var gam = math.atan2(z1.im,z1.re) * 180 / Math.PI;
        var k = bl / math.abs(z1);
        var xip = y1.re;
        var etap = y1.im;
        var s = math.sinh(etap);
        var c = math.max([0,math.cos(xip)]);
        var r = math.hypot(s,c);
        var lon = math.atan2(s,c) * 180 / Math.PI;
        var sxip = math.sin(xip);
        var tau = tauf(sxip/r,e2);
        var lat = math.atan2(tau,1+Z) * 180 / Math.PI;
        gam = gam + math.atan2(sxip*math.tanh(etap),c) * 180 / Math.PI;
        c = r != 0;
        if (c) {
            k = k * math.sqrt(e2m+e2/(1+tau**2))*math.hypot(1,tau)*r;
        }
        c = !c;
        if (c) {
            lat = 90;
            lon = 0;
            k = k*cc;
        }
        lat = lat * xisign;
        if (backside) {
            lon = 180 - lon;
        }
        lon = lon * etasign;
        lon = AngNormalize(lon+AngNormalize(origin.lon));
        if (backside) {
            gam = 180 - gam;
        }
        gam = AngNormalize(gam*xisign*etasign);
        return {lat,lon};
    }

    /**
     * Foreward transverse Mercator projection.
     * Convert (lat[°],lon[°]) coordinates from geographic position to transverse Mercator projection (x[m],y[m])
     */
    this.tranmerc_fwd = ( lat, lon )=>{
        var S = [1,1];
        var Z = 0;
        var maxpow = 6;
        var a = ellipsoid.radius;
        var e2 = math.re(ellipsoid.eccentricity**2); // verificar se e2 não pode ser um parametro global
        var f = ellipsoid.flattening; // verificar de trocar os usos de f para ellipsoid.flattening
        var e2m = 1 - e2;
        var cc = Math.sqrt(e2m) * Math.exp(eatanhe(1,e2));
        var n = f / ( 2 - f );
        var alp = alpf(n);
        var b1 = (1 - f) * (A1m1f(n) + 1);
        var a1 = b1 * a;
        lat = LatFix(lat) + Z;
        [lonlon,trash] = AngDiff(origin.lon, lon);
        lon = lonlon + Z;
        var latsign = 1 - 2 * (lat < 0);
        var lonsign = 1 - 2 * (lon < 0);
        lon = lon * lonsign;
        lat = lat * latsign;
        var backside = lon > 90;
        if (backside & lat==0) {
            latsign = -1;
        }
        if (backside) {
            lon = 180 - lon;
        }
        var sphi = null; var cphi = null; var slam = null; var clam = null;
        [sphi, cphi] = sincosdx(lat);
        [slam, clam] = sincosdx(lon);
        var tau = sphi / Math.max(Math.sqrt(2.22507385850720138309e-308), cphi);
        var taup = taupf(tau, e2);
        var xip = Math.atan2(taup, clam);
        var etap = Math.asinh(slam / Math.hypot(taup, clam));
        var gam = Math.atan2(slam * taup, clam * Math.hypot(1, taup)) * 180 / Math.PI;
        var k = Math.sqrt(e2m + e2 * cphi**2) * Math.hypot(1, tau) / Math.hypot(taup, clam);
        var c = !(lat != 90);
        if (c) {
            xip = pi/2;
            etap = 0;
            gam = lon;
            k = cc;
        }
        var c0 = Math.cos(2 * xip); var ch0 = Math.cosh(2 * etap);
        var s0 = Math.sin(2 * xip); var sh0 = Math.sinh(2 * etap);
        a = math.multiply(2,math.complex(c0 * ch0, -s0 * sh0));
        var j = maxpow;
        var y0 = math.complex(Z);
        var y1 = y0; var z0 = y0; var z1 = y0;
        if (j%2) {
            y0 = math.sum(y0, alp);
            z0 = math.sum(z0, 2*j * alp);
            j = j - 1;
        }
        for (let index = j; index > 0; index = index-2) {
            y1 = math.sum(math.subtract(math.multiply(a,y0),y1),alp[index-1]);
            z1 = math.sum(math.subtract(math.multiply(a,z0),z1),2*index*alp[index-1]);
            y0 = math.sum(math.subtract(math.multiply(a,y1),y0),alp[index-2]);
            z0 = math.sum(math.subtract(math.multiply(a,z1),z0),2*(index-1)*alp[index-2]);
        }
        a = math.divide(a,2);
        z1 = math.sum(math.subtract(1,z1),math.multiply(z0,a));
        a = math.complex(s0 * ch0, c0 * sh0);
        y1 = math.sum(math.complex(xip, etap),math.multiply(y0,a));
        gam = gam - Math.atan2(z1.im, z1.re)*180/Math.PI;
        k = k * (b1 * Math.abs(z1));
        var xi = y1.re; var eta = y1.im;
        if (backside) {
            xi = Math.PI - xi;
        }
        var y = a1 * xi * latsign;
        var x = a1 * eta * lonsign;
        if (backside) {
            gam = 180 - gam;
        }
        gam = AngNormalize(gam * latsign * lonsign);
        if (origin.lat==0) {
            y0 = 0;
        } else {
            var sbet0; var cbet0;
            [sbet0, cbet0] = sincosdx(LatFix(origin.lat));
            [sbet0, cbet0] = norm2((1-f) * sbet0, cbet0);
            y0 = a1 * (Math.atan2(sbet0, cbet0) + SinCosSeries(true, sbet0, cbet0, C1f(n)));
        }
        y = y - y0;
        return {x,y};
    }

    /**
     * convert to the UTM/UPS system to geographical coordinates, (lat,lon).
     * The input is (x,y) = (easting,northing), the zone which is either the UTM zone or 0 for UPS , and a hemisphere selector, isnorth (0 for the southern hemisphere, 1 for the northern)
     */
    this.utmups_inv = (x,y,zone,isnorth)=>{
        var Z = 0;
        x = x + Z;
        y = y + Z;
        zone = Math.floor(zone) + Z;
        isnorth = Boolean(isnorth + Z);
        Z = NaN;
        var geopos = {lat:NaN,lon:NaN};
        var gam = Z; var k = Z;
        var utm = zone > 0 & zone <= 60;
        if (utm) {
            geopos = utm_inv(zone, isnorth, x, y);
        }
        var ups = zone == 0;
        if (ups) {
            geopos = ups_inv(isnorth, x, y);
        }
        return geopos
    }

    /**
     * convert from geographical coordinates, (lat,lon), to the UTM/UPS system. The output is (x,y) = (easting,northing), zone which is either the UTM zone or 0 for UPS, and a hemisphere selector, isnorth (0 for the southern hemisphere, 1 for the northern). If setzone = -1 (the default), the standard choice is made between UTM and UPS and, if UTM, the standard zone is picked (the Norway and Svalbard exceptions are honored)
     */
    this.utmups_fwd = (lat,lon)=>{
        var Z = 0;
        var setzone = -1;
        lat = lat + Z;
        lon = lon + Z;
        var isnorth = lat >= 0;
        var zone = StandardZone(lat, lon, setzone);
        Z = NaN;
        var geopos = {x:NaN,y:NaN};
        var gam = Z;
        var k = Z;
        var utm = zone > 0;
        if (utm) {
            geopos = utm_fwd(zone, isnorth, lat, lon);
        }
        var ups = zone == 0;
        if (ups) {
            geopos = ups_fwd(isnorth, lat, lon);
        }
        if (Number.isNaN(geopos.x)) {
            zone = -4;
            isnorth = false;
        }
        geopos.zone = zone;
        geopos.isnorth = isnorth;
        return geopos;
    }

    /**
     * converts from local cartesian coordinates, x, y, z, centered at origin lat, lon, h to geodetic coordinates, lat, lon, h.
     * Latitudes and longitudes are in degrees;
     * h, h0, x, y, z are in meters.
     * x, y, z must be scalars.
     * lat0, lon0, h0 must be scalars.
     */
    this.loccart_inv = (x, y, z)=>{
        var S = 0;
        var num = 0;
        var Z = 0;
        x = x + Z;
        y = y + Z;
        z = z + Z;
        var X0 = NaN; var Y0 = NaN; var Z0 = NaN; var M0 = NaN;
        var geocent = geocent_fwd(origin.lat, origin.lon, origin.h);
        X0 = geocent.X; Y0 = geocent.Y; Z0 = geocent.Z; M0 = geocent.M;
        var r = math.multiply([x, y, z],math.transpose(M0));
        var X = math.add( math.subset( r, math.index( 0 ) ), X0 );
        var Y = math.add( math.subset( r, math.index( 1 ) ), Y0 );
        var Z = math.add( math.subset( r, math.index( 2 ) ), Z0 );
        var lat = NaN; var lon = NaN; var h = NaN; var M = NaN;
        [lat , lon , h , M] = geocent_inv(X, Y, Z);
        return {lat,lon,h};
    }

    this.loccart_fwd = (lat,lon,h)=>{
        var S = [1,1];
        var num = 1;
        var Z = 0;
        lat = lat + Z;
        lon = lon + Z;
        h = h + Z;
        var xyz = geocent_fwd(origin.lat, origin.lon, origin.h);
        var X0 = xyz.X; var Y0 = xyz.Y; var Z0 = xyz.Z; var M0 = xyz.M;
        xyz = geocent_fwd(lat , lon , h);
        var X = xyz.X; var Y = xyz.Y; var Z = xyz.Z; var M = xyz.M;
        r = math.multiply([X-X0, Y-Y0, Z-Z0],M0);
        x = r[0]; y = r[1]; z = r[2];
        return {x,y,z}
    }
}
{ // Private acessory functions

    /**
     * Conversion from geographic to geocentric coordinates
     * @param {*} lat 
     * @param {*} lon 
     * @param {*} h 
     */
    function geocent_fwd(lat, lon, h){
        var z = 0;
        lat = LatFix(lat) + z;
        lon = lon + z;
        h = h + z;
        var a = ellipsoid.radius;
        var e2 = math.re(ellipsoid.eccentricity**2);
        var e2m = 1 - e2;
        var slam = NaN; var clam = NaN;
        var sphi = NaN; var cphi = NaN;
        [slam, clam] = sincosdx(lon);
        [sphi, cphi] = sincosdx(lat);
        var n = a/Math.sqrt(1 - e2 * sphi**2);
        var Z = (e2m * n + h) * sphi;
        var X = (n + h) * cphi;
        var Y = X * slam;
        var X = X * clam;
        var M = GeoRotation(sphi, cphi, slam, clam);
        return {X,Y,Z,M};
    }

    function swap(x,y){ return [y,x]; }

    function cbrtx(x){
        var y = Math.abs(x)**(1/3);
        return x<0 ? -y : y;
    }

    function geocent_inv(X, Y, Z){
        var z = 0;
        X = X + z;
        Y = Y + z;
        Z = Z + z;
        var a = ellipsoid.radius;
        var e2 = math.re(ellipsoid.eccentricity**2);
        var e2m = 1 - e2;
        var e2a = Math.abs(e2);
        var e4a = e2**2;
        var maxrad = 2 * a / (2.22044604925031308085e-16);
        var R = math.hypot(X, Y);
        var slam = Y / R;
        var clam = X / R;
        if (R==0) {
            slam = 0;
            clam = 1;
        }
        var h = Math.hypot(R, Z);
        if (e4a==0) {
            // Treat the spherical case.  Dealing with underflow in the general case with e2 = 0 is difficult.  Origin maps to N pole same as with ellipsoid.
            var Z1 = Z;
            Z1 = h==0 ? 1 : Z1;
            [sphi, cphi] = norm2(Z1, R);
            h = h - a;
        } else {
            // Treat prolate spheroids by swapping R and Z here and by switching the arguments to phi = atan2(...) at the end.
            var p = (R / a)**2;
            var q = e2m * (Z / a)**2;
            var r = (p + q - e4a) / 6;
            if (e2<0) {
                [p, q] = swap(p, q);
            }
            // Avoid possible division by zero when r = 0 by multiplying equations for s and t by r^3 and r, resp.
            var S = e4a * p * q / 4; // S = r^3 * s
            var r2 = r**2;
            var r3 = r * r2;
            var disc = S * (2 * r3 + S);
            var u = r;
            var fl2 = disc >= 0;
            if (fl2) {
                var T3 = S+r3;
            } else {
                var T3 = 0;
            }
            // Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss of precision due to cancellation.  The result is unchanged because of the way the T is used in definition of u. T3 = (r * t)^3
            T3 = T3 + (1 - 2 * (T3 < 0)) * Math.sqrt(fl2 ? disc : 0);
            // N.B. cbrtx always returns the real root.  cbrtx(-8) = -2.
            T = cbrtx(T3);
            if (fl2) {
                u = u + T + cvmgt(r2 / T, 0, T != 0);
            }
            // T is complex, but the way u is defined the result is real.
            if (!fl2) {
                var ang = Math.atan2(Math.sqrt(-disc), -(S + r3));
                // There are three possible cube roots.  We choose the root which avoids cancellation (disc < 0 implies that r < 0).
                u = u + 2 * r * Math.cos(ang / 3);
            }
            // guaranteed positive
            v = Math.sqrt(u**2 + e4a * q);
            // Avoid loss of accuracy when u < 0.  Underflow doesn't occur in e4 * q / (v - u) because u ~ e^4 when q is small and u < 0. u+v, guaranteed positive
            var uv = u + v;
            fl2 = u < 0;
            if (fl2) {
                uv = e4a * q / (v - u);
            }
            // Need to guard against w going negative due to roundoff in uv - q.
            var w = Math.max(0, e2a * (uv - q) / (2 * v));
            var k = uv / (Math.sqrt(uv + w**2) + w);
            if (e2>=0) {
                var k1 = k; var k2 = k + e2;
            } else {
                var k1 = k - e2; var k2 = k;
            }
            var sphi; var cphi;
            [sphi, cphi] = norm2(Z / k1, R / k2);
            h = (1 - e2m / k1) * Math.hypot(k1 * R / k2, Z);
            // Deal with exceptional inputs
            var c = e4a * q == 0 & r <= 0;
            if (c) {
                // This leads to k = 0 (oblate, equatorial plane) and k + e^2 = 0 (prolate, rotation axis) and the generation of 0/0 in the general formulas for phi and h.  using the general formula and division by 0 in formula for h.  So handle this case by taking the limits:
                // f > 0: z -> 0, k      ->   e2 * sqrt(q)/sqrt(e4 - p)
                // f < 0: R -> 0, k + e2 -> - e2 * sqrt(q)/sqrt(e4 - p)
                var zz = e4a - p; var xx = p;
                if (e2<0) {
                    [zz, xx] = swap(zz, xx);
                }
                zz = Math.sqrt(zz / e2m);
                xx = Math.sqrt(xx);
                var H = Math.hypot(zz, xx);
                sphi = zz / H;
                cphi = xx / H;
                sphi = c & Z < 0 ? - sphi : sphi ;
                h = - a *  H / e2a;
                if (e2>=0) {
                    h = c ? e2m * h : h;
                }
            }
        }
        var far = h > maxrad;
        if (far) {
            // We really far away (> 12 million light years); treat the earth as a point and h, above, is an acceptable approximation to the height. This avoids overflow, e.g., in the computation of disc below.  It's possible that h has overflowed to inf; but that's OK. Treat the case X, Y finite, but R overflows to +inf by scaling by 2.
            R = Math.hypot(X/2, Y/2);
            slam = Y / R; 
            slam = (far & R == 0) ? 0 : slam;
            clam = X / R; 
            clam = (far & R == 0) ? 1 : clam;
            H = Math.hypot(Z/2, R);
            sphi = Z/2 / H;
            cphi = R / H;
        }
        lat = Math.atan2(sphi, cphi)*180/Math.PI;
        lon = Math.atan2(slam, clam)*180/Math.PI;
        M = GeoRotation(sphi, cphi, slam, clam);
        return [lat,lon,h,M];
    }

    function GeoRotation(sphi, cphi, slam, clam){
        var S = [1,1];
        var M = new Array(3).fill(0).map((x)=>{ return new Array(3).fill(0); });
        var sphi = sphi; var cphi = cphi;
        var slam = slam; var clam = clam;
        // Local X axis (east) in geocentric coords
        M[0][0] = -slam;         M[1][0] =  clam;         M[2][0] = 0;
        // Local Y axis (north) in geocentric coords
        M[0][1] = -clam * sphi;  M[1][1] = -slam * sphi;  M[2][1] = cphi;
        // Local Z axis (up) in geocentric coords
        M[0][2] =  clam * cphi;  M[1][2] =  slam * cphi;  M[2][2] = sphi;
        return M;
    }
    
    function StandardZone(lat, lon, setzone){
        var INVALID = -4;
        var UTM = -2;
        var MINZONE = 0;
        var MAXZONE = 60;
        var UPS = 0;
        var zone = Math.floor(setzone) + 0;
        if (!(zone>=INVALID & zone <= MAXZONE)) {
            zone = INVALID;
        }
        var g = zone < MINZONE & zone != INVALID;
        var c = Math.abs(lat) <= 90 & Number.isFinite(lon);
        if (g & !c) {
            zone = INVALID;
        }
        g = g & c;
        c = zone == UTM | (lat >= -80 & lat < 84);
        var u = g & c;
        var ilon = (Math.floor(lon) + 180)%360 - 180;
        var z = Math.floor((ilon + 186) / 6);
        // Norway exception
        var exception = z == 31 & Math.floor(lat / 8) == 7 & ilon >= 3;
        if (exception) { z = 32; }
        // Svalbard exception
        exception = lat >= 72 & ilon >= 0 & ilon < 42;
        if (exception) { z = 2 * Math.floor((ilon + 183)/12) + 1; }
        if (u) { zone = z; }
        if (g & !c) { zone = UPS; }
        return zone;
    }

    /**
     * Forward UTM projection
     * @param {*} zone 
     * @param {*} isnorth 
     * @param {*} lat 
     * @param {*} lon 
     */
    function utm_fwd(zone, isnorth, lat, lon){
        var lon0 = -183 + 6 * Math.floor(zone);
        var lat0 = 0;
        var bad = !(Math.abs((lon - lon0 + 180)%360 - 180) <= 60);
        var fe = 5e5;
        var fn = 100e5 * (1-isnorth);
        var k0 = 0.9996;
        var oldOrigin = {lat:origin.lat+0,lon:origin.lon+0,h:origin.h+0};
        origin.lat = lat0;
        origin.lon = lon0;
        geopos = me.tranmerc_fwd(lat, lon);
        x = geopos.x; y = geopos.y; gam = NaN; k = NaN;
        origin = oldOrigin;
        x = x * k0;
        y = y * k0;
        k = k * k0;
        bad = bad | !(Math.abs(x) <= 5e5 & y >= -91e5 & y <= 96e5);
        x = x + fe; y = y + fn;
        if (bad) {
            x = NaN; y = NaN; gam = NaN; k = NaN;
        }
        return {x,y};
    }

    /**
     * Inverse UTM projection
     * @param {*} zone 
     * @param {*} isnorth 
     * @param {*} x 
     * @param {*} y 
     */
    function utm_inv(zone, isnorth, x, y){
        var lon0 = -183 + 6 * Math.floor(zone);
        var lat0 = 0;
        var fe = 5e5;
        var fn = 100e5 * (1-isnorth);
        var k0 = 0.9996;
        x = x - fe;
        y = y - fn;
        var geopos = {lat:NaN,lon:NaN};
        var gam = NaN;
        var k = NaN;
        var bad = !(Math.abs(x) <= 5e5 & y >= -91e5 & y <= 96e5);
        if (!bad) {
            x = x / k0;
            y = y / k0;
            geopos = me.tranmerc_inv(x, y);
            k = k * k0;
        } 
        return geopos;
    }

    /**
     * Forward UPS projection
     * @param {*} isnorth 
     * @param {*} lat 
     * @param {*} lon 
     */
    function ups_fwd(isnorth, lat, lon){
        var fe = 20e5;
        var fn = 20e5;
        var k0 = 0.994;
        var x = NaN; var y = NaN; var gam = NaN; var k = NaN;
        [x, y] = polarst_fwd(isnorth, lat, lon);
        x = x * k0;
        y = y * k0;
        k = k * k0;
        var lim = (13 - 5 * isnorth) * 1e5;
        var bad = !(Math.abs(x) <= lim & Math.abs(y) <= lim);
        x = x + fe;
        y = y + fn;
        if (bad) {
            x = NaN; y = NaN; gam = NaN; k = NaN;
        }
        return {x,y};
    }

    /**
     * Inverse UPS projection
     * @param {*} isnorth 
     * @param {*} x 
     * @param {*} y 
     */
    function ups_inv(isnorth, x, y){
        var fe = 20e5;
        var fn = 20e5;
        var k0 = 0.994;
        x = x - fe;
        y = y - fn;
        var lim = (13 - 5 * isnorth) * 1e5;
        var geopos = {lat:NaN,lon:NaN};
        var gam = NaN;
        var k = NaN;
        var bad = !(Math.abs(x) <= lim & Math.abs(y) <= lim);
        if (!bad) {
            x = x / k0;
            y = y / k0;
            geopos = polarst_inv(isnorth, x, y);
            k = k * k0;
        }
        return geopos;
    }

    /**
     * Inverse polar stereographic projection
     * @param {*} isnorth 
     * @param {*} x 
     * @param {*} y 
     */
    function polarst_inv(isnorth, x, y){
        var Z = 0;
        var a = ellipsoid.radius;
        var e2 = ellipsoid.eccentricity**2;
        var e2m = 1-e2;
        var c = Math.sqrt(e2m) * Math.exp(eatanhe(1, e2));
        var isnorth = 2 * Boolean(isnorth) - 1;
        var rho = Math.hypot(x, y);
        var t = rho / (2 * a / c);
        var taup = (1 / t - t) / 2;
        var tau = tauf(taup, e2);
        var lat = Math.atan(tau)*180/Math.PI;
        if (rho==0) {
            lat = 90;
        }
        lat = isnorth * lat;
        var lon = Math.atan2(x, -isnorth * y)*180/Math.PI;
        return {lat,lon};
    }

    /**
     * Forward polar stereographic projection
     * @param {*} isnorth 
     * @param {*} lat 
     * @param {*} lon 
     */
    function polarst_fwd(isnorth, lat, lon){
        var Z = 0;
        var overflow = 1/eps^2;
        var a = ellipsoid.radius;
        var e2 = ellipsoid.eccentricity**2;
        var e2m = 1-e2;
        var c = Math.sqrt(e2m) * Math.exp(eatanhe(1, e2));
        isnorth = 2 * Boolean(isnorth) - 1;
        var lat = LatFix(lat) * isnorth;
        var tau = Math.tan(lat)*180/Math.PI;
        if (Math.abs(lat)==90) {
            tau = Math.sign(lat) * overflow;
        }
        var taup = taupf(tau, e2);
        var rho = Math.hypot(1, taup) + Math.abs(taup);
        if (taup>=0) {
            rho = lat!=90 ? 1/rho : 0; //cvmgt(1/rho,0,lat!=90);
        }
        rho = rho * (2 * a / c);
        [x, y] = sincosdx(lon);
        x = rho * x;
        y = -isnorth * rho * y;
        return [x,y];
    }

    /**
     * Conditional merge
     * @param {*} x 
     * @param {*} y 
     * @param {*} p 
     */
    function cvmgt(x,y,p){ return p ? x : y; }
    
    /**
     * Convert flattening to eccentricity
     * @param {*} f flattening
     */
    function flat2ecc ( f ) {
        return Math.sqrt( f * ( 2 - f ));
    }

    /**
     * Convevrt eccentricity to flattening
     * @param {*} e eccentricity
     */
    function ecc2flat ( e ) {
        return ( e**2 ) / ( 1 + Math.sqrt( 1 - ( e**2 )));
    }

    /**
     * returns e*atanh(e*x) where e = sqrt(e2)
     * @param {*} x 
     * @param {*} e2 
     */
    function eatanhe ( x, e2 ) {
        var e = Math.sqrt(Math.abs(e2));
        if ( e2 >= 0 ) {
            return e * Math.atanh( e * x );
        } else {
            return -e * Math.atan( e * x );
        }
    }
    
    function betf ( n ) {
        var betcoeff = [
            384796, -382725, -6720, 932400, -1612800, 1209600, 2419200, 
            -1118711, 1695744, -1174656, 258048, 80640, 3870720, 
            22276, -16929, -15984, 12852, 362880, 
            -830251, -158400, 197865, 7257600, 
            -435388, 453717, 15966720, 
            20648693, 638668800 ];
        var bet = new Array(6).fill(0);
        var o = 1;
        var d = n;
        for (let index = 1; index <= 6; index++) {
            var m = 6 - index;
            var valPoly = math.subset( betcoeff, math.index( math.range( o-1, o+m )));
            try { valPoly = polyval( valPoly, n ); } catch (error) { }
            bet[index-1] = d * valPoly / math.subset( betcoeff, math.index(o + m) );
            o = o + m + 2;
            d = d * n;
        }
        return bet;
    }

    function alpf(n){
        var alpcoeff = [ 31564, -66675, 34440, 47250, -100800, 75600, 151200, -1983433, 863232, 748608, -1161216, 524160, 1935360, 670412, 406647, -533952, 184464, 725760, 6601661, -7732800, 2230245, 7257600, -13675556, 3438171, 7983360, 212378941, 319334400, ];
        var alp = new Array(6).fill(0);
        var o = 1;
        var d = n;
        for (let index = 1; index <= 6; index++) {
            var m = 6 - index;
            var valPoly = math.subset( alpcoeff, math.index( math.range( o-1, o+m )));
            try { valPoly = polyval( valPoly, n ); } catch (error) { }
            alp[index-1] = d * valPoly / alpcoeff[o + m];
            o = o + m + 2;
            d = d * n;
        }
        return alp;
    }

    function A1m1f ( epsi ) {
        var coeff = [ 1, 4, 64, 0, 256 ];
        var eps2 = epsi ** 2;
        var t = polyval(math.subset( coeff, math.index(math.range(0,4))), eps2 ) / math.subset( coeff, math.index(4));
        return ( t + epsi ) / ( 1 - epsi );
    }

    /**
     * Compute sine and cosine with argument in degrees
     * @param {*} x angle in degrees
     */
    function sincosdx ( x ) {
        var radX = x * Math.PI / 180;
        return [ Math.sin( radX ), Math.cos( radX ) ];
    }

    /**
     * returns x is it is in the range [-90, 90]; otherwise it returns NaN.
     * @param {*} x latitude
     */
    function LatFix ( x ) {
        if (Math.abs(x) < 90) {
            return x;
        } else {
            return NaN;
        }
    }

    /**
     * normalize x and y so that x^2 + y^2 = 1
     * @param {*} x 
     * @param {*} y 
     */
    function norm2 ( x, y ) {
        var r = Math.hypot( x, y );
        return [ x/r, y/r ];
    }

    /**
     * Evaluate a sine or cosine series using Clenshaw summation
     * y = SINCOSSERIES(sinp, sinx, cosx, c) evaluate
     *  y = sum(c[i] * sin( 2*i    * x), i, 1, n), if  sinp
     *  y = sum(c[i] * cos((2*i-1) * x), i, 1, n), if ~sinp
     * where n is the size of c.  x is given via its sine and cosine in sinx
     * and cosx.  sinp is a scalar.  sinx, cosx, and y are K x 1 arrays.  c is
     * a K x N array.
     * @param {*} sinp 
     * @param {*} sinx 
     * @param {*} cosx 
     * @param {*} c 
     */
    function SinCosSeries ( sinp, sinx, cosx, c ) {
        var sizeVec = math.size(c)[0];
        var n = sizeVec;
        var ar = 2 * ( cosx - sinx ) * ( cosx + sinx );
        var y1 = 0;
        var len1 = math.size(c)[0];
        if (n%2) {
            var y0 = c.subset( math.index(math.range(0,len1)), n );
            n = n - 1;
        } else {
            var y0 = y1;
        }
        for (let index = n; index >= 1; index = index - 2) {
            y1 = ar * y0 - y1 + c[index-1];
            y0 = ar * y1 - y0 + c[index-2];
        }
        if (sinp) {
            return 2 * sinx * cosx * y0;
        } else {
            return cosx * (y0 - y1);
        }
    }

    /**
     * returns tangent of phi in terms of taup the tangent of chi.  e2, the square of the eccentricity, is a scalar; taup can be any shape.
     * @param {*} taup 
     * @param {*} e2 
     */
    function tauf ( taup, e2 ) {
        var numit = 5;
        var e2m = 1-e2;
        var tau = taup / e2m;
        var stol = 0.1 * math.sqrt(2.2204460493*10**(-16)) * math.max([1,math.abs(taup)]);
        var g = !math.isZero(1/tau);
        for (let index = 1; index <= numit; index++) {
            if (!g) {
                break;
            }
            var tau1 = math.hypot(1,tau);
            var sig = math.sinh(eatanhe(tau/tau1, e2));
            var taupa = math.hypot(1,sig)*tau-sig*tau1;
            var dtau = (taup - taupa) * (1 + e2m * tau**2) / (e2m * tau1 * math.hypot(1, taupa));
            tau = tau + dtau;
            g = g & (math.abs(dtau) >= stol);
        }
        return tau;
    }

    /**
     * Reduce angle to range (-180, 180]
     * @param {*} x 
     */
    function AngNormalize ( x ) {
        x %= 360;
        if (x>180) {
            x -= 360;
        }
        if (x<=-180) {
            x += 360;
        }
        return x;
    }

    function C1f ( epsi ) {
        var nC1 = 6;
        var coeff = [ 
            -1, 6, -16, 32, 
            -9, 64, -128, 2048, 
            9, -16, 768, 
            3, -5, 512, 
            -7, 1280, 
            -7, 2048, ];
        var C1 = new Array(nC1).fill(0);
        var eps2loc = epsi**2;
        var d = epsi;
        var o = 1;
        for (let index = 1; index <= nC1; index++) {
            var m = Math.floor((nC1 - index) / 2);
            var poly = math.subset( coeff, math.index(math.range(o-1, o + m)));
            try { poly = polyval( poly, eps2loc); } catch (error) { }
            C1[index-1] = d * poly / math.subset( coeff, math.index(o + m));
            o = o + m + 2;
            d = d * epsi;
        }
        return C1;
    }

    /**
     * Compute angle difference accurately
     * computes z = y - x, reduced to (-180,180].  d = round(z) and e = z - round(z).  x and y can be any compatible shapes
     * @param {*} x
     * @param {*} y
     */
    function AngDiff(x,y){
        var d = null; var t = null; var e = null;
        [d, t] = sumx(AngNormalize(-x), AngNormalize(y));
        d = AngNormalize(d);
        if (d==180 & t>0) {
            d = -180;
        }
        [d, e] = sumx(d, t);
        return [d,e];
    }

    /**
     * Error free sum
     * returns the rounded sum u + v in s and the error in t, such that s + t = u + v, exactly
     * @param {*} u 
     * @param {*} v 
     */
    function sumx(u,v){
        var s = u + v;
        var up = s - v;
        var vpp = s - up;
        var up = up - u;
        var vpp = vpp - v;
        var t = -(up + vpp);
        return [s,t];
    }

    /**
     * returns tangent of chi in terms of tau the tangent of phi.  e2, the square of the eccentricity, is a scalar.
     * @param {*} tau 
     * @param {*} e2 
     */
    function taupf(tau,e2){
        var tau1 = Math.hypot(1, tau);
        var sig = Math.sinh( eatanhe( tau / tau1, e2 ) );
        return (Math.hypot(1, sig) * tau - sig * tau1);
    }


}
    init();
}
module.exports = new geographic() ;