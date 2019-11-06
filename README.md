# geographic

Implementation of geographic conversions.

This library implements a subset of the [Matlab toolbox](https://www.mathworks.com/matlabcentral/fileexchange/50605-geographiclib) for [GeographicLib](https://geographiclib.sourceforge.io/). The functions are not vectorized as in Matlab implementation.

# Installation

``` 
$ npm install geographic
```

# Usage

In [node](https://nodejs.org/en/), do
``` javascript
var geographic = require("geographic");
```

# Examples

## Prepare for use

On `require` the object is loaded with the default ellipsoid (WGS84) and the projection origin on 
* lat = 0 °
* lon = 0 °
* h = 0 m

``` javascript
var geo = require("geographic");
var lat; var lon;
```

## Ellipsoid

The ellipsoid property store the radius, flattening and eccentricity.

### defaultellipsoid method

This method set the ellipsoid property to WGS84.
* radius = 6378137 m
* flattening = 1 / 298.257223563
* eccentricity => correspondent to flattening

``` javascript
geo.defaultellipsoid();
```

### setellipsoid

Set the ellipsoid to desired value

``` javascript
geo.setellipsoid(radius,flattening);
```

### getellipsoid

Return the current ellipsoid

``` javascript
var elip = geo.getellipsoid();
// elip <= {ellipsoid.radius, ellipsoid.flattening, ellipsoid.eccentricity}
```
## Projection origin

The ```origin ``` property store the projection origin.

### setorigin

Set the origin position

``` javascript
geo.setorigin(lat,lon,heigth);
```

### getorigin

Return the current origin

``` javascript
var ori = geo.getorigin();
// ori <= {origin.lat,origin.lon,origin.h}
```

## Conversions

When a conversion method is called it use the ellipsoid and origin already stored on the respective property, therefore this properties must be set before the conversion methods if the user don't intend to use the default values.

### Transverse Mercator inverse

This method convert the position from [transverse mercator](https://en.wikipedia.org/wiki/Transverse_Mercator_projection) projection to geographic position.

```javascript
var pos = geo.tranmerc_inv(5534,2065);
console.log( "lat: " + pos.lat + "°  lon: " + pos.lon + "°" );
//lat: -29.993820151321263°  lon: -50.077351253287304°
```

### Transverse Mercator forward

This method convert the position from geographic position projection to [transverse mercator](https://en.wikipedia.org/wiki/Transverse_Mercator_projection).

```javascript
var pos = geo.tranmerc_fwd(-29.9938,-50.0774);
console.log("x: " + pos.x + "m  y: " + pos.y + "m" );
//x: 5529.297435284063m  y: 2067.236174341291m
```

### UTM/UPS inverse

This method convert the position from [UTM](https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system)/[UPS](https://en.wikipedia.org/wiki/Universal_polar_stereographic_coordinate_system) projection to geographic position.

```javascript
//geo.utmups_inv(UTM.x[m],UTM.y,UTM.zone,isnorth)
var pos = geo.utmups_inv(583447.0615085221,6679518.672292533,22,false);
console.log("lat: " + pos.lat + "°  lon: " + pos.lon + "°" );
//lat: -59.89491145601121°  lon: -48.6432811522374°
```

### UTM/UPS forward

This method convert the position from geographic position projection to [UTM](https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system)/[UPS](https://en.wikipedia.org/wiki/Universal_polar_stereographic_coordinate_system).

```javascript
var pos = geo.utmups_fwd(-30.012461,-50.134703);
console.log("x: " + pos.x + "   y: " + pos.y + "   zone: " + pos.zone + "  isnorth: " + pos.isnorth );
//x: 583447.0615085221   y: 6679518.672292533   zone: 22  isnorth: false
```

### Local cartesian inverse

This method convert the position from [local cartesian](https://en.wikipedia.org/wiki/Local_tangent_plane_coordinates) projection to geographic position.

```javascript
//geo.loccart_inv(x[m],y[m],z[m]);
//x => west to east y => south to north
pos = geo.loccart_inv(5534,2065,0);
console.log("lat: " + pos.lat + "°  lon: " + pos.lon + "°°  h: " + pos.h + "m" );
//lat: -29.993820151982113°  lon: -50.077351263450474°  h: 2.73447461843128m
```

### Local cartesian forward

This method convert the position from geographic position to [local cartesian](https://en.wikipedia.org/wiki/Local_tangent_plane_coordinates) projection.

```javascript
pos = geo.loccart_fwd(-29.993820151982113,-50.077351263450474,0);
console.log("x: " + pos.x + "m  y: " + pos.y + "m  z: " + pos.z + "m");
//x: 5533.997629416378m  y: 2064.99911094806m  z: -2.7344734460590416m
```
