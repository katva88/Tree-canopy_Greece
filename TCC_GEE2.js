var geometry = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.MultiPoint(),
    geometry2 = 
    /* color: #d63000 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[22.729630462829043, 40.78728324180753],
          [22.729630462829043, 40.59151819060351],
          [23.140244476500918, 40.59151819060351],
          [23.140244476500918, 40.78728324180753]]], null, false),
    forest_h = ee.Image("projects/ee-avatitsi/assets/Forest_height_2019_NAFR_gr_laea"),
    maes = ee.Image("projects/ee-avatitsi/assets/life_maes_reclass2"),
    geometry3 = 
    /* color: #d63000 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[22.56855837991539, 38.28205157052546],
          [22.56855837991539, 37.59321695598985],
          [23.84022586038414, 37.59321695598985],
          [23.84022586038414, 38.28205157052546]]], null, false),

var tcc_points = ee.FeatureCollection("projects/ee-avatitsi/assets/TCC_ALL_etrs_CORINE_LIFEMAP_wgs84"),
    table1 = ee.FeatureCollection("projects/ee-avatitsi/assets/aktogrammh_wgs84"),
    grid = ee.FeatureCollection("projects/ee-avatitsi/assets/grid_greece_wgs84"),
    table = ee.FeatureCollection("projects/ee-avatitsi/assets/aktogrammh_buffer_wgs84"),
    S2_median = ee.Image("projects/ee-avatitsi/assets/S2_median"),
    S2_median_test = ee.Image("projects/ee-avatitsi/assets/S2_median_test"),
    tcc_points_30x30 = ee.FeatureCollection("projects/ee-avatitsi/assets/TCC_ALL_etrs_CORINE_LIFEMAP_30x30_wgs84"),
    S2_VIs = ee.Image("projects/ee-avatitsi/assets/VIs_composite"),
    geometry = 
    /* color: #d63000 */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[19.625306594900298, 39.84340758368463],
          [19.625306594900298, 39.602581327695354],
          [20.009828079275298, 39.602581327695354],
          [20.009828079275298, 39.84340758368463]]], null, false),
    PV_AOI2 = ee.FeatureCollection("projects/ee-avatitsi/assets/PV_AOI2");
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////shapefile//////////////////////////////////////////////////////////////////////////////
var shp = ee.FeatureCollection(table);
Map.addLayer(shp, {},'aktogrammi',false);

var points = ee.FeatureCollection(tcc_points);
Map.addLayer(points, {},'tcc_points',true)
var featID= 52742008;//just check for burned areas
var selectPoint = points.filter(ee.Filter.eq('id', featID)).first();
print(selectPoint,'point excluded');
var pointGeo = selectPoint.geometry();
//Map.centerObject(pointGeo, 15);

var points30x30 = ee.FeatureCollection(tcc_points_30x30);
Map.addLayer(points, {},'tcc_points30x30',false)

// var grid= ee.FeatureCollection(grid);
// Map.addLayer(grid, {},'grid',false)
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////define dates/////////////////////////////////
var start20 = ee.Date('2020-07-01');
var end20= ee.Date('2020-09-30');
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////cloud mask///////////////////////////////////
function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask).divide(10000);
}
////////////////////////////////////////////////////////////////////////////////
/////////////////////Sentinel-2//collections 2020/////////////////////////////
var S2_2020 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
                .filterBounds(table)        
                .filterDate(start20,end20)
                .filterMetadata('CLOUDY_PIXEL_PERCENTAGE','less_than',10)
                .map(maskS2clouds)
                .median();

var S2_2020_median=S2_2020.select('B[2-8]','B8A','B11','B12').clip(table);
                  
print(S2_2020_median)

var vis = {min: 0.02, max: 0.3, bands: ['B8', 'B4', 'B3']};
Map.addLayer(S2_2020_median, vis, 's2',true)

// Export.image.toAsset({
//   image:S2_2020_median,
//   description:'S2_composite_2020',
//   assetId: 'projects/ee-avatitsi/assets/S2_2020_median',
//   region:table,
//   scale:20,
//   maxPixels:1e13
// });


// /////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////Vegetation indices//////////////////////////////////////////
//NDVI=(NIR-RED)/(NIR+RED)
//EVI=2.5*(NIR-RED)/(1+NIR+6*RED-7.5*BLUE)
//NDWI=NIR-SWIR/NIR+SWIR
//IPVI=NIR/(NIR+RED)
//RVI=NIR/RED
//DVI=(NIR-RED)
//AVI=[NIR*(1-RED)*(NIR-RED)]1/3
//ΒΙ=(SWIR+RED)-(NIR+BLUE)/(SWIR+RED)+(NIR+BLUE)
//SI=[(1-BLUE)*(1-GREEN)*(1-RED)]1/3
//RENDVI=(NIR-RedEdge)/(NIR+Red Edge)
//S2_2020_median

var NDVI = S2_2020_median.expression('(NIR-RED)/ (NIR+RED)', {
          'NIR': S2_2020_median.select('B8'),
          'RED': S2_2020_median.select('B4')
});

var EVI = S2_2020_median.expression('2.5*(NIR-RED)/(1+NIR+6*RED-7.5*BLUE)',{
          'NIR': S2_2020_median.select('B8'),
          'RED': S2_2020_median.select('B4'),
          'BLUE':S2_2020_median.select('B2')
});

var NDWI = S2_2020_median.expression('(GREEN-NIR) / (GREEN+NIR)', {
          'NIR': S2_2020_median.select('B8'),
          'GREEN': S2_2020_median.select('B3')
});

var IPVI = S2_2020_median.expression('NIR/ (NIR+RED)', {
          'NIR': S2_2020_median.select('B8'),
          'RED': S2_2020_median.select('B4')
});

var RVI = S2_2020_median.expression('NIR/RED', {
          'NIR': S2_2020_median.select('B8'),
          'RED': S2_2020_median.select('B4')
});

var DVI = S2_2020_median.expression('NIR-RED', {
          'NIR': S2_2020_median.select('B8'),
          'RED': S2_2020_median.select('B4')
});

var AVI = S2_2020_median.expression('pow(NIR*(1-RED)*(NIR-RED),1/3)', {
          'NIR': S2_2020_median.select('B8'),
          'RED': S2_2020_median.select('B4')
});

var BI = S2_2020_median.expression('(SWIR+RED)-(NIR+BLUE)/(SWIR+RED)+(NIR+BLUE)', {
          'NIR': S2_2020_median.select('B8'),
          'RED': S2_2020_median.select('B4'),
          'BLUE':S2_2020_median.select('B2'),
          'SWIR':S2_2020_median.select('B11')
});

var RENDVI= S2_2020_median.expression('(NIR-RedEdge)/(NIR+RedEdge)',{
            'NIR':S2_2020_median.select('B8'),
            'RedEdge':S2_2020_median.select('B6')
});

var SI = S2_2020_median.expression('pow((1-BLUE)*(1-GREEN)*(1-RED),1/3)', {
          'BLUE':S2_2020_median.select('B2'),
          'GREEN':S2_2020_median.select('B3'),
          'RED': S2_2020_median.select('B4')
});

var VIs = ee.Image([NDVI.rename('NDVI'),EVI.rename('EVI').toFloat(),NDWI.rename('NDWI'),IPVI.rename('IPVI'),RVI.rename('RVI'),DVI.rename('DVI').toFloat(),AVI.rename('AVI'),BI.rename('BI'),SI.rename('SI'),RENDVI.rename('RENDVI')])
print(VIs)
Map.addLayer(VIs.select('NDVI'), {min: -1, max: 1}, 'NDVI',false);

// Export.image.toAsset({
//   image:VIs,
//   description:'VIs_composite2_2020',
//   assetId: 'projects/ee-avatitsi/assets/VIs_composite2_2020',
//   region:table,
//   scale:20,
//   maxPixels:1e13
// });

var lai= VIs.expression(('3.618 * EVI- 0.0118'), {
        'EVI':VIs.select('EVI')
});

print(lai, 'lai')
Map.addLayer(lai, {min: 0, max: 3}, 'lai',false);
// ////////////////////////////////////////////////////////////////////////////////////////////////////
var srtm=ee.Image("USGS/SRTMGL1_003")
// Calculate slope. Units are degrees, range is [0,90).
var slope = ee.Terrain.slope(srtm);

// Calculate aspect. Units are degrees where 0=N, 90=E, 180=S, 270=W.
var aspect = ee.Terrain.aspect(srtm);

var elevation_slope_aspect=srtm.toFloat().addBands(slope).addBands(aspect)

var elevation_slope_aspect_gr=elevation_slope_aspect.clip(table)
print("srtm",elevation_slope_aspect_gr)
Map.addLayer(elevation_slope_aspect_gr,{min:0,max:3000},'elev-sl-asp')


// Export.image.toAsset({
//   image:elevation_slope_aspect_gr,
//   description:'elevation_slope_aspect_gr',
//   assetId: 'projects/ee-avatitsi/assets/elevation_slope_aspect',
//   region:table,
//   scale:20,
//   maxPixels:1e13
// });
/////////////////////////////////////////////////////////////////////////////////////////////////////
// var targetProjection = S2_2020_median.projection();
// maes = maes.reproject({
//     crs: targetProjection
// }).rename('MAES_Name').clip(geometry2);

// forest_h = forest_h.reproject({
//     crs: targetProjection
// }).rename('forest_h').clip(geometry2);


var stack_all= S2_2020_median.toFloat()
.addBands(VIs.toFloat())
.addBands(elevation_slope_aspect_gr.toFloat())
.addBands(lai.rename('lai').toFloat())
.addBands(maes.rename('MAES_Name').toFloat())
.addBands(forest_h.rename('forest_h').toFloat())
.clip(table)
//.clip(geometry3);

Map.addLayer(stack_all,{},'stack_test')

print(stack_all)

Export.image.toDrive({
  image:stack_all.clip(subsetD),
  description:'stack_all_D_laea_2020',
  folder:'TCC',
  region:subsetD,
  scale:20,
  crs:'EPSG:3035',
  maxPixels:1e13
});
// Export.image.toAsset({
//   image:stack_all,
//   description:'stack_all_clip',
//   assetId: 'projects/ee-avatitsi/assets/stack_all_clip',
//   region:geometry3,
//   scale:20,
//   maxPixels:1e13
// });







