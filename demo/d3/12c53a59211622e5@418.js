import define1 from "./1713dcdd0861b14d@1366.js";

function _1(md){return(
md`# D3 Projections

The big list — with projections from [d3-geo](https://github.com/d3/d3-geo), [d3-geo-projection](https://github.com/d3/d3-geo-projection) and [d3-geo-polygon](https://github.com/d3/d3-geo-polygon).

To use in another notebook:
> \`\`\`{javascript}
import {projectionInput} from "@fil/d3-projections"
\`\`\`

or, with your own \`d3\`:
> \`\`\`{javascript}
import {projectionInput} with {d3} from "@fil/d3-projections"
\`\`\`

`
)}

function _projection(projectionInput2,URLSearchParams,html){return(
projectionInput2({
  name: "projection",
  value:
    new URLSearchParams(html`<a href>`.search).get("projection") || "Mercator"
})
)}

function _isea9r(FileAttachment){return(
FileAttachment("isea9r.geojson").json()
)}

function _isea9r_centroids(FileAttachment){return(
FileAttachment("isea9r_centroids.geojson").json()
)}

function _isea9r_1(FileAttachment){return(
FileAttachment("isea9r_1.geojson").json()
)}

function _isea4r_1(FileAttachment){return(
  FileAttachment("isea4r_1.geojson").json()
  )}      
  

function _isea3h(FileAttachment){return(
FileAttachment("isea3h.geojson").json()
)}

function _isea7h_1(FileAttachment){return(
  FileAttachment("isea7h_1.geojson").json()
  )}

function _isea3h_1(FileAttachment){return(
  FileAttachment("isea3h_1.geojson").json()
  )}
  
function _point(FileAttachment){return(
FileAttachment("point.geojson").json()
)}

function _landPoints(FileAttachment){return(
FileAttachment("land.geojson").json()
)}

function _icosahedron(FileAttachment){return(
FileAttachment("icosahedron.geojson").json()
)}

function _icosahedron_pole(FileAttachment){return(
FileAttachment("icosahedron_pole.geojson").json()
)}

function _fuller(FileAttachment){return(
FileAttachment("fuller.geojson").json()
)}

function _fuller_point(FileAttachment){return(
FileAttachment("fuller_point.geojson").json()
)}

function _3(map,projection,isea9r,isea9r_centroids,isea9r_1,isea4r_1,isea3h,isea7h_1,isea3h_1,point,landPoints,icosahedron,icosahedron_pole,fuller,fuller_point){return(
map(projection, { ...projection.options, isea9r, isea9r_centroids, isea9r_1, isea4r_1, isea3h, isea7h_1, isea3h_1,point, landPoints, icosahedron, icosahedron_pole, fuller, fuller_point })
)}

function _4(md){return(
md`---


## Appendix

Edit the *projections* array below to define a new projection.`
)}

function _projections(d3){return(
[
  { name: "Airocean", value: d3.geoAirocean },
  { name: "Airy’s minimum error", value: d3.geoAiry },
  { name: "Aitoff", value: d3.geoAitoff },
  { name: "American polyconic", value: d3.geoPolyconic },
  {
    name: "Armadillo",
    value: d3.geoArmadillo,
    options: { clip: { type: "Sphere" } }
  },
  { name: "August", value: d3.geoAugust },
  { name: "azimuthal equal-area", value: d3.geoAzimuthalEqualArea },
  { name: "azimuthal equidistant", value: d3.geoAzimuthalEquidistant },
  { name: "Baker dinomic", value: d3.geoBaker },
  { name: "Berghaus’ star", value: d3.geoBerghaus,
    options: { clip: { type: "Sphere" } } },
  { name: "Bertin’s 1953", value: d3.geoBertin1953 },
  { name: "Boggs’ eumorphic", value: d3.geoBoggs },
  { name: "Boggs’ eumorphic (interrupted)", value: d3.geoInterruptedBoggs,
    options: { clip: { type: "Sphere" } } },
  { name: "Bonne", value: d3.geoBonne },
  { name: "Bottomley", value: d3.geoBottomley },
  { name: "Bromley", value: d3.geoBromley },
  { name: "Butterfly (gnomonic)", value: d3.geoPolyhedralButterfly },
  { name: "Butterfly (Collignon)", value: d3.geoPolyhedralCollignon },
  { name: "Butterfly (Waterman)", value: d3.geoPolyhedralWaterman },
  { name: "Cahill-Keyes", value: d3.geoCahillKeyes },
  { name: "Collignon", value: d3.geoCollignon },
  // {name: "conic conformal", value: d3.geoConicConformal}, // Not suitable for world maps.
  { name: "conic equal-area", value: d3.geoConicEqualArea },
  { name: "conic equidistant", value: d3.geoConicEquidistant },
  { name: "Craig retroazimuthal", value: d3.geoCraig },
  { name: "Craster parabolic", value: d3.geoCraster },
  { name: "Cox", value: d3.geoCox },
  { name: "cubic", value: d3.geoCubic },
  { name: "cylindrical equal-area", value: d3.geoCylindricalEqualArea },
  {
    name: "cylindrical stereographic",
    value: d3.geoCylindricalStereographic
  },
  { name: "dodecahedral", value: () => d3.geoDodecahedral().parents([-1,0,6,8,1,2,4,5,9,10,11,7]).angle(121).rotate([75,0,-32]) },
  { name: "Eckert I", value: d3.geoEckert1 },
  { name: "Eckert II", value: d3.geoEckert2 },
  { name: "Eckert III", value: d3.geoEckert3 },
  { name: "Eckert IV", value: d3.geoEckert4 },
  { name: "Eckert V", value: d3.geoEckert5 },
  { name: "Eckert VI", value: d3.geoEckert6 },
  { name: "Eisenlohr conformal", value: d3.geoEisenlohr },
  { name: "Equal Earth", value: d3.geoEqualEarth },
  { name: "Equirectangular (plate carrée)", value: d3.geoEquirectangular },
  { name: "Fahey pseudocylindrical", value: d3.geoFahey },
  { name: "flat-polar parabolic", value: d3.geoMtFlatPolarParabolic },
  { name: "flat-polar quartic", value: d3.geoMtFlatPolarQuartic },
  { name: "flat-polar sinusoidal", value: d3.geoMtFlatPolarSinusoidal },
  { name: "Foucaut’s stereographic equivalent", value: d3.geoFoucaut },
  { name: "Foucaut’s sinusoidal", value: d3.geoFoucautSinusoidal },
  { name: "general perspective", value: d3.geoSatellite },
  // {name: "Gilbert’s two-world", value: d3.geoGilbert}, // https://github.com/d3/d3-geo-projection/issues/165
  { name: "Gingery", value: d3.geoGingery,
    options: { clip: { type: "Sphere" } } },
  { name: "Ginzburg V", value: d3.geoGinzburg5 },
  { name: "Ginzburg VI", value: d3.geoGinzburg6 },
  { name: "Ginzburg VIII", value: d3.geoGinzburg8 },
  { name: "Ginzburg IX", value: d3.geoGinzburg9 },
  { name: "Goode’s homolosine", value: d3.geoHomolosine},
  {
    name: "Goode’s homolosine (interrupted)",
    value: d3.geoInterruptedHomolosine,
    options: { clip: { type: "Sphere" } } 
  },
  { name: "gnomonic", value: d3.geoGnomonic },
  { name: "Gringorten square", value: d3.geoGringorten },
  { name: "Gringorten quincuncial", value: d3.geoGringortenQuincuncial },
  { name: "Guyou square", value: d3.geoGuyou },
  { name: "Hammer", value: d3.geoHammer },
  { name: "Hammer retroazimuthal", value: d3.geoHammerRetroazimuthal,
    options: { clip: { type: "Sphere" } } },
  { name: "HEALPix", value: d3.geoHealpix,
    options: { clip: { type: "Sphere" } } },
  { name: "Hill eucyclic", value: d3.geoHill },
  { name: "Hufnagel pseudocylindrical", value: d3.geoHufnagel },
  { name: "icosahedral", value: d3.geoIcosahedral },
  { name: "Imago", value: d3.geoImago },
  { name: "Kavrayskiy VII", value: d3.geoKavrayskiy7 },
  { name: "Lagrange conformal", value: d3.geoLagrange },
  { name: "Larrivée", value: d3.geoLarrivee },
  { name: "Laskowski tri-optimal", value: d3.geoLaskowski },
  // {name: "Littrow retroazimuthal", value: d3.geoLittrow}, // Not suitable for world maps.
  { name: "Loximuthal", value: d3.geoLoximuthal },
  { name: "Mercator", value: d3.geoMercator },
  { name: "Miller cylindrical", value: d3.geoMiller },
  { name: "Mollweide", value: d3.geoMollweide },
  {
    name: "Mollweide (Goode’s interrupted)",
    value: d3.geoInterruptedMollweide,
    options: { clip: { type: "Sphere" } }
  },
  {
    name: "Mollweide (interrupted hemispheres)",
    value: d3.geoInterruptedMollweideHemispheres,
    options: { clip: { type: "Sphere" } }
  },
  { name: "Natural Earth", value: d3.geoNaturalEarth1 },
  { name: "Natural Earth II", value: d3.geoNaturalEarth2 },
  { name: "Nell–Hammer", value: d3.geoNellHammer },
  { name: "Nicolosi globular", value: d3.geoNicolosi },
  { name: "orthographic", value: d3.geoOrthographic },
  { name: "Patterson cylindrical", value: d3.geoPatterson },
  { name: "Peirce quincuncial", value: d3.geoPeirceQuincuncial },
  { name: "rectangular polyconic", value: d3.geoRectangularPolyconic },
  { name: "Robinson", value: d3.geoRobinson },
  { name: "sinusoidal", value: d3.geoSinusoidal },
  { name: "sinusoidal (interrupted)", value: d3.geoInterruptedSinusoidal,
    options: { clip: { type: "Sphere" } } },
  { name: "sinu-Mollweide", value: d3.geoSinuMollweide },
  {
    name: "sinu-Mollweide (interrupted)",
    value: d3.geoInterruptedSinuMollweide,
    options: { clip: { type: "Sphere" } }
  },
  { name: "stereographic", value: d3.geoStereographic },
  { name: "Lee’s tetrahedal", value: d3.geoTetrahedralLee },
  { name: "Times", value: d3.geoTimes },
  { name: "Tobler hyperelliptical", value: d3.geoHyperelliptical },
  { name: "transverse Mercator", value: d3.geoTransverseMercator },
  { name: "Van der Grinten", value: d3.geoVanDerGrinten },
  { name: "Van der Grinten II", value: d3.geoVanDerGrinten2 },
  { name: "Van der Grinten III", value: d3.geoVanDerGrinten3 },
  { name: "Van der Grinten IV", value: d3.geoVanDerGrinten4 },
  { name: "Wagner IV", value: d3.geoWagner4 },
  { name: "Wagner VI", value: d3.geoWagner6 },
  { name: "Wagner VII", value: d3.geoWagner7 },
  {
    name: "Werner",
    value: d3.geoBonne ? () => d3.geoBonne().parallel(90) : null
  },
  { name: "Wiechel", value: d3.geoWiechel },
  { name: "Winkel tripel", value: d3.geoWinkel3 }
].filter(p => p.value)
)}

function _projectionInput(Inputs,projections){return(
function projectionInput({ name, value }) {
  return Inputs.select(
    new Map(
      projections.map(({ name, value, options }) => [
        name,
        // TODO: Select forces us to instantiate every projection… see projectionInput2 below for tests
        // see https://github.com/observablehq/inputs/issues/128#issuecomment-829275682 for a solution
        // TODO: apply it here…
        Object.assign(value(), { options: { name, ...options } })
      ])
    ),
    {
      key: value,
      label: name
    }
  );
}
)}

function _d3(require){return(
require("d3-array@3", "d3-geo@3", "d3-geo-projection@4", "d3-geo-polygon@1.8")
)}

function _9(md){return(
md`---
*These tests are trying to avoid instantiating all the projections while still allowing to use Select(…, value: "Mercator")*`
)}

function _projectionInput2(d3,projections,Inputs){return(
function projectionInput2({ name, value }) {
  const m = d3.index(projections, (d) => d.name);
  return Inputs.select(
    projections.map(({ name }) => name),
    {
      value, // this breaks because of valueof :(
      valueof: (name) =>
        Object.assign(m.get(name).value(), {
          options: { name, ...m.get(name).options }
        }),
      label: name
    }
  );
}
)}

function _p2(projectionInput2){return(
projectionInput2({ name: "projection", value: "Mercator" })
)}

function _12(p2){return(
p2
)}

function _13(p2){return(
p2.options
)}

export default function define(runtime, observer) {
  const main = runtime.module();
  function toString() { return this.url; }
  const fileAttachments = new Map([
    ["isea9r.geojson", {url: new URL("./files/isea9r.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["isea9r_centroids.geojson", {url: new URL("./files/isea9r_centroids.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["isea9r_1.geojson", {url: new URL("./files/isea9r_1.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["isea4r_1.geojson", {url: new URL("./files/isea4r_1.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["isea3h.geojson", {url: new URL("./files/isea3h.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["isea7h_1.geojson", {url: new URL("./files/isea7h_1.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["isea3h_1.geojson", {url: new URL("./files/isea3h_1.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["point.geojson", {url: new URL("./files/point.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["land.geojson", {url: new URL("./files/land.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["icosahedron.geojson", {url: new URL("./files/icosahedron.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["icosahedron_pole.geojson", {url: new URL("./files/icosahedron_pole.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["fuller.geojson", {url: new URL("./files/fuller.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["fuller_point.geojson", {url: new URL("./files/fuller_point.geojson", import.meta.url), mimeType: "application/json", toString}]
  ]);
  main.builtin("FileAttachment", runtime.fileAttachments(name => fileAttachments.get(name)));
  main.variable(observer()).define(["md"], _1);
  main.variable(observer("viewof projection")).define("viewof projection", ["projectionInput2","URLSearchParams","html"], _projection);
  main.variable(observer("projection")).define("projection", ["Generators", "viewof projection"], (G, _) => G.input(_));
  main.variable(observer("isea9r")).define("isea9r", ["FileAttachment"], _isea9r);
  main.variable(observer("isea9r_centroids")).define("isea9r_centroids", ["FileAttachment"], _isea9r_centroids);
  main.variable(observer("isea9r_1")).define("isea9r_1", ["FileAttachment"], _isea9r_1);
  main.variable(observer("isea4r_1")).define("isea4r_1", ["FileAttachment"], _isea4r_1);
  main.variable(observer("isea3h")).define("isea3h", ["FileAttachment"], _isea3h);
  main.variable(observer("isea7h_1")).define("isea7h_1", ["FileAttachment"], _isea7h_1);
  main.variable(observer("isea3h_1")).define("isea3h_1", ["FileAttachment"], _isea3h_1);
  main.variable(observer("point")).define("point", ["FileAttachment"], _point);
  main.variable(observer("landPoints")).define("landPoints", ["FileAttachment"], _landPoints);
  main.variable(observer("icosahedron")).define("icosahedron", ["FileAttachment"], _icosahedron);
  main.variable(observer("icosahedron_pole")).define("icosahedron_pole", ["FileAttachment"], _icosahedron_pole);
  main.variable(observer("fuller")).define("fuller", ["FileAttachment"], _fuller);
  main.variable(observer("fuller_point")).define("fuller_point", ["FileAttachment"], _fuller_point);
  main.variable(observer()).define(["map","projection","isea9r","isea9r_centroids","isea9r_1","isea4r_1","isea3h","isea7h_1","isea3h_1","point","landPoints","icosahedron","icosahedron_pole","fuller","fuller_point"], _3);
  main.variable(observer()).define(["md"], _4);
  main.variable(observer("projections")).define("projections", ["d3"], _projections);
  main.variable(observer("projectionInput")).define("projectionInput", ["Inputs","projections"], _projectionInput);
  const child1 = runtime.module(define1);
  main.import("map", child1);
  main.variable(observer("d3")).define("d3", ["require"], _d3);
  main.variable(observer()).define(["md"], _9);
  main.variable(observer("projectionInput2")).define("projectionInput2", ["d3","projections","Inputs"], _projectionInput2);
  main.variable(observer("viewof p2")).define("viewof p2", ["projectionInput2"], _p2);
  main.variable(observer("p2")).define("p2", ["Generators", "viewof p2"], (G, _) => G.input(_));
  main.variable(observer()).define(["p2"], _12);
  main.variable(observer()).define(["p2"], _13);    
  return main;
}
