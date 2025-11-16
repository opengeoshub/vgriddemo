import define1 from "./e1ede52dc01ead30@252.js";

function _1(md){return(
md`<div style="color: grey; font: 13px/25.5px var(--sans-serif); text-transform: uppercase;"><h1 style="display: none;">World map (canvas)</h1><a href="https://d3js.org/">D3</a> › <a href="/@d3/gallery">Gallery</a></div>

# World map (canvas)

Compare to [SVG](/@d3/world-map-svg).`
)}

function _projection(projectionInput,URLSearchParams,location){return(
projectionInput({value: new URLSearchParams(location.search).get("projection") || "orthographic"})
)}

function _map(DOM,width,height,d3,projection,outline,graticule,land,healpix,healpixZCurve,healpix0,centroids,centroids1)
{
  const context = DOM.context2d(width, height);
  const path = d3.geoPath(projection, context);
  context.save();
  context.beginPath(), path(outline), context.clip(), context.fillStyle = "#fff", context.fillRect(0, 0, width, height);
  // context.beginPath(), path(graticule), context.strokeStyle = "#ccc", context.stroke();
  context.beginPath(), path(land), context.fillStyle = "#c2c8ca", context.fill();
   context.beginPath(), path(healpix), context.strokeStyle = "#ff0000", context.lineWidth = 2, context.stroke();
  // context.beginPath(), path(healpixZCurve), context.strokeStyle = "#0000ff", context.lineWidth = 1.5, context.setLineDash([5, 5]), context.stroke();
  // Draw healpix_0 as polylines with blue color
  if (healpix0 && healpix0.features && Array.isArray(healpix0.features) && healpix0.features.length > 0) {
    context.strokeStyle = "#008000"; // Blue color
    context.lineWidth = 4;
    context.setLineDash([]); // Solid line (no dash)
    
    healpix0.features.forEach(feature => {
      if (feature && feature.geometry) {
        context.beginPath();
        path(feature);
        context.stroke();
      }
    });
  }
  
  // Draw labels from centroids (label only, no point markers)
  if (centroids && centroids.features && Array.isArray(centroids.features)) {
    context.font = "18px sans-serif";
    context.textAlign = "center";
    context.textBaseline = "middle";
    centroids.features.forEach(feature => {
      if (!feature || !feature.geometry || feature.geometry.type !== "Point") return;
      const coords = feature.geometry.coordinates;
      const projected = projection(coords);
      if (!projected) return;
      const [x, y] = projected;
      const label = feature.properties && feature.properties.dggal_healpix ? feature.properties.dggal_healpix : "";
      if (!label) return;
      // white halo for readability
      context.fillStyle = "#008000";
      for (let dx = -1; dx <= 1; dx++) {
        for (let dy = -1; dy <= 1; dy++) {
          if (dx === 0 && dy === 0) continue;
          context.fillText(label, x + dx, y + dy);
        }
      }
      // main text
      context.fillStyle = "#000000";
      context.fillText(label, x, y);
    });
  }
  
  // Draw labels from centroids_1 (red labels)
  if (centroids1 && centroids1.features && Array.isArray(centroids1.features)) {
    context.font = "18px sans-serif";
    context.textAlign = "center";
    context.textBaseline = "middle";
    centroids1.features.forEach(feature => {
      if (!feature || !feature.geometry || feature.geometry.type !== "Point") return;
      const coords = feature.geometry.coordinates;
      const projected = projection(coords);
      if (!projected) return;
      const [x, y] = projected;
      const label = feature.properties && feature.properties.dggal_healpix ? feature.properties.dggal_healpix : "";
      if (!label) return;
      context.fillStyle = "#ff0000";
      context.fillText(label, x, y);
    });
  }
  
    
  // Draw labels with white halo
  // context.font = "14px sans-serif";
  // context.textAlign = "center";
  // context.textBaseline = "middle";
  // healpix.features.forEach(feature => {
  //   const centroid = d3.geoCentroid(feature.geometry);
  //   const [x, y] = projection(centroid);
  //   if (x != null && y != null && feature.properties && feature.properties.dggal_healpix) {
  //     const label = feature.properties.dggal_healpix;
  //     // Draw white halo by drawing text multiple times with offsets
  //     context.fillStyle = "#ffffff";
  //     context.lineWidth = 1;
  //     context.strokeStyle = "#ffffff";
  //     for (let dx = -1; dx <= 1; dx++) {
  //       for (let dy = -1; dy <= 1; dy++) {
  //         if (dx !== 0 || dy !== 0) {
  //           context.strokeText(label, x + dx, y + dy);
  //         }
  //       }
  //     }
  //     // Draw the main text
  //     context.fillStyle = "#000000";
  //     context.fillText(label, x, y);
  //   }
  // });
  
  // Draw healpix_0 labels with red color
  // if (healpix0 && healpix0.features && Array.isArray(healpix0.features) && healpix0.features.length > 0) {
  //   context.font = "20px sans-serif";
  //   context.textAlign = "center";
  //   context.textBaseline = "middle";
  //   healpix0.features.forEach(feature => {
  //     if (feature && feature.geometry && feature.properties) {
  //       const centroid = d3.geoCentroid(feature.geometry);
  //       const [x,y] = projection(centroid);
  //       if (x != null && y != null && feature.properties.dggal_healpix) {
  //         const label = feature.properties.dggal_healpix;
  //         // Draw the text in red without halo
  //         context.fillStyle = "#ff0000";
  //         context.fillText(label, x, y);
  //       }
  //     }
  //   });
  // }
  
  context.restore();
  context.beginPath(), path(outline), context.strokeStyle = "#000", context.stroke();
  return context.canvas;
}


function _height(d3,projection,width,outline)
{
  const [[x0, y0], [x1, y1]] = d3.geoPath(projection.fitWidth(width, outline)).bounds(outline);
  const dy = Math.ceil(y1 - y0), l = Math.min(Math.ceil(x1 - x0), dy);
  projection.scale(projection.scale() * (l - 1) / l).precision(0.2);
  return dy;
}


function _outline(){return(
{type: "Sphere"}
)}

function _graticule(d3){return(
d3.geoGraticule10()
)}

function _land(topojson,world){return(
topojson.feature(world, world.objects.land)
)}

function _world(FileAttachment){return(
FileAttachment("land-50m.json").json()
)}

function _healpix(FileAttachment){return(
FileAttachment("healpix.geojson").json()
)}

function _healpixZCurve(FileAttachment){return(
FileAttachment("healpix_2_z_curve.geojson").json()
)}

function _healpix0(FileAttachment){return(
FileAttachment("healpix_0.geojson").json()
)}

function _centroids(FileAttachment){return(
FileAttachment("centroids.geojson").json()
)}

function _centroids1(FileAttachment){return(
FileAttachment("centroids_1.geojson").json()
)}

function _d3(require){return(
require("d3-geo@3", "d3-geo-projection@4")
)}

export default function define(runtime, observer) {
  const main = runtime.module();
  function toString() { return this.url; }
  const fileAttachments = new Map([
    ["land-50m.json", {url: new URL("./files/efcaaf9f0b260e09b6afeaee6dbc1b91ad45f3328561cd67eb16a1754096c1095f70d284acdc4b004910e89265b60eba2706334e0dc84ded38fd9209083d4cef.json", import.meta.url), mimeType: "application/json", toString}],
    ["healpix.geojson", {url: new URL("./files/healpix.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["healpix_2_z_curve.geojson", {url: new URL("./files/healpix_2_z_curve.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["healpix_0.geojson", {url: new URL("./files/healpix_0.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["centroids.geojson", {url: new URL("./files/centroids.geojson", import.meta.url), mimeType: "application/json", toString}],
    ["centroids_1.geojson", {url: new URL("./files/centroids_1.geojson", import.meta.url), mimeType: "application/json", toString}]
  ]);
  main.builtin("FileAttachment", runtime.fileAttachments(name => fileAttachments.get(name)));
  main.variable(observer()).define(["md"], _1);
  main.variable(observer("viewof projection")).define("viewof projection", ["projectionInput","URLSearchParams","location"], _projection);
  main.variable(observer("projection")).define("projection", ["Generators", "viewof projection"], (G, _) => G.input(_));
  main.variable(observer("map")).define("map", ["DOM","width","height","d3","projection","outline","graticule","land","healpix","healpixZCurve","healpix0","centroids","centroids1"], _map);
  main.variable(observer("height")).define("height", ["d3","projection","width","outline"], _height);
  main.variable(observer("outline")).define("outline", _outline);
  main.variable(observer("graticule")).define("graticule", ["d3"], _graticule);
  main.variable(observer("land")).define("land", ["topojson","world"], _land);
  main.variable(observer("world")).define("world", ["FileAttachment"], _world);
  main.variable(observer("healpix")).define("healpix", ["FileAttachment"], _healpix);
  main.variable(observer("healpixZCurve")).define("healpixZCurve", ["FileAttachment"], _healpixZCurve);
  main.variable(observer("healpix0")).define("healpix0", ["FileAttachment"], _healpix0);
  main.variable(observer("centroids")).define("centroids", ["FileAttachment"], _centroids);
  main.variable(observer("centroids1")).define("centroids1", ["FileAttachment"], _centroids1);
  main.variable(observer("d3")).define("d3", ["require"], _d3);
  const child1 = runtime.module(define1);
  main.import("projectionInput", child1);
  return main;
}
