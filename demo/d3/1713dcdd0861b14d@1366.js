// https://observablehq.com/@fil/base-map@1366
function _1(md){return(
md`# Base map

Usage:
\`\`\`javascript
import { map } from "@fil/base-map"
map(projection, [opts])
\`\`\`

\`projection\` can be a d3 projection or a string, such as \`"orthographic"\`.

The options are:
\`\`\`javascript
{
  svg: false, // defaults to canvas
  rotate: [0,0], // initial rotation
  inertia: true, // allow mouse rotation with d3-inertia
  clip: [GeoJSON], // for example, clip to Sphere
  background: "white", // background color
}
\`\`\`

\`map()\` returns a graphical object (svg or canvas), with properties that can be read and overloaded: \`projection\`, \`path\` and \`render\`.

**Important:** Don't trust me. I'm still testing this. The API may change at any time without warning. If you need to import this base map, please fork it so your notebooks can be safe.
`
)}

function _EXAMPLE(map,d3){return(
map(d3.geoOrthographic(), {
  svg: true,
  rotate: [0, -90],
  inertia: true
})
)}

function _EXAMPLE2(map,d3){return(
map(d3.geoAirocean(), {
  show_equator: true,
  show_structure: true,
  background: "white"
})
)}

function _4(md){return(
md`------
_Tech zone_`
)}

function _map(d3,width,DOM,mapsvg,mapcanvas,invalidation){return(
function (projection, opts) {
  if (!projection) {
    projection = d3.geoOrthographic();
  } else if (typeof projection == "string") {
    const name = "geo" + projection[0].toUpperCase() + projection.slice(1);
    if (name in d3) projection = d3[name]();
    else throw `Error: "${name}" not found in d3.`;
  }

  opts = { inertia: true, zoom: true, ...opts };

  var c = this;
  if (!c && projection.rotate && opts.rotate) {
    projection.rotate(opts.rotate);
  }

  var path = d3.geoPath(projection);
  var bounds = path.bounds({ type: "Sphere" }),
    ratio = (bounds[1][1] - bounds[0][1]) / (bounds[1][0] - bounds[0][0]),
    height = opts.height || width * Math.sqrt(ratio * 0.6);

  if (opts.fitExtent !== false)
    projection.fitExtent(
      [
        [2, 2],
        [width - 2, height - 2]
      ],
      { type: "Sphere" }
    );

  if (opts.svg) {
    if (c && c.getContext) c = false;
    c = c || DOM.svg(width, height);
    mapsvg(projection, c, opts);
  } else {
    if (c && !c.getContext) c = false;
    if (!c) c = DOM.context2d(width, height).canvas;
    mapcanvas(projection, c, opts);
  }

  c.projection = projection;

  if (projection.invert && opts.inertia && d3.geoInertiaDrag) {
    let sel = d3
      .select(c)
      .style("cursor", "-webkit-grab")
      .style("cursor", "-moz-grab")
      .style("cursor", "grab");
    d3.geoInertiaDrag(
      sel,
      function () {
        c.render();
      } /* allow overloading map.render */,
      c.projection
    );
    invalidation.then(() => sel.on(".drag", null));
  }

  /* not ready
  if (opts.zoom) {
    d3.select(c).call(zoom) 
  }
  */

  return c;
}
)}

function _mapsvg(d3,graticule,land){return(
function(projection, c, opts) {
  const svg = d3.select(c);
  const path = d3.geoPath(projection);

  if (opts.clip) {
    var defs = svg.append("defs");
    defs
      .append("path")
      .datum(opts.clip)
      .attr("id", "clipOpt")
      .attr("d", path);
    defs
      .append("clipPath")
      .attr("id", "clip")
      .append("use")
      .attr("xlink:href", "#clipOpt");
  }

  svg
    .append("path")
    .datum({ type: "Sphere" })
    .attr("stroke", "black")
    .attr("fill", "#fefef2");

  svg
    .append("path")
    .datum(graticule)
    .attr("stroke", "black")
    .attr("stroke-width", 0.2)
    .attr("fill", "none");

  var l = land.features
    ? land
    : {
        type: "FeatureCollection",
        features: land
      };

  svg
    .append("g")
    .selectAll("path")
    .data(l.features)
    .enter()
    .append("path")
    .attr("id", d => d.properties.name)
    .attr("fill", "black");

  function render() {
    svg
      .selectAll("path")
      .attr("d", path)
      .attr("clip-path", "url(#clip)");
  }

  c.path = path;
  (c.render = render)();

  return c;
}
)}

function _mapcanvas(d3,width,graticule,land){return(
function mapcanvas(
  projection,
  c,
  {
    background = "#fff",
    clip = false,
    keepcanvas = false,
    show_structure = false,
    show_sphere = true,
    show_equator = false,
    show_land = true,
    show_graticule = false,
    fill = () => "#c2c8ca",
    isea9r = null,
    isea9r_centroids = null,
    isea9r_1 = null,
    isea4r_1 = null,
    isea3h = null,
    isea7h_1 = null,
    isea3h_1 = null,
    point = null,
    landPoints = null,
    icosahedron = null,
    icosahedron_pole = null,
    fuller = null,
    fuller_point = null,
    ...opts
  } = {}
) {
  const context = c.getContext("2d");

  var path = d3.geoPath(projection, context);

  if (clip) context.beginPath(), path(clip), context.clip();

  function render() {
    if (!keepcanvas) context.clearRect(0, 0, width, c.height);
    if (show_sphere) {
      context.beginPath();
      path({ type: "Sphere" });
      if (background) {
        context.fillStyle = background;
        context.fill();
      }

      if (show_graticule) {
        context.beginPath(),
          path(graticule),
          (context.strokeStyle = "#ccc"),
          context.stroke();
      }
    }
    if (show_land) {
      var l = land.features ? land.features : [land];
      l.forEach((f, i) => {
        context.beginPath(),
          path(f),
          (context.fillStyle = fill(f, i)),
          context.fill();
      });
    }
    if (isea9r && isea9r.features && Array.isArray(isea9r.features)) {
      context.strokeStyle = "#ff0000";
      context.lineWidth = 4
      context.setLineDash([]);
      isea9r.features.forEach(feature => {
        if (feature && feature.geometry) {
          context.beginPath();
          path(feature);
          context.stroke();
        }
      });
    }
    // if (isea3h && isea3h.features && Array.isArray(isea3h.features)) {
    //   context.strokeStyle = "#ff0000";
    //   context.lineWidth = 4;
    //   context.setLineDash([]);
    //   isea3h.features.forEach(feature => {
    //     if (feature && feature.geometry) {
    //       context.beginPath();
    //       path(feature);
    //       context.stroke();
    //     }
    //   });
    // }
    // if (isea7h_1 && isea7h_1.features && Array.isArray(isea7h_1.features)) {
    //   context.strokeStyle = "#0000ff";
    //   context.lineWidth = 2;
    //   context.setLineDash([]);
    //   isea7h_1.features.forEach(feature => {
    //     if (feature && feature.geometry) {
    //       context.beginPath();
    //       path(feature);
    //       context.stroke();
    //     }
    //   });
    // }
    // if (isea3h_1 && isea3h_1.features && Array.isArray(isea3h_1.features)) {
    //   context.strokeStyle = "#0000ff";
    //   context.lineWidth = 2;
    //   context.setLineDash([]);
    //   isea3h_1.features.forEach(feature => {
    //     if (feature && feature.geometry) {
    //       context.beginPath();
    //       path(feature);
    //       context.stroke();
    //     }
    //   });
    // }


    // // if (icosahedron && icosahedron.features && Array.isArray(icosahedron.features)) {
    //   context.strokeStyle = "#ff0000";
    //   context.lineWidth = 4;
    //   context.setLineDash([]);
    //   icosahedron.features.forEach(feature => {
    //     if (feature && feature.geometry) {
    //       context.beginPath();
    //       path(feature);
    //       context.stroke();
    //     }
    //   });
    // }
    // if (icosahedron_pole && icosahedron_pole.features && Array.isArray(icosahedron_pole.features)) {
    //   context.strokeStyle = "#ff0000";
    //   context.lineWidth = 2;
    //   context.setLineDash([]);
    //   icosahedron_pole.features.forEach(feature => {
    //     if (feature && feature.geometry) {
    //       context.beginPath();
    //       path(feature);
    //       context.stroke();
    //     }
    //   });
    // }
    // if (fuller && fuller.features && Array.isArray(fuller.features)) {
    //   context.strokeStyle = "#ff0000";
    //   context.lineWidth = 4;
    //   context.setLineDash([]);   
    //   fuller.features.forEach(feature => {
    //     if (feature && feature.geometry) {
    //       context.beginPath();
    //       path(feature);
    //       context.stroke();
    //     }
    //   });
    // }
    // if (isea9r_1 && isea9r_1.features && Array.isArray(isea9r_1.features)) {
    //   context.strokeStyle = "#0000ff";
    //   context.lineWidth = 2;
    //   context.setLineDash([]);
    //   isea9r_1.features.forEach(feature => {
    //     if (feature && feature.geometry) {
    //       context.beginPath();
    //       path(feature);
    //       context.stroke();
    //     }
    //   });
    // }
    if (isea4r_1 && isea4r_1.features && Array.isArray(isea4r_1.features)) {
      context.strokeStyle = "#0000ff";
      context.lineWidth = 2;
      context.setLineDash([]);
      isea4r_1.features.forEach(feature => {
        if (feature && feature.geometry) {
          context.beginPath();
          path(feature);
          context.stroke();
        }
      });
    }

    // if (point && point.features && Array.isArray(point.features)) {
    //   point.features.forEach(feature => {
    //     if (feature && feature.geometry && feature.geometry.type === "Point") {
    //       const [x, y] = projection(feature.geometry.coordinates);
    //       if (x !== undefined && y !== undefined && isFinite(x) && isFinite(y)) {
    //         // Draw circle marker
    //         context.beginPath();
    //         context.arc(x, y, 10, 0, 2 * Math.PI);
    //         context.fillStyle = "#0000ff";
    //         context.fill();
    //         context.strokeStyle = "#ffffff";
    //         context.lineWidth = 2;
    //         context.stroke();
            
    //         // Draw hollow text label positioned top-right with spacing
    //         if (feature.properties && feature.properties.label) {
    //           const spacing = 13; // Spacing between marker and label
    //           context.font = "18px sans-serif";
    //           context.textAlign = "left";
    //           context.textBaseline = "top";
    //           context.strokeStyle = "#000000";
    //           context.lineWidth = 1;
    //           context.strokeText(feature.properties.label, x + spacing, y - spacing*2);
    //         }
    //       }
    //     }
    //   });
    // }
    // if (landPoints && landPoints.features && Array.isArray(landPoints.features)) {
    //   landPoints.features.forEach(feature => {
    //     if (feature && feature.geometry && feature.geometry.type === "Point") {
    //       const [x, y] = projection(feature.geometry.coordinates);
    //       if (x !== undefined && y !== undefined && isFinite(x) && isFinite(y)) {
    //         // Draw circle marker
    //         context.beginPath();
    //         context.arc(x, y, 10, 0, 2 * Math.PI);
    //         context.fillStyle = "green";
    //         context.fill();
    //         context.strokeStyle = "#ffffff";
    //         context.lineWidth = 2;
    //         context.stroke();
            
    //         // Draw hollow text label positioned top-right with spacing
    //         if (feature.properties && feature.properties.label) {
    //           const spacing = 13; // Spacing between marker and label
    //           context.font = "18px sans-serif";
    //           context.textAlign = "left";
    //           context.textBaseline = "top";
    //           context.strokeStyle = "#000000";
    //           context.lineWidth = 1;
    //           context.strokeText(feature.properties.label, x + spacing, y - spacing);
    //         }
    //       }
    //     }
    //   });
    // }
    // if (fuller_point && fuller_point.features && Array.isArray(fuller_point.features)) {
    //   fuller_point.features.forEach(feature => {
    //     if (feature && feature.geometry && feature.geometry.type === "Point") {
    //       const [x, y] = projection(feature.geometry.coordinates);
    //       if (x !== undefined && y !== undefined && isFinite(x) && isFinite(y)) {
    //         // Draw circle marker (same style as landPoints but no label)
    //         context.beginPath();
    //         context.arc(x, y, 10, 0, 2 * Math.PI);
    //         context.fillStyle = "blue";
    //         context.fill();
    //         context.strokeStyle = "#ffffff";
    //         context.lineWidth = 2;
    //         context.stroke();
    //       }
    //     }
    //   });
    // }
    // if (isea9r_centroids && isea9r_centroids.features && Array.isArray(isea9r_centroids.features)) {
    //   context.fillStyle = "#000000";
    //   context.font = "10px sans-serif";
    //   context.textAlign = "center";
    //   context.textBaseline = "middle";
    //   isea9r_centroids.features.forEach(feature => {
    //     if (feature && feature.geometry && feature.geometry.type === "Point" && feature.properties && feature.properties.zoneID) {
    //       const [x, y] = projection(feature.geometry.coordinates);
    //       if (x !== undefined && y !== undefined && isFinite(x) && isFinite(y)) {
    //         context.fillText(feature.properties.zoneID, x, y);
    //       }
    //     }
    //   });
    // }
    if (show_equator)
      context.beginPath(),
        path(
          d3
            .geoGraticule()
            .step([0, 100])
            .extent([
              [-179.99, -25],
              [179.99, 25]
            ])()
        ),
        (context.strokeStyle = "brown"),
        context.stroke();
    if (show_sphere || show_structure || !show_land)
      context.beginPath(),
        path({ type: "Sphere" }),
        (context.strokeStyle = "#000"),
        context.stroke();

    // Polyhedral projections expose their structure as projection.tree()
    // To draw them we need to cancel the rotate
    if (show_structure && projection.tree) {
      var rotate = projection.rotate();
      projection.rotate([0, 0, 0]);

      // run the tree of faces to get all sites and folds
      var sites = [],
        folds = [],
        i = 0;
      function recurse(face) {
        var site = d3.geoCentroid({
          type: "MultiPoint",
          coordinates: face.face
        });
        site.id = face.id || i++;
        sites.push(site);
        if (face.children) {
          face.children.forEach(function (child) {
            folds.push({
              type: "LineString",
              coordinates: child.shared.map((e) =>
                d3.geoInterpolate(e, face.centroid)(1e-5)
              )
            });
            recurse(child);
          });
        }
      }
      recurse(projection.tree());

      // sites & numbers
      context.beginPath(),
        path.pointRadius(10)({ type: "MultiPoint", coordinates: sites }),
        (context.fillStyle = "white"),
        (context.strokeStyle = "black"),
        context.fill(),
        context.stroke();
      sites.forEach((site) => {
        (context.textAlign = "center"),
          (context.fillStyle = "black"),
          (context.font = "16px Georgia"),
          (context.textBaseline = "middle"),
          context.fillText(
            site.id,
            projection(site)[0],
            projection(site)[1] - 1
          );
      });

      // folding lines
      folds.forEach((fold) => {
        context.beginPath(),
          (context.lineWidth = 0.5),
          context.setLineDash([3, 4]),
          (context.strokeStyle = "#888"),
          path(fold),
          context.stroke(),
          context.setLineDash([]);
      });

      // restore the projection’s rotate
      projection.rotate(rotate);
    }

    if (false)
      d3.select(context.canvas).on("mousemove", function (event) {
        var gni = projection.invert(d3.pointer(event, this));
        context.beginPath(),
          path.pointRadius(2)({ type: "Point", coordinates: gni }),
          (context.fillStyle = "green"),
          (context.strokeStyle = "black"),
          context.fill(),
          context.stroke();
      });
  }

  context.canvas.path = path;
  (context.canvas.render = render)();

  return context.canvas;
}
)}

function _d3(require){return(
require("d3-selection@3", "d3-geo@3", "d3-fetch@3", "d3-geo-projection@4", "d3-geo-polygon@1.8", "d3-inertia@0.4")
)}

function _land(d3){return(
d3.json(
  "https://unpkg.com/visionscarto-world-atlas@0.0.6/world/110m_land.geojson"
)
)}

function _graticule(d3){return(
d3.geoGraticule()()
)}

function _topojson(){return(
null
)}

export default function define(runtime, observer) {
  const main = runtime.module();
  main.variable(observer()).define(["md"], _1);
  main.variable(observer("EXAMPLE")).define("EXAMPLE", ["map","d3"], _EXAMPLE);
  main.variable(observer("EXAMPLE2")).define("EXAMPLE2", ["map","d3"], _EXAMPLE2);
  main.variable(observer()).define(["md"], _4);
  main.variable(observer("map")).define("map", ["d3","width","DOM","mapsvg","mapcanvas","invalidation"], _map);
  main.variable(observer("mapsvg")).define("mapsvg", ["d3","graticule","land"], _mapsvg);
  main.variable(observer("mapcanvas")).define("mapcanvas", ["d3","width","graticule","land"], _mapcanvas);
  main.variable(observer("d3")).define("d3", ["require"], _d3);
  main.variable(observer("land")).define("land", ["d3"], _land);
  main.variable(observer("graticule")).define("graticule", ["d3"], _graticule);
  main.variable(observer("topojson")).define("topojson", _topojson);
  return main;
}
