import define1 from "./cd89ef62a7430a3c@66.js";

function _1(md){return(
md`# HEALPix

A <b>H</b>ierarchical <b>E</b>qual <b>A</b>rea iso<b>L</b>atitude <b>Pix</b>elisation of a 2-sphere.`
)}

function _lobes(html)
{
  const form = html`<form>
  <input name=i type=range min=1 max=12 value=4 step=1 style="width:180px;">
  <output style="font-size:smaller;font-style:oblique;" name=o></output>
</form>`;
  form.i.oninput = () => form.o.value = `${form.value = form.i.valueAsNumber} lobes`;
  form.i.oninput();
  return form;
}


function _projection(d3){return(
d3.geoHealpix()
)}

function _4(map){return(
map
)}

export default function define(runtime, observer) {
  const main = runtime.module();
  main.variable(observer()).define(["md"], _1);
  main.variable(observer("viewof lobes")).define("viewof lobes", ["html"], _lobes);
  main.variable(observer("lobes")).define("lobes", ["Generators", "viewof lobes"], (G, _) => G.input(_));
  main.variable(observer("projection")).define("projection", ["d3"], _projection);
  main.variable(observer()).define(["map"], _4);
  const child1 = runtime.module(define1).derive(["projection"], main);
  main.import("d3", child1);
  main.import("map", child1);
  return main;
}
