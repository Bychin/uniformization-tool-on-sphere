const SERVER_URL = "http://127.0.0.1:8080"
const ISOLINE_API = "/api/isoline"
const STATS_API = "/api/stats"
  
const ISOLINE_AREA_RATIO = [0.9, 0.7, 0.5, 0.3, 0.1];


let sphere = new Sphere(CLASSIC_SPHERE, DRAW_MODE.LINE_STRIP, 36, 1, {x:0, y:0, z:0});
let graphic = new Graphic(sphere);

function main() {
    graphic.tick();
}