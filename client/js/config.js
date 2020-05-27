const SERVER_URL = "http://127.0.0.1:8080"
const ISOLINE_API = "/api/isoline"
const STATS_API = "/api/stats"

const SPIRAL_SPHERE = 'spiral';
const CLASSIC_SPHERE = 'classic';

const DRAW_MODE = {
  POINTS: 'points',
  LINE_STRIP: 'line_strip',
  LINE_LOOP: 'line_loop',
  TRIANGLE_STRIP: 'triangle_strip'
}

const FPS_LIMIT = 30;

const ISOLINE_AREA_RATIO = [0.9, 0.7, 0.5, 0.3, 0.1];
const MAX_POINTS_IN_ISOLINE = 64;

const COLORS = {
  BLACK: [0., 0., 0., 1.],
  RED: [1., 0., 0., 1.],
  LIGHT_RED: [0.5, 0., 0., 0.4],
  BLUE: [0., 0., 0.9, 0.8],
  LIGHT_BLUE: [0., 208/256., 255/256., 1.],
}

const LIGHT_COLOR = [0.8, 0.8, 0.8];
const LIGHT_POS = [5.0, 0.0, 0.5];
const AMBIENT_LIGHT = [0.2, 0.2, 0.2];

const DEFAULT_SPHERE_ANGLE_DEG = 0;

const DEFAULT_CAMERA_POS = {
  x: 5,
  y: 0,
  z: 0
}

const DEFAULT_CAMERA_LOOK_AT = {
  x: 0,
  y: 0,
  z: 0
}

const DEFAULT_CAMERA_UP = {
  x: 0,
  y: 0,
  z: 1
}
